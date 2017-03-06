#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include <igl/sortrows.h>
#include <igl/slice.h>
#include <Eigen/Dense>

#include <glm/glm.hpp>

#include "Types.h"
#include "Laplacian.h"

namespace DGP
{
    /* return bundary of mesh 
     *
     * Args:
     *      F: Faces, [F, 3] MatrixXi
     * Returns:
     *      B: Boundarys, [B, 2] MatrixXi
     */
    BMat boundary(const FMat& F)
    {
        int f_num = F.innerSize();

        /* halfedge info, (sort(u, v), (u, v)) for (u, v) \in Edges */
        BMat HalfEdges(f_num * 3, 4);
		for (int i = 0; i < f_num; i++)
		{
			glm::ivec3 v_id = glm::ivec3(F(i, 0), F(i, 1), F(i, 2));
			for (int j = 0; j < 3; j++)
			{
				glm::ivec2 h_id = glm::ivec2(v_id[j], v_id[(j + 1) % 3]);
				HalfEdges(i * 3 + j, 0) = (h_id[0] < h_id[1]) ? h_id[0] : h_id[1];
				HalfEdges(i * 3 + j, 1) = (h_id[0] < h_id[1]) ? h_id[1] : h_id[0];
                HalfEdges(i * 3 + j, 2) = h_id[0];
                HalfEdges(i * 3 + j, 3) = h_id[1];
			}
		}

        /* sort and get dual halfedges */
        BMat sortedHalf, I;
        igl::sortrows(HalfEdges, true, sortedHalf, I);

        /* label dual halfedges */
        iVec dual_label = iVec::Ones(f_num * 3);
        
        int p = 1;
        while(p < dual_label.size())
        {
            /* if two halfedges are dual, label, then p += 2 */
            if (sortedHalf(p-1, 0) == sortedHalf(p, 0) && sortedHalf(p-1, 1) == sortedHalf(p, 1))
            {
                dual_label(p-1) = 0;
                dual_label(p) = 0;
                p += 2;
            }
            else
            { /* else jump through */
                p ++;
            }
        }
        
        BMat bound = BMat(dual_label.sum(), 2);
        int count = 0;
        for(int i = 0; i < dual_label.size(); i++)
            if (dual_label(i) == 1)
            {
                bound(count, 0) = sortedHalf(i, 2); 
                bound(count, 1) = sortedHalf(i, 3); 
				count++;
            }

        return bound;
    }

	/* build Area matrix of a disk like 2d manifold 
	 *
	 *  Args:
	 *	    V: matrix of vertices positions, V * 3 double
	 *      F: matrix of triangle faces,     F * 3 int
     *
     *  Returns:
     *
     *      A: matrix of area computing, V * V double
     *         usage:  area = -i/4 * (z^*T A z )
	 */
	SpMat buildAreaMatrix(const VMat& V, const FMat& F)
	{
		int v_num = V.innerSize();

		BMat B = boundary(F);
		int b_num = B.innerSize();

		/* build Area matrix */
		std::vector<T> area_coeff;
		area_coeff.clear();
		area_coeff.reserve(b_num * 2);

		for (int i = 0; i < b_num; i++)
		{
			area_coeff.push_back(T(B(i, 0), B(i, 1), 1));
			area_coeff.push_back(T(B(i, 1), B(i, 0), -1));
		}

		SpMat A(v_num, v_num);
		A.setFromTriplets(area_coeff.begin(), area_coeff.end());
		return A;
	}

    /* build comfromal transform energy from disk topology manifold to complex plane
     *  E_C = E_D - A
     *  E_D = 1/2 <\laplace z, z> =>  E_D = (Z^*T) d*d Z 
     *  A = -i/2 * \int d\bar{z} \up dz = - i/4 \sum_{ij}(\bar{z}_i z_j - \bar{z}_j z_i)
     *
	 *  Args:
	 *	    V: matrix of vertices positions, V * 3 double
	 *      F: matrix of triangle faces,     F * 3 int
     *
     *  Returns:
     *      E_c: conformal transform energy, V * V sparse complex
     */
    SpMatC buildParameterizationEnergy(const VMat& V, const FMat& F)
    {
		SpMatC E = (Laplacian(V, F)).cast<complex>() + complex(0, 0.5) * buildAreaMatrix(V, F).cast<complex>();
        return E;
    }

	/* return the redidual between y and its principle eigen direction */
	double residual(const SpMatC& A, const cVec& y)
	{
		cVec r = A*y - (y.dot(A * y)) * y;
		return r.norm();
	}

    /* return eigen vector correspond to smallest eigen value of Matrix  
     * 
     *  Args:
     *      L: matrix
     *
     *  Returns:
     *      v: eigen vector with smallest eigen value
     */
    cVec smallestEigen(const SpMatC& A, const cVec& x, const double e)
    {
        Eigen::SimplicialLDLT<SpMatC> solver;
        solver.compute(A);

        const cVec one = cVec::Ones(x.size());
        cVec y = x - x.mean() * one;

        while ( residual(A, y) > e)
        {
            y = solver.solve(y);
            y -= y.mean() * one;
            y /= y.norm();
        }

        return y;
    }


    /* Compute conformal parameterization for a disk like 2d manifold
     *
	 *  Args:
	 *	    V: matrix of vertices positions, V * 3 double
	 *      F: matrix of triangle faces,     F * 3 int
     *      e: tolarance error when getting smallest eigen vector
     *
     *  Returns:
     *      T: matrix of vertices texture cooridnates, V * 2 double
     */
    TMat conformalParameterization(const VMat& V, const FMat& F, const double e)
    {
        /* get conformal transform energy */
        SpMatC E = buildParameterizationEnergy(V, F);

        /* assign boudary vertices to a disk */
        BMat B = boundary(F);
        //if (!B.innerSize() || B(0,0) != B(B.innerSize()-1, 1))
		if (!B.innerSize())
        {
            printf("Not a disk like mesh!\n");
            return TMat();
        }

        cVec x = cVec::Zero(V.innerSize());
        for(int i = 0; i < B.innerSize(); i++)
        {
            double theta = 2 * M_PI * i / B.innerSize();
            x(B(i, 0)) = complex(cos(theta), sin(theta));
        }

        /* get smallest eigen vector */
        x = smallestEigen(E, x, e);

        /* return texture coordinate */
        TMat T = TMat::Zero(V.innerSize(), 2);
        T.col(0) = x.real();
        T.col(1) = x.imag();

        return T;
    }

	TMat leastSquareConformalParameterization(const VMat& V, const FMat& F, const iVec& b, const TMat& bc)
	{
		int v_num = V.innerSize();
		int p_num = b.size();
		int f_num = v_num - p_num;

		/* get conformal transform energy */
		SpMatC E = buildParameterizationEnergy(V, F);
		
		/* get permutation matrix to make all fixed vertices labeled last */
		SpMatC P(v_num, v_num);
		P.setIdentity();
		
		iVec b_I, b_sort;
		igl::sortrows(b, true, b_sort, b_I);

		cVec bc_complex = bc.col(0).cast<complex>() + complex(0, 1) * bc.col(1).cast<complex>();
		cVec bc_sort;

		igl::slice(bc_complex, b_I, bc_sort);

		for (int i = 0; i < p_num; i++)
		{
			int src = b_sort(p_num - i - 1);
			int targ = v_num - i - 1;
			SpVecC tmp = P.col(targ);
			P.col(targ) = P.col(src);
			P.col(src) = tmp;
		}

		/* get least square equation */
		SpMatC A = P.transpose() * E * P;
		cVec y = A.rightCols(p_num) * bc_sort;
		y = y.topRows(f_num);
		A = A.topLeftCorner(f_num, f_num);

		Eigen::SimplicialLDLT<SpMatC> solver;
		solver.compute(A);

		cVec x_tmp = solver.solve(y);

		cVec x(v_num);
		x << x_tmp, bc_sort;

		x = P * x;

		TMat texture(v_num, 2);
		texture.col(0) = x.real();
		texture.col(1) = x.imag();

		return texture;
	}

}
