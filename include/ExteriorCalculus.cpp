#include <vector>

#include "Types.h"
#include "ExteriorCalculus.h"

#include "igl/doublearea.h"
#include "igl/sortrows.h"

namespace {
    
    double cot(const DGP::Vec& a, const DGP::Vec& b)
    {
        double cos = a.dot(b);
		DGP::Vec cross(3);
		cross << a(1)*b(2) - a(2)*b(1), a(2)*b(0) - a(0)*b(2), a(0)*b(1) - a(1)*b(0);
        double sin = cross.norm();
        return cos / sin;
    }
}

namespace DGP
{
    IGL_INLINE FMat dualHalfEdge(const FMat& F)
    {
        /* each row is (e_start, e_targ, tri_id, v_rank, v_id) */
        FMat halfedge_pairs(6 * F.rows(), 2 + 1 + 1 + 1);

        int count = 0;
        for(int i = 0; i < F.rows(); i++)
            for(int j = 0; j < 3; j++)
            {
                int e_start = i * 3 + j;
                int e_targ = i * 3 + (j + 1) % 3;
                int v_rank = (j + 2) % 3;
                int v_id = F(i, v_rank);

                halfedge_pairs.row(count++) << e_start, e_targ, i, v_rank, v_id;
                halfedge_pairs.row(count++) << e_targ, e_start, i, v_rank, v_id;
            }

        FMat sorted_hp, I;
        igl::sortrows(halfedge_pairs, true, sorted_hp, I);

        FMat dualHE = - FMat::Ones(3 * F.rows(), 3);
        int iter = 1;
        while (iter < sorted_hp.rows())
        {
            if (sorted_hp(iter, 0) == sorted_hp(iter - 1, 0) && 
                sorted_hp(iter, 1) == sorted_hp(iter - 1, 1))
            {
                int item_id0 = 3 * sorted_hp(iter, 2) + sorted_hp(iter, 3);
                dualHE.col(item_id0) << sorted_hp.col(iter-1).tail(3);

                int item_id1 = 3 * sorted_hp(iter, 2) + sorted_hp(iter, 3);
                dualHE.col(item_id1) << sorted_hp.col(iter).tail(3);

                iter += 2;
            }
            else { iter ++; }
        }
        
        return dualHE;
    }

    IGL_INLINE SpMat HodgeStar0Form(const VMat& V, const FMat& F)
    {
        Vec area;
        igl::doublearea(V, F, area);
        area /= 3.0;
    
        Vec v_area = Vec::Zero(V.rows());
        for (int i = 0; i < F.rows(); i++)
            for(int j = 0; j < 3; j++)
                v_area(F(i,j)) += area(i);

        SpMat hs0(V.rows(), V.rows());
        hs0 = v_area.asDiagonal();

        return hs0;
    }

    IGL_INLINE SpMat HodgeStar1Form(const VMat& V, const FMat& F)
    {
        Vec he_cot(3 * F.rows());

        int count = 0;
        for(int i = 0; i < F.rows(); i++)
        {
            Vec v0 = V.row(F(i, 0));
            Vec v1 = V.row(F(i, 1));
            Vec v2 = V.row(F(i, 2));

            he_cot(count++) = cot(v1-v0, v2-v0);
            he_cot(count++) = cot(v2-v1, v0-v1);
            he_cot(count++) = cot(v0-v2, v1-v2);
        }

        SpMat hs1(3 * F.rows(), 3 * F.rows());
        hs1 = he_cot.asDiagonal();

        return hs1;
    }

    IGL_INLINE SpMat HodgeStar2Form(const VMat& V, const FMat& F)
    {
        Vec area;
        igl::doublearea(V, F, area);
        area = 2.0 / area.array();
        
        SpMat hs2(V.rows(), V.rows());
        hs2 = area.asDiagonal();

        return hs2;
    }

    IGL_INLINE SpMat ReverseHodgeStar0Form(const VMat& V, const FMat& F)   /* *0^{-1} */
    {
        Vec area;
        igl::doublearea(V, F, area);
        area /= 3.0;
    
        Vec v_area = Vec::Zero(V.rows());
        for (int i = 0; i < F.rows(); i++)
            for(int j = 0; j < 3; j++)
                v_area(F(i,j)) += area(i);

        v_area = 1.0 / v_area.array();
        SpMat rhs0(V.rows(), V.rows());
        rhs0 = v_area.asDiagonal();

        return rhs0;
    }

    IGL_INLINE SpMat ReverseHodgeStar1Form(const VMat& V, const FMat& F)   /* *1^{-1} */
    {
        FMat dualHE = dualHalfEdge(F);
        
        std::vector<T> rhs1_coeff;
        rhs1_coeff.clear();
        rhs1_coeff.reserve(3 * F.rows());

        for (int i = 0; i < F.rows(); i++)
            for(int j = 0; j < 3; j++)
            {
                int he_id = 3 * i + j;
                if (dualHE(he_id, 0) != -1)
                {
                    int targ_he_id = dualHE(he_id, 1) * 3 + dualHE(he_id, 2);

                    iVec v_id = F.row(i);
                    Vec a = V.row(v_id((j+1)%3)) - V.row(v_id(j));
                    Vec b = V.row(v_id((j+2)%3)) - V.row(v_id(j));
                    double ratio = 1.0 / cot(a, b);

                    rhs1_coeff.push_back(T(targ_he_id, he_id, ratio));
                }
            }

        SpMat rhs1(3 * F.rows(), 3 * F.rows());
        rhs1.setFromTriplets(rhs1_coeff.begin(), rhs1_coeff.end());

        return rhs1;
    }

    IGL_INLINE SpMat ReverseHodgeStar2Form(const VMat& V, const FMat& F)
    {
        Vec area;
        igl::doublearea(V, F, area);
        area /= 2.0;
        
        SpMat rhs2(V.innerSize(), V.innerSize());
        rhs2 = area.asDiagonal();

        return rhs2;
    }

    IGL_INLINE SpMat ExteriorDerivative0Form(const VMat& V, const FMat& F)
    {
        std::vector<T> d0_coeff;
        d0_coeff.clear();
        d0_coeff.reserve(F.rows() * 6);

        for (int i = 0; i < F.rows(); i++)
        {
            int v0, v1, v2;
            v0 = F(i, 0);
            v1 = F(i, 1);
            v2 = F(i, 2);
            
            int e0, e1, e2;
            e0 = 3 * i;
            e1 = 3 * i + 1;
            e2 = 3 * i + 2;

            /* e0 = 1/2 * (v2 - v1) */
            d0_coeff.push_back(T(e0, v1, -0.5));
            d0_coeff.push_back(T(e0, v2,  0.5));

            /* e1 = 1/2 * (v0 - v2) */
            d0_coeff.push_back(T(e1, v2, -0.5));
            d0_coeff.push_back(T(e1, v0,  0.5));

            /* e2 = 1/2 * (v1 - v0) */
            d0_coeff.push_back(T(e2, v0, -0.5));
            d0_coeff.push_back(T(e2, v1,  0.5));
        }

        SpMat D0(F.rows() * 3, V.rows());
        D0.setFromTriplets(d0_coeff.begin(), d0_coeff.end());

        return D0;
    }

    IGL_INLINE SpMat ExteriorDerivative1Form(const VMat& V, const FMat& F)
    {
        std::vector<T> d1_coeff;
        d1_coeff.clear();
        d1_coeff.reserve(6 * F.rows());

        FMat dualHE = dualHalfEdge(F);

        for(int i = 0; i < F.rows(); i++)
            for(int j = 0; j < 3; j++)
            {
                int he_id = 3 * i + j;
                d1_coeff.push_back(T(i, he_id, 1));
                
                /* if there is opposite triangle */
                if (dualHE(he_id, 0) != -1)
                {
                    d1_coeff.push_back(T(dualHE(he_id, 0), he_id, -1));
                }
                else
                {
                    /* Half edge on the boundary, add this to replace 
                     * missing dual halfedge due to boudnary.
                     * Thus making D1 * D0 == 0 every where 
                     */
                    d1_coeff.push_back(T(i, he_id, 1));
                }
            }

        SpMat D1(F.rows(), 3 * F.rows());
        D1.setFromTriplets(d1_coeff.begin(), d1_coeff.end());

        return D1;
    }
}
