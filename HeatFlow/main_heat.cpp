#include <stdlib.h>
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/jet.h>

#include "../include/Types.h"
#include "../include/ExteriorCalculus.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;
double STEP = 0.01;

namespace DGP
{

	Vec heatTrans(const VMat& V, const FMat& F, const Vec& heat, double h)
	{
		SpMat d0 = ExteriorDerivative0Form(V, F);
		SpMat hs1 = HodgeStar1Form(V, F);
		SpMat dual_d0 = DualExteriorDerivative0Form(V, F);
		SpMat hs0 = HodgeStar0Form(V, F);


		SpMat Lap = dual_d0 * hs1 * d0;
		SpMat H = hs0 - h * Lap;
		Eigen::SimplicialLDLT<SpMat> solver;
		solver.compute(H);

		Vec new_heat = solver.solve(hs0 * heat);
		return new_heat;
	}

	/* return |\nabla U| on each face */
	Vec normFactor(const VMat& V, const FMat& F, const Vec& U)
	{
		Vec normfactor = Vec::Zero(F.rows());

		for (int i = 0; i < F.rows(); i++)
		{
			iVec v_id = F.row(i);
			double ea = U(v_id(1)) - U(v_id(0));
			double eb = U(v_id(2)) - U(v_id(0));

			Vec a_tmp = V.row(v_id(1)) - V.row(v_id(0));
            Eigen::Vector3d a;
			a << a_tmp(0), a_tmp(1), a_tmp(2);

			Vec b_tmp = V.row(v_id(2)) - V.row(v_id(0));
			Eigen::Vector3d b;
			b << b_tmp(0), b_tmp(1), b_tmp(2);

            double a_norm = a.norm();
            double b_norm = b.norm();

            double cos = a.dot(b) / (a_norm * b_norm);
            double sin = a.cross(b).norm() / (a_norm * b_norm);

            double ea_n = ea / a_norm;
            double eb_n = eb / b_norm;

            normfactor(i) = sqrt( ea_n * ea_n + eb_n * eb_n - 2 * ea_n * eb_n * cos) / sin;
		}

        return normfactor;
	}

    Vec solveGradientField(const VMat& V, const FMat& F, const Vec& U)
    {
        /* build 1-form gradient field */
        Vec he(3 * F.rows());
        Vec normfactor = normFactor(V, F, U);
        for(int i = 0; i < F.rows(); i++)
        {
            iVec v_id = F.row(i);
            he(i * 3    ) = (U(v_id(2)) - U(v_id(1))) / normfactor(i);
            he(i * 3 + 1) = (U(v_id(0)) - U(v_id(2))) / normfactor(i);
            he(i * 3 + 2) = (U(v_id(1)) - U(v_id(0))) / normfactor(i);
        }

	    SpMat d0 = ExteriorDerivative0Form(V, F);
		SpMat hs1 = HodgeStar1Form(V, F);
		SpMat dual_d0 = DualExteriorDerivative0Form(V, F);

        Eigen::SimplicialLDLT<SpMat> solver;
        solver.compute(dual_d0 * hs1 * d0);

        return solver.solve(dual_d0 * hs1 * he);
    }

    Vec solveGeodesicDistance(const VMat& V, const FMat& F, const int origin)
    {
        Vec heat = Vec::Zero(V.rows());
        heat(origin) = 1;

        heat = heatTrans(V, F, heat, 1e-2);

        return solveGradientField(V, F, heat);
    }
}

int main(int argc, char *argv[])
{

	if (argc < 2)
	{
		printf("usage: %s objfile\n", argv[0]);
		return 0;
	}

	std::string filename(argv[1]);

	igl::readOBJ(filename, V, F);

	/* compute heat gradient, then norm it */

    DGP::Vec geodec = DGP::solveGeodesicDistance(V, F, rand() % V.rows());
	double norm = (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
	geodec *= 8 / norm;

	DGP::Vec pattern = (geodec.array() - geodec.array().round()).abs();
	
    DGP::VMat C;
    igl::jet(pattern, true, C);

	igl::viewer::Viewer viewer;

	viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);
	viewer.launch();
}
