#include <stdlib.h>
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>

#include "../include/Types.h"
#include "../include/HodgeDecomposition.h"

Eigen::MatrixXd V_ori;
Eigen::MatrixXd V;
Eigen::MatrixXi F;
double STEP = 0.01;

namespace DGP
{

    VMat smoothMesh_DO(const VMat& V, const FMat& F, double step)
    {
        SpMat d0 = ExteriorDerivative0Form(V, F);
        SpMat hs1 = HodgeStar1Form(V, F);
        SpMat dual_d0 = DualExteriorDerivative0Form(V, F); 
        SpMat hs0 = HodgeStar0Form(V, F);

        
        SpMat Lap = dual_d0 * hs1 * d0;
        SpMat H = hs0 - step * Lap;
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.compute(H);

        VMat new_V = solver.solve(hs0 * V);
        return new_V;
    }
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	printf("press key: %c, id: %d\n", key, key);

	if (key == ' ')
	{
		V = DGP::smoothMesh_DO(V, F, STEP);
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		//viewer.core.align_camera_center(V, F);
	}

	if (key == 'R')
	{
		V = V_ori;
		viewer.data.clear();
		viewer.data.set_mesh(V, F);
		//viewer.core.align_camera_center(V, F);
	}

	return 0;
}

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        printf("usage: %s objfile\n", argv[0]);
		return 0;
    }

    std::string filename(argv[1]);

    igl::readOBJ(filename, V_ori, F);
	V = V_ori;

    igl::viewer::Viewer viewer;

    viewer.data.set_mesh(V, F);
	viewer.callback_key_down = &key_down;
    viewer.launch();
}
