#include <stdlib.h>
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/lscm.h>
#include <igl/viewer/Viewer.h>

#include "../include/Types.h"
#include "../include/Parameterization.h"

DGP::VMat V;
DGP::FMat F;
DGP::TMat T;
DGP::TMat V_uv;
DGP::TMat T_uv;

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	printf("press key: %c, id: %d\n", key, key);

	if (key == '1')
	{
		viewer.data.set_uv(V_uv);
	}

	if (key == '2')
	{
		viewer.data.set_uv(T);
	}

	if (key == '3')
	{
		viewer.data.set_uv(T_uv);
	}


	return 0;
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("usage: %s objfile\n", argv[0]);
        return -1;
    }

    /* read in objfile */
    igl::readOBJ(argv[1], V, F);
    
	/* myself program */
    T = DGP::conformalParameterization(V, F, 1e-5);
	T *= 100;

	/* standard program */
	DGP::iVec bnd, b(2, 1);
	igl::boundary_loop(F, bnd);
	b(0) = bnd(0);
	b(1) = bnd(bnd.size() / 2);
	DGP::VMat bc(2, 2);
	bc << 0, 0, 1, 0;
	printf("fix ids: %d, %d\n", b(0), b(1));

	igl::lscm(V, F, b, bc, V_uv);
	V_uv *= 5;

	/* myself program 2 */
	T_uv = DGP::leastSquareConformalParameterization(V, F, b, bc);
	T_uv *= 5;


    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V, F);
    viewer.data.set_uv(T);
	viewer.core.show_texture = true;
	viewer.callback_key_down = &key_down;
    viewer.launch(true);
}
