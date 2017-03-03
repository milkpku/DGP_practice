#include <stdlib.h>
#include <stdio.h>

#define _USE_MATH_DEFINES
#include <cmath>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>

#include "../include/Types.h"
#include "../include/Parameterization.h"

DGP::VMat V;
DGP::FMat F;
DGP::TMat T;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("usage: %s objfile\n", argv[0]);
        return -1;
    }

    /* read in objfile */
    igl::readOBJ(argv[1], V, F);
    
    T = DGP::conformalParameterization(V, F, 0.01);
	T *= 50;

    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V, F);
    viewer.data.set_uv(T);
	viewer.core.show_texture = true;
    viewer.launch(true);
}
