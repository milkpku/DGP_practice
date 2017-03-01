#include <string>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>

#include "Laplacian.h"

Eigen::MatrixXd V_ori;
Eigen::MatrixXd V;
Eigen::MatrixXi F;


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	printf("press key: %c\n", key);

	if (key == ' ')
	{
		V = DGP::smoothMesh(V, F, 0.1);
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
