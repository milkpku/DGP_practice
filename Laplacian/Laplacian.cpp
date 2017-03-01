#define M_PI 3.1415926
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;


int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        printf("usage: %s objfile\n", argv[0]);
		return 0;
    }

    std::string filename(argv[1]);

    igl::readOBJ(filename, V, F);

    igl::viewer::Viewer viewer;

    viewer.data.set_mesh(V, F);
    viewer.launch();
}
