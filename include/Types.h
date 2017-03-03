#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace DGP
{
    /* geometry property */
    typedef Eigen::MatrixXd VMat;       /* vertex point position */
    typedef Eigen::MatrixXi FMat;       /* triangle vertex id triplet */
    typedef Eigen::MatrixXi BMat;       /* boundary halfedge pairs */
    typedef Eigen::MatrixXd TMat;       /* texture coordinate */

    /* linear algebra types */
    typedef Eigen::VectorXd Vec;                /* double vector */
    typedef Eigen::VectorXi iVec;               /* int vector */
    typedef Eigen::VectorXcd cVec;              /* double complex vector */
    typedef Eigen::SparseMatrix<double> SpMat;  /* double sparse matrix*/
    typedef Eigen::Triplet<double> T;           /* double Spmat triplets*/

    typedef std::complex<double> complex;       /* double complex */
    typedef Eigen::SparseMatrix<complex> SpMatC;/* double complex SpMat */
}
