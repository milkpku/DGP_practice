#pragma once

#include "../include/Parameterization.h"
#include <gtest/gtest.h>

namespace DGP {
	TEST(Parameterization, boundary) {

	}


	TEST(Parameterization, residual) {
		
        for(int i = 0; i < 100; i++)
        {
            int N = rand() % 30 + 10;
		
            SpMatC A = Eigen::MatrixXcd::Random(N, N).sparseView();
		
            Vec u = Vec::Random(N);
		    Vec v = Vec::Random(N);

            cVec y = u.cast<complex>() + complex(0, 1) * v.cast<complex>();
            cVec y_ = u.cast<complex>() - complex(0, 1) * v.cast<complex>();

            double e = residual(A, y);

            cVec res = A*y - (y_.transpose() * A * y) * y;
            double L2 = res.real().squaredNorm() + res.imag().squaredNorm();

            ASSERT_LT(abs(e*e - L2), 1e-6);
        }

	}

    TEST(Parameterization, area)
    {
        VMat V = VMat::Zero(4, 2);
        FMat F(2, 3);
        F << 0, 1, 2, 0, 2, 3;

        cVec z(4);
        z << complex(0, 0), complex(1, 0), complex(1, 1), complex(0, 1);

        SpMatC A = buildAreaMatrix(V, F).cast<complex>();

        ASSERT_EQ(complex(1, 0), complex(0, -0.25) * z.dot(A * z));
    }
}
