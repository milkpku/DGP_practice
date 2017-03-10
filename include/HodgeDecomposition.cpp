#include <stdlib.h>
#include <stdio.h>

#include "Types.h"
#include "HodgeDecomposition.h"

namespace DGP
{
    void initHodgeDecomposition(const VMat& V, const FMat& F)
    {
        d0 = ExteriorDerivative0Form(V, F);
        d1 = ExteriorDerivative1Form(V, F);
        dual_d0 = DualExteriorDerivative0Form(V, F);
        dual_d1 = DualExteriorDerivative1Form(V, F);

        hs0 = HodgeStar0Form(V, F);
        hs1 = HodgeStar1Form(V, F);
        hs2 = HodgeStar2Form(V, F);
        rhs0 = ReverseHodgeStar0Form(V, F);
        rhs1 = ReverseHodgeStar1Form(V, F);
        rhs2 = ReverseHodgeStar2Form(V, F);
    }

    /* d*d a = d*w */
    Vec computeZeroFormPotential(const Vec& W)
    {
        SpMat L = dual_d0 * hs1 * d0;
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.compute(L);

        return solver.solve( dual_d0 * hs1 * W)
    }

    Vec computeTwoFormPotensial(const Vec& W)
    {
        SpMat L = d1 * rhs1 * dual_d1;
        Eigen::SimplicialLDLT<SpMat> solver;
        solver.compute(L);

        return solver.solve(d1 * W);
    }

    Vec extractHarmonicPart(const Vec& W)
    {
        Vec alpha = computeZeroFormPotential(W);
        Vec dual_beta = computeTwoFormPotensial(W);

        Vec gamma = W - d0 * alpha - dual_d1 * rhs1 * dual_beta;

        return gamma;
    }
}
