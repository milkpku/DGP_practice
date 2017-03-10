#pragma once

#include "Types.h"
#include "ExteriorCalculus.h"

namespace DGP
{

    static SpMat d0;
    static SpMat d1;
    static SpMat dual_d0;
    static SpMat dual_d1;

    static SpMat hs0;
    static SpMat hs1;
    static SpMat hs2;
    static SpMat rhs0;
    static SpMat rhs1;
    static SpMat rhs2;

    void initHodgeDecomposition(const VMat& V, const FMat& F);

    Vec computeZeroFormPotential(const Vec& H);

    Vec computeTwoFormPotensial(const Vec& H);

    Vec extractHarmonicPart(const Vec& H);
}
