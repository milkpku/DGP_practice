#pragma once
/*
 * Discrete Differential Operators based on Halfedge concept
 *
 *         V  --- d0 -->  E  -- d1 -->  F
 *         |              |             | 
 *      *0 |           *1 |          *2 |
 *         v              v             v
 *         V* <-- ~d0 --- E* <-- ~d1 -- F*
 *
 *   V: R^|V|, row i represents 0-form value on v_i.
 *   E: R^|E| = R^|3F|, row 3*i+j represents 1/2 of integration of 1-form along 
 *      the j th halfedge of orinted triangle i, in its interior.
 *   F: R^|F|, row i represents integration of 2-form on triangle i.
 *   
 *   V*: R^|V|, row i represents integration of 2-form on dual area of vertex i.
 *   E*: R^|E| = R^|3F|, row 3*i+j represents 1/2 of integration of 1-form along
 *      the dual edge of j th halfedge of orinted triangle i, in the edge's source's
 *      dual area.
 *   F*: R^|F|, row i represents 0-form value on dual vertex of triangle i.
 *
 *
 */

#include "Types.h"

#define IGL_INLINE inline

namespace DGP
{

    IGL_INLINE SpMat HodgeStar0Form(const VMat& V, const FMat& F);          /* *0 */
    IGL_INLINE SpMat HodgeStar1Form(const VMat& V, const FMat& F);          /* *1 */
    IGL_INLINE SpMat HodgeStar2Form(const VMat& V, const FMat& F);          /* *2 */

    IGL_INLINE SpMat ReverseHodgeStar0Form(const VMat& V, const FMat& F);   /* *0^{-1} */
    IGL_INLINE SpMat ReverseHodgeStar1Form(const VMat& V, const FMat& F);   /* *1^{-1} */
    IGL_INLINE SpMat ReverseHodgeStar2Form(const VMat& V, const FMat& F);   /* *2^{-1} */

    IGL_INLINE SpMat ExteriorDerivative0Form(const VMat& V, const FMat& F); /* d0 */
    IGL_INLINE SpMat ExteriorDerivative1Form(const VMat& V, const FMat& F); /* d1 */

    IGL_INLINE SpMat DualExteriorDerivative0Form(const VMat& V, const FMat& F) /* ~d0 */
    {
        return 2 * ExteriorDerivative0Form(V, F).transpose();
    }

    IGL_INLINE SpMat DualExteriorDerivative1Form(const VMat& V, const FMat& F) /* ~d1 */
    {
        return 0.5 * ExteriorDerivative1Form(V, F).transpose();
    }

}
