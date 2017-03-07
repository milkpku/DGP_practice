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
 *   V: Represents of 0-form on vertex. 
 *      R^|V|, Row i represents 0-form value on v_i.
 *
 *   E: Represents of 1-form on halfedge.
 *      R^|E| = R^|3F|, Row 3*i+j represents 1/2 of integration of 1-form along 
 *      the j th halfedge of orinted triangle i, in its interior.
 *
 *   F: Represents of 2-form on triangle face.
 *      R^|F|, row i represents integration of 2-form on triangle i.
 *   
 *   V*: Represents of 2-form on dual vertex.
 *      R^|V|, row i represents integration of 2-form on dual area of vertex i.
 *
 *   E*: Represents of 1-form on dual halfedge.
 *      R^|E| = R^|3F|, row 3*i+j represents 1/2 of integration of 1-form along
 *      the dual edge of j th halfedge of orinted triangle i, in the edge's source's
 *      dual area.
 *
 *   F*: Represents of 0-form on dual triangle face.
 *      R^|F|, row i represents 0-form value on dual vertex of triangle i.
 *
 *   d0:  d operator on 0-form represented by vertex, results in 1-form represented by halfedge. R^{|E|*|V|}.
 *   d1:  d operator on 1-form represented by halfedge, results in 2-form represented by triangle face. R^{|F|*|E|}.
 *   ~d0: d operator on 1-form represented by dual halfedge in dual graph, R^{|V|*|E|}.
 *   ~d1: d operator on 0-form represented by dual face in dual graph, R^{|E|*|F|}.
 *
 *   *0: hodge-star operator on 0-form represented by vertex, results in 2-form represented by
 *       dual vertex in dual graph, R^{|V|*|V|}.
 *   *1: hodge-star operator on 1-form represented by halfedge, results in 1-form represented by
 *       dual halfedge in dual graph, R^{|E|*|E|}.
 *   *2: hodge-star operator on 2-form represented by triangle face, results in 0-form represented by
 *       dual triangle face in dual graph, R^{|F|*|F|}.
 *   *0^{-1}: hodge-star operator on 2-form represented by dual vertex in dual graph, results in 0-form
 *       
 */

#include "Types.h"

#define IGL_INLINE inline

namespace DGP
{
    IGL_INLINE FMat dualHalfEdge(const FMat& F); 

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
        return -2 * ExteriorDerivative0Form(V, F).transpose();
    }

    IGL_INLINE SpMat DualExteriorDerivative1Form(const VMat& V, const FMat& F) /* ~d1 */
    {
        return 0.5 * ExteriorDerivative1Form(V, F).transpose();
    }

}

#ifndef IGL_STATIC_LIBRARY
#include "HodgeDecomposition.cpp"
#endif
