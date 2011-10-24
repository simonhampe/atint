/*
 T his program is free s*oftware; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA  02110-1301, USA.
 
 ---
 Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 This file contains the definition of the generalized refinement function
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"

#ifndef ATINT_REFINE_h_
#define ATINT_REFINE_h_

namespace polymake { namespace atint { 

  struct RefinementResult {
	perl::Object complex;
	Matrix<Rational> rayRepFromX;
	Matrix<Rational> rayRepFromY;
	Matrix<Rational> linRepFromX;
	Matrix<Rational> linRepFromY;
	Vector<int> associatedRep;
  }; 
  
  /**
  @brief This is a multi-purpose function used to refine polyhedral complexes and rational functions along each other. It is not intended for direct use by the user but should be wrapped by other function according to each purpose. Its precise behavior is the following: It takes two WeightedComplex objects, X and Y, of which we assume that |X| is contained in |Y|. Furthermore, if Y has homog. coordinates, then so does X (the converse need not be true). It then computes a weighted complex X' that has the same support as X, but each cone of which is contained in a cone of Y. It can also compute the following data: For each ray (lineality space generator) of X' it can compute an appropriate linear representation in rays (and linspace generators) of X and/or Y, which are in a cone containing that ray. This can then be uses to compute function values of a function defined on X or Y on the new complex X'. Furthermore, if X uses homogeneous coordinates, it can compute for each directional ray the index of a vertex sharing a cone with that ray.
  @param perl::Object X A WeightedComplex (not necessarily with defined TROPICAL_WEIGHTS)
  @param perl::Object Y A WeightedComplex (not necessarily with defined TROPICAL_WEIGHTS), such that |X| is contained in |Y|
  @param bool repFromX Whether a representation of the new rays in the generators of X should be computed
  @param bool repFromY Whether a representation of the new rays in the generators of Y should be computed
  @param bool computeAssoc Whether an affine ray sharing a cone should be found for each directional ray
  @param bool refine Whether we actually need to refine X (true) or whether X is already fine in Y (false)
  @return RefinementResult A struct containing the following data (the rep values are only defined, if the corresponding boolean flag was set):
  1) perl::Object complex: The new refined complex X'. It has homog. coordinates, iff X has and it has weights, iff X does
  2) Matrix<Rational> rayRepFromX, linRepFromX: A row in these matrices corresponds to a ray/lin.space generator of X' and gives the coefficients in a linear representation in terms of the rays/lin.space generators of X
  3) Matrix<Rational> rayRepFromY, linRepFromY: Same as 2),only for Y
  4) Vector<int> associatedRep Gives for each ray the index of a vertex sharing a cone. For each vertex i, associatedRep[i] = i
  */
  RefinementResult refinement(perl::Object X, perl::Object Y, bool repFromX, bool repFromY,bool computeAssoc,bool refine);

}}

#endif // ATINT_REFINE_h_