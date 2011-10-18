/*
 This program is free software; you can redistribute it and/or
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
 */

#include "polymake/client.h"
#include "polymake/Array.h"

#ifndef ATINT_DIVISOR_H
#define ATINT_DIVISOR_H

namespace polymake { namespace atint {

/**
 @brief This function can be used in two ways (which are chosen via the boolean value isFunction): You can either give it a weighted complex and a rational function or a weighted complex and a polyhedral complex. In the first case it computes a refinement of the complex such that the function is cellwise linear and a new function that is just the restriction of the old one onto the new refinement. In the second case it assumes the variety is contained in the complex and computes a refinement. The function assumes that both variety and container have the same non-homogeneous dimension. I.e. a complex in homog. coordinates in R^3 (i.e. with non-homog. dimension 4) and a fan in non-homog. coordinates in R^3 would not be compatible, the fan would need to be homogenized first.
 @param perl::Object variety A WeightedComplex object. It needn't actually have defined TROPICAL_WEIGHTS. In this case the resulting variety will not have any weights either.
 @param perl::Object container Either a RationalFunction object or a WeightedComplex object. In either case, the function assumes that the variety is contained in the domain of the function or the given complex.
 @param bool isFunction This value indicates whether container is a RationalFunction (true) or a WeightedComplex (false)
 @param bool isMinMaxFunction This value is only relevant if isFunction is true. It indicates whether the rational function container is of type MinMaxFunction
 @return A perl::ListReturn object. In the first case (RationalFunction) it returns first the refined variety and second the new rational function. In the second case it only returns the refined variety.
 */
perl::ListReturn refine(perl::Object variety, perl::Object container, bool isFunction, bool isMinMaxFunction) ; 
  
/**
@brief Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete (i.e. its support is the whole ambient space) to compute the intersection correctly.
@param fan An arbitrary polyhedral fan
@param completeFan A complete polyhedral fan, in non-homog. coordinates
@return WeightedComplex The intersection of both fans (whose support is equal to the support of fan). The 
resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, 
the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned.
*/
perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan);

/**
  @brief Takes  as input a tropical fan / tropical variety and an array of rational values. The array length should coincide with the number of [[CMPLX_RAYS]] of the fan plus the dimension of the lineality space and will be interpreted as a rational function, where each value has been assigned to the rays given by $fan->CMPLX_RAYS and to the generators given by $fan->LINEALITY_SPACE (in that order). Missing values will be filled up by 0's, superfluous ones will be ignored. The function will then compute the corresponding Weil divisor and return it as a tropical variety given as a fan. The fan uses homogeneous coordinates, if and only the input fan does. 
  @param WeightedComplex A tropical variety on which the divisor is computed.
  @param Vector<Rational> values An array of rational values that define an integer affine map on the fan. 
  @return The divisor of the function defined by values on the given fan, as a tropical variety.
*/
perl::Object divisorByValueVector(perl::Object fan, Vector<Rational> values);

/**
  @brief Computes the divisor of a MinMaxFunction on a given tropical variety. The result will be in homogeneous coordinates, whether the tropical variety uses them or not. The function should be given on the affine coordinates of the variety, NOT the homogeneous ones.
  @param WeightedComplex fan A tropical variety, on which the divisor is computed
  @param MinMaxFunction A function whose DOMAIN should be equal to the affine coordinate space of the variety, i.e. AMBIENT_DIM-1, if the variety uses homogeneous coordinates, AMBIENT_DIM otherwise.
  @return The corresponding divisor as a tropical variety in homogeneous coordinates.
*/
perl::Object divisorByPLF(perl::Object fan, perl::Object);

/**
  @brief Computes the function value of a min-max function at a given point
  @param Matrix<Rational> functionMatrix The function matrix of the min-max-function. Each row corresponds to a function, the last column contains the constant coefficients
  @param Vector<Rational> point The point at which the function should be evaluated
  @param bool uses_min True, if we take the mininum over all function rows, otherwise we take the maximum
  @param bool uses_homog Whether point is given in homogeneous coordinates
*/
inline Rational functionValue(Matrix<Rational> functionMatrix, Vector<Rational> point, bool uses_min, bool uses_homog);

}}

#endif
