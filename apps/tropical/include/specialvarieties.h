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

#ifndef ATINT_SPECIAL_VARIETIES_H
#define ATINT_SPECIAL_VARIETIES_H

/**
 * @brief Creates the linear tropical space L^n_k. This tropical fan is defined in the following way: 
 * As rays we take -e_i,i=1,...,n, where e_i is the i-th standard basis vector of R^n and 
 * e_0 = e_1 + ... + e_n. As maximal cones we take the cones generated by rays {e_i, i in S}, where
 * S runs over all k-element subsets of {0,..,n}.
 * @param n The ambient dimension of the fan.
 * @param k The dimension of the fan (should be smaller equal n, otherwise an error is thrown).
 * If k is zero, this creates the fan consisiting of the origin in R^n
 * @return A PolyhedralFan object representing L^n_k
 */
perl::Object tropical_lnk(const int &n, const int &k);

/**
@brief Creates the bergman fan of a given matroid fan.
@param matroid::Matroid m A matroid
@param bool modOutLineality Optional argument. If set to TRUE, the lineality space is divided out before returning the 
fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE
@param int projectionCoordinate Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is 0.
@return fan::PolyhedralFan The bergman fan of the matroid, possibly with the lineality space divided out
*/
// perl::Object bergman_fan_via_polytope(perl::Object polytope, int rank, bool modOutLineality = false, int projectionCoordinate = 0);

#endif