/*
 This program is free s*oftware; you can redistribute it and/or
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
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"

#ifndef ATINT_MODULI_H_
#define ATINT_MODULI_H_

namespace polymake { namespace atint {

/**
 @brief Takes an index n and computes a matrix that assigns to all pairs (i,j) for i != j in {0,...,n-1} their coordinate index in R^(n over 2), when ordered lexicographically
 @param int n
 @returns Matrix<int> The entry at position (i,j) (or (j,i)) contains the coordinate index of that pair
 */
Matrix<int> pair_index_map(int n);

//Documentation see perl wrapper
Integer count_mn_cones(int n);

/**
 @brief Does exactly the same as count_mn_cones, but returns an int. This will only work for n <= 12, since larger values produce too large integer. In that case you have to use count_mn_cones
 */
int count_mn_cones_int(int n);

}}

#endif // ATINT_MODULI_H_