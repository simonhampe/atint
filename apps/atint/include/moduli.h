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
#include "polymake/Set.h"
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
 @brief Does exactly the same as count_mn_cones, but returns an int. This will only work for n <= 12, since larger values produce too large integers. In that case you have to use count_mn_cones
 */
int count_mn_cones_int(int n);

/**
 @brief Takes a Pruefer sequence encoding a combinatorial type of n-marked rational curve and decodes it into the edge partitions of the corresponding graph
 @param Vector<int> seq The Pruefer sequence. Should be of length n + (no of bounded edges -1) and should contain only entries in (n,..) and each entry should occur at least twice. 
 @param int  The number of leafs. If not given (or set to something negative), the function assumes that the Pruefer sequence is ordered and that the first entry is hence equal to the number of leaves.
 @return Vector<Set<int> > A list of the partitions each edge induces. The leaves are given with indices (0,..,n-1) and each set is given such that it doesn't contain (n-1).
 */
Vector<Set<int> > decodePrueferSequence(const Vector<int> &seq, int n = -1);

}}

#endif // ATINT_MODULI_H_