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
	Copyright (C) 2011 - 2015, Simon Hampe <simon.hampe@googlemail.com>

	Functions to compute moduli spaces of rational curves.
	*/

#ifndef POLYMAKE_ATINT_MODULI_RATIONAL_H
#define POLYMAKE_ATINT_MODULI_RATIONAL_H

#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"

namespace polymake { namespace tropical {

	/**
	  @brief Takes a Pruefer sequence encoding a combinatorial type of n-marked rational curve and 
	  decodes it into the edge partitions of the corresponding graph
	  @param Vector<int> seq The Pruefer sequence. Should be of length 
	  n + (no of bounded edges -1) and should contain only entries in (n,..) and each entry 
	  should occur at least twice. 
	  @param int  The number of leafs. If not given (or set to something negative), the function assumes 
	  that the Pruefer sequence is ordered and that the first entry is hence equal to the number of leaves.
	  @return Vector<Set<int> > A list of the partitions each edge induces. 
	  The leaves are given with indices (0,..,n-1) and each set is given such that it doesn't contain (n-1).
	  */
	Vector<Set<int> > decodePrueferSequence(const Vector<int> &pseq, int n);

}}

#endif
