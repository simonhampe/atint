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

Contains the definitions for basicoperations.cc
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"

#ifndef ATINT_DIVISOR_H
#define ATINT_DIVISOR_H

namespace polymake { namespace atint{
  
  /**
    @brief Takes a matrix and returns the row indices where the first coordinate is nonzero and where the first coordinate is zero in  two different sets
    @param Matrix<Rational> m The matrix whose rows we consider
    @param Set<int> affine A reference to a set where the indices of the affine rows are stored
    @param Set<int> directional A reference to a set where the indices of the directional rows are stored
    @param bool uses_homog When set to false, all rows are directional and affine only contains the index -1
    @return std::pair<Set<int>, Set<int> > The first set contains the row indices of rows that start with a nonzero entry, the second set is the complement
  */
  inline std::pair<Set<int>, Set<int> > separateRays(Matrix<Rational> m, Set<int> &affine, Set<int> &directional, bool uses_homog);
  
}}

#endif