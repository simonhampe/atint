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
#include "polymake/IncidenceMatrix.h"

#ifndef ATINT_CONVERTERS_H
#define ATINT_CONVERTERS_H

namespace polymake { namespace atint{

  /**
   * @brief Converts an IncidenceMatrix to a Vector<Set<int> >
   */
  inline Vector<Set<int> > incmatrixToVector(const IncidenceMatrix<> &m) {
    Vector<Set<int> > result;
    for(int r = 0; r < m.rows(); r++) {
	result |= m.row(r);
    }
    return result;
  }
  
}}

#endif