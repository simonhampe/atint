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

	These functions wrap the basic convex hull functionality of ppl, taking into account
	tropical homogeneous coordinates

*/


#ifndef POLYMAKE_ATINT_MISC_TOOLS_H
#define POLYMAKE_ATINT_MISC_TOOLS_H

#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"


namespace polymake { namespace tropical {

	/**
	  @brief Takes a matrix and returns the row indices where the first coordinate is nonzero and where the first coordinate is zero in  two different sets
	  @param Matrix<Rational> m The matrix whose rows we consider
	  @return std::pair<Set<int>, Set<int> > The first set contains the row indices of rows that start with a zero entry, the second set is the complement
	  */
	std::pair<Set<int>, Set<int> > far_and_nonfar_vertices(const Matrix<Rational> &m);

	/*
	 * @brief Takes a polyhedral complex and returns [[CONES]] summarized into one single incidence matrix.
	 * @param PolyhedralComplex
	 * @return IncidenceMatrix<>
	 */
	IncidenceMatrix<> all_cones_as_incidence(perl::Object complex);

	/*
	 * @brief Converts an incidence matrix to a Vector<Set<int> >
	 */
	Vector<Set<int> > incMatrixToVector(const IncidenceMatrix<> &i);

}}

#endif
