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

	Defines a function to compute the minimal interior cells of a polyhedral subdivision.
	*/

#ifndef POLYMAKE_ATINT_MINIMAL_INTERIOR_H
#define POLYMAKE_ATINT_MINIMAL_INTERIOR_H

#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/tropical/solver_def.h"

namespace polymake { namespace tropical{

	/*
	 * @brief Computes the minimal interior cells of a polyhedral subdivision of a polyhedron.
	 * @param Matrix<Rational> vertices The vertices of the subdivision.
	 * @param IncidenceMatrix<> polytopes The subdivision cells 
	 * @param IncidenceMatrix<> The minimal interior cells of thei subdivision.
	 * @param solver<Rational> a convex hull solver to compute the facets of the cell
	 */
	IncidenceMatrix<> minimal_interior(const Matrix<Rational> &vertices, 
			const IncidenceMatrix<> &polytopes, solver<Rational> &sv);

}}

#endif
