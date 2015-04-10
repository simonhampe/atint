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

	Basic functionality concerning lattices.
	*/


#ifndef POLYMAKE_ATINT_LATTICE_H
#define POLYMAKE_ATINT_LATTICE_H

#include "polymake/Rational.h"
#include "polymake/Matrix.h"

namespace polymake { namespace tropical {

	/*
	 * @brief Takes a rational matrix and makes each row integer by 
	 * multiplying it with an appropriate integer
	 */
	Matrix<Integer> make_rowwise_integer(const Matrix<Rational> &m); 

	/*
	 * @brief Computes the lattice basis of a cone, given in not-tropically-homogeneous coordinates and whose 
	 * dimension is known.
	 * @return A lattice basis, given as row vectors of a matrix
	 */
	Matrix<Integer> lattice_basis_of_cone(const Matrix<Rational> &rays, const Matrix<Rational> &lineality, int dim);

}}

#endif
