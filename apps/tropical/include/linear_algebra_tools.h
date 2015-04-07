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

	This file includes functions that perform basic linear algebra
	computations
	*/


#ifndef POLYMAKE_TROPICAL_LINEAR_ALGEBRA_TOOLS_H
#define POLYMAKE_TROPICAL_LINEAR_ALGEBRA_TOOLS_H

#include "polymake/Rational.h"
#include "polymake/Matrix.h"

namespace polymake { namespace tropical {

	/**
	  @brief  This method takes a set of row indices for [[SEPARATED_VERTICES]] and a vector that is supposed 
	  to be in the affine span of these row vectors and the lineality space (see description 
	  of [[LATTICE_NORMAL_FCT_VECTOR]]). It then computes the corresponding representation in these vectors
	  @param s a set of row indices of [[SEPARATED_VERTICES]]
	  @param v a vector supposed to lie in the affine span of [[SEPARATED_VERTICES]] + [[LINEALITY_SPACE]]
	  @param ambient_dim The ambient dimension of the cycle
	  @param rays The matrix of [[SEPARATED_VERTICES]]
	  @param linealitySpace A matrix of generators of the lineality space
	  @param lineality_dim The dimension of the lineality space
	  @return A vector of length [[SEPARATED_VERTICES]]->rows() + [[LINEALITY_DIM]] with linear coefficients 
	  of a representation in the generators chosen via s. The last elements always refer to the lineality space.
	  @throw std::runtime_error If the vector is not in the affine span of the given vectors
	  */
	Vector<Rational> functionRepresentationVector(const Set<int> &rayIndices, const Vector<Rational> &v,
			int ambient_dim, bool uses_homog, 
			const Matrix<Rational> &rays,
			const Matrix<Rational> &linealitySpace,
			int lineality_dim); 

}}

#endif
