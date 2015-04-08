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

	These functions create certain special cycles. For documentation see
	the corresponding perl wrappers.
	*/


#ifndef POLYMAKE_ATINT_SPECIAL_CYCLES_H
#define POLYMAKE_ATINT_SPECIAL_CYCLES_H

#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Integer.h"

namespace polymake { namespace tropical {

	template <typename Addition>
		perl::Object empty_cycle(int ambient_dim); 


	template <typename Addition>
		perl::Object point_collection(Matrix<Rational> m, Vector<Integer> weights);


	template <typename Addition> 
		perl::Object uniform_linear_space(const int n, const int k);


	template <typename Addition>
		perl::Object halfspace_subdivision(Rational a, Vector<Rational> g, Integer weight);


	template <typename Addition>
		perl::Object projective_torus(int n, Integer weight=1); 
}}

#endif
