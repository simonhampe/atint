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

	This file provides functionality to compute certain special tropical varieties
	*/


#include "polymake/client.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {

	// PERL WRAPPER -------------------------------------------

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
									"# Creates the empty cycle in a given ambient dimension"
									"# (i.e. it will set the property [[PROJECTIVE_AMBIENT_DIM]]."
									"# @param Int ambient_dim The ambient dimension"
									"# @tparam Addition Max or Min"
									"# @return Cycle The empty cycle",
									"empty_cycle<Addition>($)");
	
	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
							"# Creates a cycle consisting of a collection of points"
							"# with given weights"
							"# @param Matrix<Rational> points The points, in tropical homogeneous coordinates"
							"# (though not with leading ones for vertices)."
							"# @param Vector<Integer> weights The list of weights for the points"
							"# @tparam Addition Max or Min"
							"# @return Cycle The point collection.",
							"point_collection<Addition>(Matrix<Rational>, Vector<Integer>)");

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
			"# Creates the linear space of the uniform matroid of rank k+1 on n+1 variables."
			"# @param Int n The ambient (projective) dimension."
			"# @param Int k The (projective dimension of the fan."
			"# @tparam Addition A The tropical addition (min or max)"
			"# @return Cycle A tropical linear space.",
			"uniform_linear_space<Addition>($,$)");       

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
									"# Creates a subdivision of the tropical projective torus"
									"# along an affine hyperplane into two halfspaces."
									"# This hyperplane is defined by an equation gx = a"
									"# @param Rational a The constant coefficient of the equation"
									"# @param Vector<Rational> g The linear coefficients of the equation"
									"# Note that the equation must be homogeneous in the sense that (1,..1)"
									"# is in its kernel, i.e. all entries of g add up to 0."
									"# @param Integer The (constant) weight this cycle should have"
									"# @tparam Addition Max or Min"
									"# @return Cycle The halfspace subdivision",
									"halfspace_subdivision<Addition>($,Vector<Rational>,$)");

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
									"# Creates the tropical projective torus of a given dimension."
									"# In less fancy words, the cycle is the complete complex"
									"# of given (tropical projective) dimension n, i.e. R<sup>n</sup>"
									"# @param Int n The tropical projective dimension."
									"# @param Integer w The weight of the cycle. Optional and 1 by default."
									"# @tparam Addition Max or Min."
									"# @return Cycle The tropical projective torus.",
									"projective_torus<Addition>($;$=1)");
}}

