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

	Reconstructs the matroid from a bergman fan.
	*/


#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {

	template <typename Addition>
	perl::Object matroid_from_fan(perl::Object cycle) {
		//Find rank and ground set 
		int ambient_dim = cycle.give("PROJECTIVE_AMBIENT_DIM");
		int n = ambient_dim+1;
		int dim = cycle.give("PROJECTIVE_DIM");
		int r = dim+1;

		if(dim == ambient_dim) {
			return CallPolymakeFunction("matroid::uniform_matroid",n,n);
		}

		//Take all r-sets and check if they are a basis
		Array<Set<int> > rset = all_subsets_of_k( sequence(0,n),r);
		Vector<Set<int> > bases;
		Matrix<Rational> unitm = unit_matrix<Rational>(n);
		for(int b = 0; b < rset.size(); b++) {
			perl::Object hp = affine_linear_space<Addition>(unitm.minor(~rset[b],All));
			perl::Object inter = CallPolymakeFunction("intersect", cycle, hp);
			if(CallPolymakeFunction("is_empty",inter)) bases |= rset[b];
		}
		perl::Object result("matroid::Matroid");
			result.take("N_ELEMENTS") << n;
			result.take("BASES") << bases;
		return result;
	}

	UserFunctionTemplate4perl("# @category Matroids"
			"# Takes the bergman fan of a matroid and reconstructs the corresponding matroid"
			"# The fan has to be given in its actual matroid coordinates, not as an isomorphic"
			"# transform. The actual subdivision is not relevant."
			"# @param Cycle<Addition> A tropical cycle, the Bergman fan of a matroid"
			"# @tparam Addition Min or Max"
			"# @return matroid::Matroid",
			"matroid_from_fan<Addition>(Cycle<Addition>)");


}}
