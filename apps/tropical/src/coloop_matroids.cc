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

	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {


	// Creates a matroid of rank r on n elements such that
	// the given sets of elements is the set of coloops and the matroid
	// is a uniform matroid restricted to the complement of the coloops
	perl::Object coloop_uniform_matroid(int n, int r, Set<int> coloops) {
		Set<int> complement = sequence(0,n) - coloops;
		int missing_rank = r - coloops.size();
		Array<Set<int> > all_sets = all_subsets_of_k(complement, missing_rank);
		Vector<Set<int> > bases;
		for(int as = 0; as < all_sets.size(); as++) {
			bases |= (coloops + all_sets[as]);
		}

		perl::Object result("matroid::Matroid");
			result.take("N_ELEMENTS") << n;
			result.take("BASES") << bases;
		return result;
	}

	UserFunction4perl("",&coloop_uniform_matroid,"coloop_uniform_matroid($,$,$)");

}}
