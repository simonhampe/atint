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

	Implementations of miscelleaneous tools
	*/

#include "polymake/client.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"


namespace polymake { namespace tropical { 

	std::pair<Set<int>, Set<int> > far_and_nonfar_vertices(const Matrix<Rational> &m) {
		Set<int> nonfar;
		Set<int> far;
		for(int r = 0; r < m.rows(); r++) {
			(m(r,0) == 0 ? far : nonfar) += r;
		}
		return std::pair<Set<int>,Set<int> >(far,nonfar);
	}


	IncidenceMatrix<> all_cones_as_incidence(perl::Object complex) {
		Array<IncidenceMatrix<> > all_cones = complex.give("CONES");
		if(all_cones.size() == 0) return IncidenceMatrix<>();
		IncidenceMatrix<> result(0,all_cones[0].cols());
		for(int i = 0; i < all_cones.size(); i++) {
			result /= all_cones[i];
		}
		return result;
	}


	Vector<Set<int> > incMatrixToVector(const IncidenceMatrix<> &i) {
		Vector<Set<int> > result;
		for(int r = 0; r < i.rows(); r++) {
			result |= i.row(r);
		}
		return result;
	}


	Array<Integer> randomInteger(const int& max_arg, const int &n) {
		static Integer upperBound = 0;
		static UniformlyRandomRanged<Integer> rg(max_arg);
		if(max_arg != upperBound)  {
			rg = UniformlyRandomRanged<Integer>(max_arg);
			upperBound = max_arg;
		}
		Array<Integer> result(n);
		for(int i = 0; i < n; i++) {
			result[i] = rg.get();
		}
		return result;
	}


  UserFunction4perl("# @category Lattices"
                  "# Returns n random integers in the range 0.. (max_arg-1),inclusive"
                  "# Note that this algorithm is not optimal for real randomness:"
                  "# If you change the range parameter and then change it back, you will"
                  "# usually get the exact same sequence as the first time"
                  "# @param int max_arg The upper bound for the random integers"
                  "# @param int n The number of integers to be created"
                  "# @return Array<Integer>",
                  &randomInteger,"randomInteger($, $)"); 
}}
