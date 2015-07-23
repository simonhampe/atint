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

	Contains functions to compute the affine transform of a cycle 
	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/PowerSet.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {

	perl::ListReturn corank_one_matroids(int n_elements) {
		perl::ListReturn result;
		//Iterate the number of coloops
		for(int cl = 0; cl < n_elements -1; cl++) {
			Array<Set<int> > clsets = all_subsets_of_k(sequence(0,n_elements),cl);
			for(int s = 0; s < clsets.size(); s++) {
				Vector<Set<int> > bases;
				//Create missing uniform sets
				Set<int> remaining = sequence(0, n_elements) - clsets[s];
				Array<Set<int> > uniform_sets = all_subsets_of_k( remaining, n_elements - cl-1);
				for(int us = 0; us < uniform_sets.size(); us++) {
					bases |= (uniform_sets[us] + clsets[s]);
				}//END iterate uniform sets
				perl::Object matroid("matroid::Matroid");
				matroid.take("N_ELEMENTS") << n_elements;
				matroid.take("BASES") << bases;
				result << matroid;
			}//END iterate coloop sets
		}//END iterate number of coloops
		return result;
	}//END corank_one_matroids

	perl::ListReturn flat_as_coloop_matroids(perl::Object matroid) {
		perl::ListReturn result;
		//Extract values 
		int n_elements = matroid.give("N_ELEMENTS");
		perl::Object lattice = matroid.give("LATTICE_OF_FLATS");
		IncidenceMatrix<> flats = lattice.give("FACES");
		
		//Iterate flats
		for(int f = 0; f < flats.rows(); f++) {
			//Ignore the full set 
			if(flats.row(f).size() == n_elements) continue;

			Set<int> remaining = sequence(0,n_elements) - flats.row(f);
			Array<Set<int> > uniform_sets = all_subsets_of_k( remaining, n_elements - flats.row(f).size()-1);
			Vector<Set<int> > bases;
			for(int us = 0; us < uniform_sets.size(); us++) {
				bases |= (uniform_sets[us] + flats.row(f));
			}
			perl::Object matroid("matroid::Matroid");
				matroid.take("N_ELEMENTS") << n_elements;
				matroid.take("BASES") << bases;
			result << matroid;
		}//END iterate flats
		return result;
	}//END flat_as_coloop_matroids

	UserFunction4perl("",&corank_one_matroids, "corank_one_matroids($) : returns(@)");
	UserFunction4perl("",&flat_as_coloop_matroids, "flat_as_coloop_matroids(matroid::Matroid) : returns(@)");

}}
