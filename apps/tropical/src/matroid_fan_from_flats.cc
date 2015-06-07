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

	Computes a matroid fan in its chains-of-flats structure.
	*/

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {


	template <typename Addition>
		perl::Object matroid_fan_from_flats(perl::Object matroid) {
			//Extract properties 
			int n = matroid.give("N_ELEMENTS");
			int r = matroid.give("RANK");
			Set<int> loops = matroid.give("LOOPS");
			if(loops.size() != 0) {
				return empty_cycle<Addition>(n-1);
			}

			//Construct chains of flats 

			perl::Object flats = matroid.give("LATTICE_OF_FLATS");
			IncidenceMatrix<> faces = flats.give("FACES");
			IncidenceMatrix<> edges = flats.CallPolymakeMethod("EDGES");
			IncidenceMatrix<> nodes_in_edges = T(edges);

			//Prepare variables for chains, start with one chain with element 0
			//(that will correspond to the vertex later)
			//Each chain is an ordered list of indices, referring to FACES
			//Top element saves the index of the top rank flat of each chain. 
			Vector< Set<int> > chains_as_sets(1);
			Vector<int> top_element(1);
				chains_as_sets[0] = scalar2set(0); 
				top_element[0] = 0;

			for(int rk = 1; rk < r; rk++) {
				Set<int> flats_of_rank = flats.CallPolymakeMethod("nodes_of_dim",rk);
				Vector<Set<int> > new_chains;
				Vector<int> new_top_element;
				//For each chain, find each possible continuation.
				for(int chain = 0; chain < chains_as_sets.dim(); chain++) {
					//Find top neighbours of the top element 
					Set<int> connected_edges = nodes_in_edges.row(top_element[chain]);
					Set<int> neighbours = accumulate(rows( edges.minor(connected_edges,All)), operations::add());
						neighbours *= flats_of_rank;
					for(Entire<Set<int> >::iterator nb = entire(neighbours); !nb.at_end(); nb++) {
						new_chains |= ( chains_as_sets[chain] + *nb);
						new_top_element |= *nb;
					}
				}
				chains_as_sets = new_chains;
				top_element = new_top_element;
			}//END iterate ranks

			//Create rays 
			//The 0-th ray is the vertex
			Matrix<Rational> unitm = unit_matrix<Rational>(n);
				unitm = zero_vector<Rational>() | unitm;
			Matrix<Rational> rays(faces.rows()-1, n+1);
			rays(0,0) = 1;
			for(int f = 1; f < faces.rows()-1; f++) {
				rays.row(f) = Addition::orientation() * accumulate(rows(unitm.minor(faces.row(f),All)), operations::add());	
			}

			perl::Object result(perl::ObjectType::construct<Addition>("Cycle"));
				result.take("PROJECTIVE_VERTICES") << rays;
				result.take("MAXIMAL_POLYTOPES") << chains_as_sets;
				result.take("WEIGHTS") << ones_vector<Integer>(chains_as_sets.dim());
			return result;
		}


	UserFunctionTemplate4perl("# @category Matroids"
			"# Computes the fan of a matroid in its chains-of-flats subdivision."
			"# Note that this is potentially very slow for large matroids."
			"# @param matroid::Matroid A matroid. Should be loopfree."
			"# @tparam Addition Min or max, determines the matroid fan coordinates."
			"# @return Cycle<Addition>",
			"matroid_fan_from_flats<Addition>(matroid::Matroid)");

}}
