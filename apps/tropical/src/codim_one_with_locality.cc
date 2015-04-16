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

	Implements cartesian_product.h
	*/

#include "polymake/tropical/codim_one_with_locality.h"
#include "polymake/tropical/solver_def.h"
#include "polymake/tropical/separated_data.h"

namespace polymake { namespace tropical {

	CodimensionOneResult calculateCodimOneData(const Matrix<Rational> &rays, const IncidenceMatrix<> &maximalCones, const Matrix<Rational> &linspace, const IncidenceMatrix<>  &local_restriction) {

		//dbgtrace << "Computing all facets..." << endl;

		//First we construct the set of all facets 
		//Array<IncidenceMatrix<> > maximal_cone_incidence = fan.give("MAXIMAL_CONES_INCIDENCES");
		//Compute the rays-in-facets for each cone directly
		Vector<IncidenceMatrix<> > maximal_cone_incidence;
		for(int mc = 0; mc < maximalCones.rows(); mc++) {
			Set<int> mset = maximalCones.row(mc);
			//Extract inequalities
			//dbgtrace << "Computing facets for cone set " << mset << endl;
			Matrix<Rational> facets = solver<Rational>().enumerate_facets(
					rays.minor(mset,All),
					linspace,false,false).first;
			//dbgtrace << "Done. Checking rays..." << endl;
			//For each inequality, check which rays lie in it
			Vector<Set<int> > facetIncidences;
			for(int row = 0; row < facets.rows(); row++) {
				Set<int> facetRays;
				for(Entire<Set<int> >::iterator m = entire(mset); !m.at_end(); m++) {
					if(facets.row(row) * rays.row(*m) == 0) {
						facetRays += *m;
					}
				}
				facetIncidences |= facetRays;
			}
			//dbgtrace << "Done." << endl;
			maximal_cone_incidence |= IncidenceMatrix<>(facetIncidences);
		}

		//dbgtrace << "Check for doubles and useless facets..." << endl;

		//This will contain the set of indices defining the codim one faces
		Vector<Set<int> > facetArray;

		//This will define the codim-1-maximal-cone incidence matrix
		Vector<Set<int> > fIncones;

		for(int maxcone = 0; maxcone < maximal_cone_incidence.size(); maxcone++) {
			//This is the incidence matrix for the maximal cone indexed by maxcone
			IncidenceMatrix<> fcts = maximal_cone_incidence[maxcone];
			for(int facet = 0; facet < fcts.rows(); facet++) {
				Set<int> facetToCheck = fcts.row(facet);
				//If there is a local restriction, check if the facet is compatible
				if(local_restriction.rows() > 0) {
					if(!is_coneset_compatible(facetToCheck, local_restriction)) continue;
				}
				//Check if this facet intersects x0 = 1, otherwise go to the next one 
				//More precisely: Check if at least one of its rays has x0-coord != 0
				Vector<Rational> firstColumn = rays.minor(facetToCheck,All).col(0);
				if(firstColumn == zero_vector<Rational>(firstColumn.dim())) {
					continue;
				}
				//Otherwise check if we already have that facet and remember its index
				int fcIndex = -1;
				for(int existing = 0; existing < facetArray.dim(); existing++) {
					if(facetArray[existing] == facetToCheck) {
						fcIndex = existing;
						break;
					}
				}
				//Add the facet if necessary and add its maximal-cone indices
				if(fcIndex == -1) {
					facetArray |= facetToCheck;
					Set<int> singlecone;
					singlecone = singlecone + maxcone;
					fIncones |= singlecone;
				}
				else {
					fIncones[fcIndex] = fIncones[fcIndex] + maxcone;
				}
			}
		}

		CodimensionOneResult r;
		r.codimOneCones = IncidenceMatrix<>(facetArray);
		r.codimOneInMaximal = IncidenceMatrix<>(fIncones);
		return r;    
	}


}}

