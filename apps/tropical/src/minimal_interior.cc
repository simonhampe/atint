
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

	Computes the minimal interior cells of a subdivision of a polyhedron
	*/

#include "polymake/client.h"
#include "polymake/Set.h"
#include "polymake/Array.h"
#include "polymake/Matrix.h"
#include "polymake/SparseMatrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Map.h"
#include "polymake/linalg.h"
#include "polymake/common/hermite_normal_form.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/lattice.h"
#include "polymake/tropical/solver_def.h"
#include "polymake/tropical/linear_algebra_tools.h"
#include "polymake/tropical/LoggingPrinter.h"


namespace polymake { namespace tropical {

	using namespace atintlog::donotlog;
	// using namespace atintlog::dolog;
	// using namespace atintlog::dotrace;
	
	/*
	 * @brief Computes the set of indices such that the corresponding entries
	 * in a given vector are either all zero or nonzero.
	 * @param Vector v The vector
	 * @param bool search_zeros: If true, find all zero entries, if false all non-zero
	 * entries
	 * @return Set<int>
	 */
	template <typename T>
		Set<int> binaryFinder(const GenericVector<T> &v, bool search_zeros) {
			if(search_zeros) return indices( attach_selector(v.top(), operations::is_zero()));
			else return indices( attach_selector(v.top(), operations::non_zero())); 
		}

	/**
	 * @brief Finds all the maximal sets and returns the row indices
	 */
	Set<int> find_maximal_faces(const IncidenceMatrix<> mat) {
		Set<int> smaller_sets;
		for(int i = 0; i < mat.rows(); i++) {
			for(int j = 0; j < mat.rows(); j++) {
				if(i != j) {
					if( mat.row(i) * mat.row(j) == mat.row(i)) smaller_sets += i;
				}
			}
		}
		return sequence(0,mat.rows()) - smaller_sets; 
	}//END find_maximal_faces

	inline Set<int> pairset(int i, int j) {
		return scalar2set(i) + scalar2set(j);
	}

	IncidenceMatrix<> minimal_interior(const Matrix<Rational> &vertices,	
			const IncidenceMatrix<> &polytopes, solver<Rational> &sv) {
		IncidenceMatrix<> result(0,vertices.rows());

		//If it only has one cone, it is minimal
		if(polytopes.rows() == 1) {
			return IncidenceMatrix<>(polytopes);	
		}
	
		//Compute facet normals of the original polyhedron
		Matrix<Rational> facet_normals = sv.enumerate_facets(vertices, Matrix<Rational>(0,vertices.cols()),false,false).first;

		//The first step is computing interior and exterior codimension one cells

		//Exterior one are the maximal intersections of the facet normals with the polytopes
		IncidenceMatrix<> intersect_fn_with_polytopes(0,vertices.rows());
		for(int f = 0; f < facet_normals.rows(); f++) {
			//Find the vertices of the polytopes lying in the facet
			Set<int> facet_vertices = binaryFinder( vertices * facet_normals.row(f), true);
			for(int p = 0; p < polytopes.rows(); p++) {
				intersect_fn_with_polytopes /= facet_vertices * polytopes.row(p);
			}
		}
		IncidenceMatrix<> exterior_codim = intersect_fn_with_polytopes.minor( 
				find_maximal_faces(intersect_fn_with_polytopes),All);

		//Interior ones are the maximal intersections of the polytopes
		IncidenceMatrix<> pairwise_intersections(0,vertices.rows());
		IncidenceMatrix<> associated_pairs(0,polytopes.rows());
		for(int p1 = 0; p1 < polytopes.rows(); p1++) {
			for(int p2 = p1+1; p2 < polytopes.rows(); p2++) {
				pairwise_intersections /= (polytopes.row(p1)*polytopes.row(p2));
				associated_pairs /= pairset(p1,p2);
			}
		}
		Set<int> pairwise_maximal = find_maximal_faces(pairwise_intersections);
		IncidenceMatrix<> interior_codim = pairwise_intersections.minor(pairwise_maximal,All);
		//Tells for all polytopes, which interior codim 1 cells they contain.
		IncidenceMatrix<> interior_in_max = T(associated_pairs.minor(pairwise_maximal,All));

		//If there are no codimension one cones, return the maximal cells
		if(exterior_codim.rows() + interior_codim.rows() == 0) {
			return IncidenceMatrix<>(polytopes);
		}

		//For each maximal cone, we compute all its minimal interior faces as maximal intersections
		//of interior codimension one faces. However, we only use codim-1-faces, that we haven't used
		//in another maximal cone so far (Since any minimal face in such a face is hence also a minimal
		//face of this other maximal cone, so we already computed it).

		Set<int> markedFaces;

		for(int mc = 0; mc < polytopes.rows(); mc++) {
			//Compute all non-marked codim-1-cells of mc. If there are none left, go to the next cone
			Vector<int> nonmarked(interior_in_max.row(mc) - markedFaces);
			if(nonmarked.dim() == 0) continue;
			//dbgtrace << "Remaining interior cells are: " << nonmarked << endl;
			int k = nonmarked.dim();
			//ordered list of indices of interior codim-1-cells (in nonmarked)
			//indices != -1 correspond to codim-1-cells that we intersect to obtain a minimal face
			Vector<int> currentSet(k,-1); 
			currentSet[0] = 0;
			//Indicates the index below the next cone index (in nonmarked) we should try to add. Is always
			//larger equal than the last element != -1 in currentSet
			int lowerBound = 0;
			//Indicates the current position in currentSet we're trying to fill
			int currentPosition = 1;
			//Indicates at position i < currentPosition the intersection of the cones specified by the
			// elements currentSet[0] .. currentSet[i]
			Vector<Set<int> > currentIntersections(k);
			currentIntersections[0] = interior_codim.row(nonmarked[0]);
			//Now iterate all intersections in a backtrack algorithm:
			//If an intersection is maximal, we don't need to go any further
			//We stop when we have tried all possibilities at the first position
			while(!(currentPosition == 0 && lowerBound == k-1)) {
				//Try the next posssible index
				int j = lowerBound+1;
				//If we're already beyond k-1, we have found a maximal intersection	
				//Check if it is a minimal face, then go back one step
				if(j == k) {
					//We test, if the set is not contained in any border face and does not contain
					//any existing minimal face
					Set<int> potentialMinimal(currentIntersections[currentPosition-1]);
					bool invalid = false;
					for(int ec = 0; ec < exterior_codim.rows(); ec++) {
						if(potentialMinimal.size() == (potentialMinimal * exterior_codim.row(ec)).size()) {
							invalid = true; break;
						}
					}
					for(int mf = 0; mf < result.rows() && !invalid; mf++) {
						if(result.row(mf).size() == (result.row(mf) * potentialMinimal).size()) {
							invalid = true; break;
						}
					}
					if(!invalid) result /= potentialMinimal;
					//Go back one step
					lowerBound = currentSet[currentPosition-1];
					currentPosition--;
				}
				else {
					//Compute the intersection with the next codim 1 cell
					//(If we're at the first position, we just insert a cell)
					Set<int> intersection;
					if(currentPosition == 0) intersection = interior_codim.row(nonmarked[j]);
					else intersection = 
						currentIntersections[currentPosition-1] * interior_codim.row(nonmarked[j]);
					//If its still not empty, go forward one step
					if(intersection.size() > 0) {
						currentSet[currentPosition] = j;
						currentIntersections[currentPosition] = intersection;
						currentPosition++;
						lowerBound = j;
					}
					else { //Otherwise go up one step
						lowerBound++;
					}
				}//END if j == k
			}//END while(...)

			//Now mark all codim -1-faces of mc
			markedFaces += interior_in_max.row(mc);

		}//END iterate maximal cones

		return result;
	}//END minimal_interior
/*
	IncidenceMatrix<> test(Matrix<Rational> m, IncidenceMatrix<> i) {
		solver<Rational> sv;
		return minimal_interior(m,i,sv);
	}

	Function4perl(&test, "minimal_interior(Matrix,IncidenceMatrix)");*/

}}
