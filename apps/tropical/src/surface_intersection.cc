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

	Contains functions to compute intersection products in tropical surfaces.
	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/linear_algebra_tools.h"
#include "polymake/tropical/lattice.h"
#include "polymake/tropical/specialcycles.h"
#include "polymake/tropical/convex_hull_tools.h"

namespace polymake { namespace tropical {

	/*
	 * @brief Takes a list of rays which add up to 0, whose span dimension is one less than their number and whose maximal minors are 1. Then it also
	 * takes a list of vectors lying in that span. It will then compute for each vector the unique representation as a positive linear combination of the rays.
	 * @param Matrix<Rational> rank_one_flats, as PRIMITIVE INTEGER
	 * @param Matrix<Rational> curve_rays, as PRIMITIVE INTEGER
	 * @return Matrix<Integer> Each row corresponds to the representation of the same row in curve_rays. Each column corresponds to a row of surface_rays and contains
	 * the corresponding positive coefficient.
	 * */
	Matrix<Integer> positive_decomposition(const Matrix<Rational> &rank_one_flats, const Matrix<Rational> &curve_rays) {
		Matrix<Integer> result(curve_rays.rows(), rank_one_flats.rows());
		//Iterate curve rays 
		for(int cr = 0; cr < curve_rays.rows(); cr++) {
			//Compute a linear representation of the vector in the rays
			Vector<Rational> linRep = linearRepresentation(curve_rays.row(cr), rank_one_flats);
			//Go through its entries. For each negative entry, add the absolute value to all other entries and set this one to 0.
			for(int entry = 0; entry < linRep.dim(); entry++) {
				if(linRep[entry] < 0) linRep = linRep - (linRep[entry] * ones_vector<Rational>(linRep.dim()));
			}
			result.row(cr) = Vector<Integer>(linRep);
		}
		return result;
	}//END positive_decomposition

	/*
	 * @brief Takes a list of rays of a curve given by their positive linear representations wrt to a rank-1-flat-vectors matrix and a list of weights of those rays. It then computes the degree of that curve.
	 * @param Matrix<Rational> curve_decompositions 
	 * @param Vector<Integer> curve_weights 
	 * @return Integer
	 */
	Integer degree_via_decomposition(const Matrix<Integer> &curve_decompositions, const Vector<Integer> &curve_weights) {
		Integer deg = 0;
		for(int i = 0; i < curve_decompositions.rows(); i++) {
			deg += curve_decompositions(i,0) * curve_weights[i];
		}
		return deg;
	}//END degree_in_uniform

	/*
	 * @brief Computes the intersection multiplicity of two fan curves in a surface that is GLnZ-isomorphic
	 * to a uniform surface
	 * @param Matrix<Integer> rank_one_flats An integer matrix, whose rows are rays in a surface corresponding to the rank one flats of a matroid realization of that surface.
	 * @param Matrix<Rational> curve_a_rays The rays of the first curve (not homog, no leading coord)
	 * @param Vector<Integer> curve_a_weights The weights of the rays of the first curve, in the same order as the rays
	 * @param Matrix<Rational> curve_b_rays The rays of the second curve (not homog, no leading coord)
	 * @param Vector<Integer> curve_b_weights The weights of the rays of the second curve, in the same order as the rays.
	 */
	Integer intersection_multiplicity_via_flats( Matrix<Rational> &rank_one_flats,
																 Matrix<Rational> &curve_a_rays,
																 const Vector<Integer> &curve_a_weights,
																 Matrix<Rational> &curve_b_rays,
																 const Vector<Integer> &curve_b_weights) {
		//Make everything integer
		rank_one_flats = Matrix<Rational>(makePrimitiveInteger(rank_one_flats));
		curve_a_rays = Matrix<Rational>(makePrimitiveInteger(curve_a_rays));
		curve_b_rays = Matrix<Rational>(makePrimitiveInteger(curve_b_rays));

		Matrix<Integer> curve_a_decompositions = positive_decomposition( rank_one_flats, curve_a_rays); 
		Matrix<Integer> curve_b_decompositions = positive_decomposition( rank_one_flats, curve_b_rays);

		Integer result = degree_via_decomposition( curve_a_decompositions, curve_a_weights) *
								degree_via_decomposition( curve_b_decompositions, curve_b_weights);


		//We iterate pairs of rays of the two curves 
		for(int aray = 0; aray < curve_a_rays.rows(); aray++) {
			Vector<Integer> amap = curve_a_decompositions.row(aray);
			for(int bray = 0; bray < curve_b_rays.rows(); bray++) {
				Vector<Integer> bmap = curve_b_decompositions.row(bray);
				Integer correction = 0;
				//Iterate pairs of rank one flats 
				for(pm::Subsets_of_k_iterator<const pm::Series<int,true> &> pair = entire(all_subsets_of_k( sequence(0,rank_one_flats.rows()),2)); !pair.at_end(); pair++) {
					Vector<int> pair_as_vector(*pair);
					correction = std::max( correction,
							std::min(
								amap[pair_as_vector[0]] * bmap[pair_as_vector[1]],
								amap[pair_as_vector[1]] * bmap[pair_as_vector[0]]
								)*curve_a_weights[aray]*curve_b_weights[bray]
							);
				}//END iterate pairs of rank one flats
				result -= correction;
			}//END iterate curve B rays
		}//END iterate curve A rays

		

		return result;
	}//END intersection_multiplicity_in_uniform


	template <typename Addition>
		perl::Object intersect_in_smooth_surface(perl::Object surface, perl::Object cycle_a, perl::Object cycle_b) {
			//Extract data 
			int dim_a = cycle_a.give("PROJECTIVE_DIM");
			int dim_b = cycle_b.give("PROJECTIVE_DIM");
			int ambient_dim = surface.give("PROJECTIVE_AMBIENT_DIM");

			//Basic sanity checks 
			if(dim_a + dim_b <= 1) return empty_cycle<Addition>(ambient_dim);
			if(dim_a > 2 || dim_b > 2) throw std::runtime_error("intersect_in_smooth_surface: Cycles dimension too large.");

			//If one is full-dimensional, return the other one with multiplied weights 
			Vector<Integer> weights_a = cycle_a.give("WEIGHTS");
			Vector<Integer> weights_b = cycle_b.give("WEIGHTS");
			if(dim_a == 2) {
				return cycle_b.CallPolymakeMethod("multiply_weights",weights_a[0]);
			}
			if(dim_b == 2) {
				return cycle_a.CallPolymakeMethod("multiply_weights", weights_b[1]);
			}
			
			//From here on we know we have two curves.
			
			//Refine both curves along the surface 
			perl::Object refined_cycle_a = CallPolymakeFunction("intersect_container",cycle_a, surface);
			perl::Object refined_cycle_b = CallPolymakeFunction("intersect_container",cycle_b, surface);

			Matrix<Rational> rays_a = refined_cycle_a.give("VERTICES");
				rays_a = tdehomog(rays_a);
			Matrix<Rational> rays_b = refined_cycle_b.give("VERTICES");
				rays_b = tdehomog(rays_b);
			Matrix<Rational> lin_a = refined_cycle_a.give("LINEALITY_SPACE");
				lin_a = tdehomog(lin_a);
			Matrix<Rational> lin_b = refined_cycle_b.give("LINEALITY_SPACE");
				lin_b = tdehomog(lin_b);
			IncidenceMatrix<> cones_a = refined_cycle_a.give("MAXIMAL_POLYTOPES");
			IncidenceMatrix<> cones_b = refined_cycle_b.give("MAXIMAL_POLYTOPES");

			//Now intersect the refined versions.
			fan_intersection_result ab_intersect = cdd_fan_intersection(
					rays_a, lin_a, cones_a, rays_b, lin_b, cones_b);








		}//END intersect_in_smooth_surface


	// --------------------- PERL WRAPPERS -----------------------------



}}
