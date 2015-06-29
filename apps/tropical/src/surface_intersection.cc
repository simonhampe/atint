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

namespace polymake { namespace tropical {

	/*
	 * @brief Takes the coarse rays of a surface GLnZ-isomorphic to a uniform surface and a list of rays in that surface.
	 * It will then compute which rays span the minimal cone containing each vector and return the corresponding
	 * coefficients.
	 * @param Matrix<Rational> surface_rays, as PRIMITIVE INTEGER
	 * @param Matrix<Rational> curve_rays, as PRIMITIVE INTEGER
	 * @return Array<Map<int,Rational>> Maps the indices of the (at most two) rays spanning the minimal cone containing
	 * the vector to the respective coefficients.
	 */
	Array<Map<int, Integer> > positive_decomposition(const Matrix<Rational> &surface_rays, const Matrix<Rational> &curve_rays) {
		Array<Map<int, Integer> > result(curve_rays.rows());
		//Iterate curve rays 
		for(int cr = 0; cr < curve_rays.rows(); cr++) {
			bool found = false;
			//Go through all pairs of rays of the surface
			for(int i = 0; i < surface_rays.rows(); i++ && !found) {
				for(int j = i+1; j < surface_rays.rows(); j++ && !found) {
					//Compute a linear representation of the vector in these rays
					Vector<Rational> linRep = linearRepresentation(curve_rays.row(cr), surface_rays.minor( scalar2set(i) + scalar2set(j), All));
					//Check it it is in the span and all coefficients are positive.
					if(linRep.dim() != 0) {
						if(linRep[0] >= 0 && linRep[1] >= 0) {
							if(linRep[0] > 0) result[cr][i] = Integer(linRep[0]);
							if(linRep[1] > 0) result[cr][j] = Integer(linRep[1]);
							found = true;
						}
					}
				}//END iterate second ray
			}//END iterate first ray
		}
		return result;
	}//END positive_decomposition

	/*
	 * @brief Takes a list of rays of a fan
	 * curve in a uniform surface, given by their positive linear representations wrt to the rays and a list of weights of those rays. It then computes the degree of that curve.
	 * @param Matrix<Rational> surface_rays
	 * @param Array<Map<int,Rational> > curve_decompositions 
	 * @param Vector<Integer> curve_weights 
	 * @return Integer
	 */
	Integer degree_in_uniform(const Array<Map<int,Integer> > &curve_decompositions, const Vector<Integer> &curve_weights) {
		Integer deg = 0;
		for(int i = 0; i < curve_decompositions.size(); i++) {
			if(curve_decompositions[i].contains(0))
				deg += curve_decompositions[i][0] * curve_weights[i];
		}
		return deg;
	}//END degree_in_uniform

	/*
	 * @brief Computes the intersection multiplicity of two fan curves in a surface that is GLnZ-isomorphic
	 * to a uniform surface
	 * @param Matrix<Rational> The coarse rays of the surface, not homogeneous and without leading coordinate.
	 * @param Matrix<Rational> curve_a_rays The rays of the first curve (not homog, no leading coord)
	 * @param Vector<Integer> curve_a_weights The weights of the rays of the first curve, in the same order as the rays
	 * @param Matrix<Rational> curve_b_rays The rays of the second curve (not homog, no leading coord)
	 * @param Vector<Integer> curve_b_weights The weights of the rays of the second curve, in the same order as the rays.
	 */
	Integer intersection_multiplicity_in_uniform( Matrix<Rational> &surface_rays,
																 Matrix<Rational> &curve_a_rays,
																 const Vector<Integer> &curve_a_weights,
																 Matrix<Rational> &curve_b_rays,
																 const Vector<Integer> &curve_b_weights) {
		//Make everything integer 
		surface_rays = Matrix<Rational>(makePrimitiveInteger(surface_rays));
		curve_a_rays = Matrix<Rational>(makePrimitiveInteger(curve_a_rays));
		curve_b_rays = Matrix<Rational>(makePrimitiveInteger(curve_b_rays));

		Array<Map<int,Integer> > curve_a_decompositions = positive_decomposition( surface_rays, curve_a_rays); 
		Array<Map<int,Integer> > curve_b_decompositions = positive_decomposition( surface_rays, curve_b_rays);

		Integer result = degree_in_uniform( curve_a_decompositions, curve_a_weights) *
								degree_in_uniform( curve_b_decompositions, curve_b_weights);

		//Iterate maximal cones of the surface:
		for(pm::Subsets_of_k_iterator<const pm::Series<int,true> &> face = entire(all_subsets_of_k( sequence(0, surface_rays.rows()),2)); !face.at_end(); face++) {
			Vector<int> face_as_vector( *face);
			//Iterate rays of curve A
			for(int ray_a = 0; ray_a < curve_a_rays.rows(); ray_a++) {
				Map<int,Integer> amap = curve_a_decompositions[ray_a];
				if(amap.contains(face_as_vector[0]) && amap.contains(face_as_vector[1])) {
					//Iterate rays of curve B
					for(int ray_b = 0; ray_b < curve_b_rays.rows(); ray_b++) {
						Map<int,Integer> bmap = curve_b_decompositions[ray_b];
						if(bmap.contains(face_as_vector[0]) && bmap.contains(face_as_vector[1])) {
							result -= curve_a_weights[ray_a] * curve_b_weights[ray_b] * 
								std::min( 
										amap[face_as_vector[0]] * bmap[face_as_vector[1]],
										amap[face_as_vector[1]] * bmap[face_as_vector[0]]
										);
						}
					}//END iterate B-rays
				}
			}//END iterate A-rays
		}//END iterate faces
		return result;
	}//END intersection_multiplicity_in_uniform


	// --------------------- PERL WRAPPERS -----------------------------
	
	Function4perl( &positive_decomposition, "pd(Matrix, Matrix)");
	Function4perl( &degree_in_uniform, "diu(Array<Map<Int,Integer> >, Vector<Integer>)");
	
	Function4perl( &intersection_multiplicity_in_uniform, "imu( Matrix,Matrix, Vector<Integer>, Matrix, Vector<Integer>)");

}}
