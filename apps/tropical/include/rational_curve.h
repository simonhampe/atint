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

	Contains functions to convert rational curves in their various
	representations.
	*/

#ifndef POLYMAKE_ATINT_RATIONAL_CURVE_H
#define POLYMAKE_ATINT_RATIONAL_CURVE_H

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"

namespace polymake { namespace tropical {

	///////////////////////////////////////////////////////////////////////////////////////

	/*
	 * @brief Computes n from (n choose 2).
	 */
	int moduliDimensionFromLength(int length);



	//Documentation see perl wrapper of wrapTestFourPointCondition (does the same except that it returns
	// a vector of int)
	Vector<int> testFourPointCondition(Vector<Rational> v) ;



	/**
	  @brief Computes a rational curve (in the v_I-representation) and its graph from a given metric (or more precisely a vector equivalent to a metric). It is wrapped by curveFromMetric and graphFromMetric, which should be called instead and whose documentation can be found in the corr. perl wrappers rational_curve_from_metric and curve_graph_from_metric. Note that the order of the [[EDGES]] of the graph object that do not contain leaf vertices is the same as the order of the corresponding [[SETS]] of the curve. Furthermore, the first n vertices of the graph are the end vertices of the leave edges (in order 1 .. n)
	  */
	perl::Object curveAndGraphFromMetric(Vector<Rational> metric);


	
	//Documentation see perl wrapper
	perl::Object curveFromMetric(Vector<Rational> metric); 



	//Documentation see perl wrapper
	template <typename Addition>
		perl::Object rational_curve_from_matroid_coordinates(Vector<Rational> matroidVector) {

			matroidVector = matroidVector.slice(~scalar2set(0));

			//Convert vector to a map
			int n = moduliDimensionFromLength(matroidVector.dim())+1;
			Matrix<Rational> d(n,n);
			int index = 0;
			for(int i = 1; i < n-1; i++) {
				for(int j = i+1; j <= n-1; j++) {
					//The isomorphism is rigged for max, so we need to insert a sign here
					d(i,j) = (-Addition::orientation())*matroidVector[index];
					index++;
				}
			}

			//Now apply mapping
			Vector<Rational> metric;
			for(int i = 1; i < n; i++) {
				for(int j = i+1; j <= n; j++) {
					if(j == n) {
						metric |= 0;
					}
					else {
						metric |= (2* d(i,j));
					}
				}
			}
			//dbgtrace << metric << endl;

			return curveFromMetric(metric); 
		}

	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see perl wrapper
	template <typename Addition>
		perl::ListReturn rational_curve_list_from_matroid_coordinates(Matrix<Rational> m) {
			perl::ListReturn result;

			for(int i = 0; i < m.rows(); i++) {
				result << rational_curve_from_matroid_coordinates<Addition>(m.row(i));
			}

			return result;
		}

	///////////////////////////////////////////////////////////////////////////////////////


}}

#endif
