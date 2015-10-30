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
#include "polymake/PowerSet.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/TropicalNumber.h"
#include "polymake/Ring.h"
#include "polymake/Polynomial.h"

namespace polymake { namespace tropical {

	/*
	 * Takes a tropical polynomial F in n variables and a tuple of n tropical polynomials g_1,..,g_n living
	 * in the same ring.
	 * It then computes for each monomial m of F the evaluation m(g_1,...,g_n)
	 */
	template <typename Addition>
		Array<Polynomial<TropicalNumber<Addition> > > evaluate_monomials(Polynomial<TropicalNumber<Addition> > poly, Array<Polynomial<TropicalNumber<Addition> > > tuple) {
			typedef TropicalNumber<Addition> TNumber;
			Matrix<int> monomials = poly.monomials_as_matrix();
			Vector<TNumber> coeffs = poly.coefficients_as_vector();
			Array<Polynomial<TNumber> > result(monomials.rows());
			
			Ring<TNumber> tuple_ring = tuple[0].get_ring();

			for(int i = 0; i < monomials.rows(); i++) {
				Polynomial<TNumber> ev_poly(coeffs[i],tuple_ring);
				for(int j = 0; j< monomials.cols(); j++) {
					for(int p = 0; p < monomials(i,j); p++) {
						ev_poly *= tuple[j];
					}
				}
				result[i] = ev_poly;
			}

			return result;
		}

	//Only keeps monomials not using a given list of variables.
	template <typename Addition>
		Polynomial<TropicalNumber<Addition> > forget_coordinates(Polynomial<TropicalNumber<Addition> > poly, Set<int> coords) {
			typedef TropicalNumber<Addition> TNumber;
			Ring<TNumber> r(poly.n_vars()-coords.size());
			Set<int> using_monomials;
			for(Entire<Set<int> >::iterator c = entire(coords); !c.at_end(); c++) {
				using_monomials += support(poly.monomials_as_matrix().col(*c));
			}
			return Polynomial<TNumber>(poly.monomials_as_matrix().minor(~using_monomials,~coords), poly.coefficients_as_vector().slice(~using_monomials), r);
		}

	UserFunctionTemplate4perl("","evaluate_monomials<Addition>(Polynomial<TropicalNumber<Addition>>, Polynomial<TropicalNumber<Addition>>+)");

	UserFunctionTemplate4perl("", "forget_coordinates<Addition>(Polynomial<TropicalNumber<Addition>>,$)");

}}

