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

	Implements the pullback of RationalFunction via a Morphism
	*/

#ifndef POLYMAKE_ATINT_PULLBACK_H
#define POLYMAKE_ATINT_PULLBACK_H

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/TropicalNumber.h"
#include "polymake/Polynomial.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/morphism_composition.h"

namespace polymake { namespace tropical {

	/*
	 * @brief Computes the pull-back form of a rational function given by numerator and denominator via
	 * a morphism given by a matrix and translate.
	 * @return A pair of (numerator, denominator).
	 */
	template <typename Addition>
		std::pair<Polynomial<TropicalNumber<Addition> >, Polynomial<TropicalNumber<Addition> > >
		polynomialPullback(const Matrix<Rational> &matrix, const Vector<Rational> &translate,
								const Polynomial<TropicalNumber<Addition> > &numerator, 
								const Polynomial<TropicalNumber<Addition> > &denominator) {
			Matrix<Rational> numerator_monoms( numerator.monomials_as_matrix());
			Vector<Rational > numerator_coeffs(numerator.coefficients_as_vector());
			Matrix<Rational> denominator_monoms( denominator.monomials_as_matrix());
			Vector<Rational> denominator_coeffs(denominator.coefficients_as_vector());

			Matrix<Rational> newnum_monoms = numerator_monoms*matrix;
			Matrix<Rational> newden_monoms = denominator_monoms*matrix;

			Vector<Rational> newnum_coeffs_rational = numerator_monoms*translate + numerator_coeffs;
			Vector<Rational> newden_coeffs_rational = denominator_monoms*translate + denominator_coeffs;


			Vector<TropicalNumber<Addition> > newnum_coeffs(newnum_coeffs_rational.dim());
			for(int i = 0; i < newnum_coeffs.dim(); i++) newnum_coeffs[i] = TropicalNumber<Addition>(newnum_coeffs_rational[i]);
			Vector<TropicalNumber<Addition> > newden_coeffs(newden_coeffs_rational.dim());
			for(int i = 0; i < newden_coeffs_rational.dim(); i++) newden_coeffs[i] = TropicalNumber<Addition>(newden_coeffs_rational[i]);

						
			Ring<TropicalNumber<Addition> > newring(matrix.cols());
			Polynomial<TropicalNumber<Addition> > newnum( (Matrix<int>(newnum_monoms)), newnum_coeffs, newring);
			Polynomial<TropicalNumber<Addition> > newden( (Matrix<int>(newden_monoms)), newden_coeffs, newring);
			return std::make_pair(newnum, newden);
		}

	template <typename Addition>
		perl::Object pullback(perl::Object morphism, perl::Object function) {
			//Convert the rational function to a morphism object
			perl::Object fmorphism(perl::ObjectType::construct<Addition>("Morphism"));
			perl::Object function_domain = function.give("DOMAIN");

			bool is_function_global = function.give("IS_GLOBALLY_DEFINED");
			bool is_morphism_global = morphism.give("IS_GLOBALLY_AFFINE_LINEAR");
			//If both are global, we're quickly done
			if(is_function_global && is_morphism_global) {
				Matrix<Rational> matrix = morphism.give("MATRIX");
				Vector<Rational> translate = morphism.give("TRANSLATE");
				Polynomial<TropicalNumber<Addition> > numerator = function.give("NUMERATOR");
				Polynomial<TropicalNumber<Addition> > denominator = function.give("DENOMINATOR");
				perl::Object newfunction(perl::ObjectType::construct<Addition>("RationalFunction"));
				std::pair<Polynomial<TropicalNumber<Addition> >, Polynomial<TropicalNumber<Addition> > > newpair =
					polynomialPullback(matrix, translate, numerator, denominator);
				newfunction.take("NUMERATOR") << newpair.first;
				newfunction.take("DENOMINATOR") << newpair.second;
				return newfunction;
			}

			Vector<Rational> function_vvalues = function.give("VERTEX_VALUES");
			Vector<Rational> function_lvalues = function.give("LINEALITY_VALUES");
			Matrix<Rational> fmorph_vvalues(function_vvalues.dim(),0);
			Matrix<Rational> fmorph_lvalues(function_lvalues.dim(),0);
				fmorph_vvalues |= function_vvalues;
				fmorph_lvalues |= function_lvalues;
			fmorphism.take("DOMAIN") << function_domain;
			fmorphism.take("VERTEX_VALUES") << thomog(fmorph_vvalues,0,false);
			fmorphism.take("LINEALITY_VALUES") << thomog(fmorph_lvalues,0,false);

			//Compute the composition
			perl::Object comp = morphism_composition<Addition>(morphism,fmorphism);

			//Now convert back to rational function
			perl::Object resultDomain = comp.give("DOMAIN");
			Matrix<Rational> result_vvalues = comp.give("VERTEX_VALUES");
			Matrix<Rational> result_lvalues = comp.give("LINEALITY_VALUES");

			perl::Object result(perl::ObjectType::construct<Addition>("RationalFunction"));
				result.take("DOMAIN") << resultDomain;
				result.take("VERTEX_VALUES") << tdehomog(result_vvalues,0,false).col(0);
				result.take("LINEALITY_VALUES") << (result_vvalues.rows() > 0? tdehomog(result_lvalues,0,false).col(0) : Vector<Rational>());
			if( (function.exists("NUMERATOR") || function.exists("DENOMINATOR")) &&
				 (morphism.exists("MATRIX") || morphism.exists("TRANSLATE")) ) {
				Matrix<Rational> matrix = morphism.give("MATRIX");
				Vector<Rational> translate = morphism.give("TRANSLATE");
				Polynomial<TropicalNumber<Addition> > numerator = function.give("NUMERATOR");
				Polynomial<TropicalNumber<Addition> > denominator= function.give("DENOMINATOR");
				perl::Object newfunction(perl::ObjectType::construct<Addition>("RationalFunction"));
				std::pair<Polynomial<TropicalNumber<Addition> >, Polynomial<TropicalNumber<Addition> > > newpair =
					polynomialPullback(matrix, translate, numerator, denominator);
				result.take("NUMERATOR") << newpair.first;
				result.take("DENOMINATOR") << newpair.second;
			}
			return result;
		}



}}

#endif
