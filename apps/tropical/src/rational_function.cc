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

	Functions for the basic ruleset of RationalFunction
	*/

#include "polymake/client.h"
#include "polymake/Set.h"
#include "polymake/Array.h"
#include "polymake/Matrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Polynomial.h"
#include "polymake/TropicalNumber.h"
#include "polymake/tropical/specialcycles.h"
#include "polymake/tropical/refine.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/linear_algebra_tools.h"
#include "polymake/tropical/polynomial_tools.h"
#include "polymake/tropical/LoggingPrinter.h"


namespace polymake { namespace tropical {

	typedef std::pair<Vector<Rational>, Vector<Rational> > ValuePair;

	template <typename Addition>
		perl::Object computePolynomialDomain(const Polynomial<TropicalNumber<Addition> > &p) {
			Matrix<Rational> monoms(p.monomials_as_matrix());
			Vector<TropicalNumber<Addition> > coefs = p.coefficients_as_vector();

			if(monoms.rows() <= 1) {
				return projective_torus<Addition>(monoms.cols()-1,0);
			}

			//FIXME Same computation as in the beginning of the hypersurface client. Refactor?
			ListMatrix< Vector<Rational> > ineq;
			const TropicalNumber<Addition> zero=TropicalNumber<Addition>::zero();
			for (int i=0; i< monoms.rows(); ++i) {
				if (coefs[i]==zero)
					ineq /= unit_vector<Rational>(monoms.cols()+1,0);
				else
					ineq /= Addition::orientation()*(Rational(coefs[i])|monoms[i]);
			}

			perl::Object dome(perl::ObjectType::construct<Rational>("polytope::Polytope"));
			dome.take("INEQUALITIES") << ineq;
			dome.take("FEASIBLE") << true;
			dome.take("BOUNDED") << false;

			Matrix<Rational> vertices = dome.give("VERTICES");
			Matrix<Rational> lineality = dome.give("LINEALITY_SPACE");
			IncidenceMatrix<> polytopes = dome.give("VERTICES_IN_FACETS");
			Set<int> far_face = dome.give("FAR_FACE");

			//Find and eliminate the far face - Start from the end since it seems to be always there.
			for(int r = polytopes.rows()-1; r >=0; r++) {
				if(polytopes.row(r) == far_face) {
					polytopes = polytopes.minor(~scalar2set(r),All);
					break;
				}
			}
		
			perl::Object domain(perl::ObjectType::construct<Addition>("Cycle"));
				domain.take("VERTICES") << vertices;
				domain.take("MAXIMAL_POLYTOPES") << polytopes;
				domain.take("LINEALITY_SPACE") << lineality;
			return domain;
		}//END computePolynomialDomain

	/**
	 * @brief Computes the property [[DOMAIN]] from [[NUMERATOR]] and [[DENOMINATOR]].
	 */
	template <typename Addition>
		void computeDomain(perl::Object function) {
			Polynomial<TropicalNumber<Addition> > num = function.give("NUMERATOR");
			Polynomial<TropicalNumber<Addition> > den = function.give("DENOMINATOR");
			
			perl::Object domain_num = computePolynomialDomain(num);
			perl::Object domain_den = computePolynomialDomain(den);

			RefinementResult r = refinement(domain_num, domain_den, false,false,false,true,false);
			function.take("DOMAIN") << r.complex;
		}//END computeDomain


	/**
	 * @brief Computes properties [[VERTEX_VALUES]] and [[LINEALITY_VALUES]]
	 * from [[[DOMAIN]], [NUMERATOR]] and [[DENOMINATOR]].
	 */
	template <typename Addition>
		void computeGeometricFunctionData(perl::Object function) {
			Polynomial<TropicalNumber<Addition> > num = function.give("NUMERATOR");
			Polynomial<TropicalNumber<Addition> > den = function.give("DENOMINATOR");
			perl::Object domain = function.give("DOMAIN");
		
			//This just computes the associated vertices
			RefinementResult r = refinement(domain, domain, false,false,true,false,false);

			Matrix<Rational> separated_vertices = domain.give("SEPARATED_VERTICES");
				std::pair<Set<int>,Set<int> > vertex_list = far_and_nonfar_vertices(separated_vertices);
				separated_vertices = separated_vertices.minor(All,~scalar2set(0));
			Matrix<Rational> lineality = domain.give("LINEALITY_SPACE");
				lineality = lineality.minor(All,~scalar2set(0));
			Vector<int> assocRep = r.associatedRep;

			Vector<Rational> vertexValues(separated_vertices.rows());
			Vector<Rational> linealityValues(lineality.rows());


			//Compute values for all nonfar vertices
			for(Entire<Set<int> >::iterator v = entire(vertex_list.second); !v.at_end(); v++) {
				vertexValues[*v] = evaluate_polynomial(num, separated_vertices.row(*v)) -
										 evaluate_polynomial(den, separated_vertices.row(*v));
			}

			//For all far vertices, compute slope with respect to associated vertex
			for(Entire<Set<int> >::iterator r= entire(vertex_list.first); !r.at_end(); r++) {
				Vector<Rational> associated_vertex = separated_vertices.row(assocRep[*r]);
				vertexValues[*r] = ( evaluate_polynomial(num, associated_vertex + separated_vertices.row(*r)) -
											evaluate_polynomial(den, associated_vertex + separated_vertices.row(*r))) -
											vertexValues[assocRep[*r]];
			}

			function.take("VERTEX_VALUES") << vertexValues;

			//Same for lineality space generators - we use a fixed vertex as base vertex
			int base_vertex_index = *(vertex_list.second.begin());
			Vector<Rational> base_vertex = separated_vertices.row(base_vertex_index);
			for(int l = 0; l < lineality.rows(); l++) {
				linealityValues[l] = ( evaluate_polynomial(num, base_vertex + lineality.row(l)) -
											  evaluate_polynomial(den, base_vertex + lineality.row(l)) ) -
											vertexValues[base_vertex_index];
			}

			function.take("LINEALITY_VALUES") << linealityValues;
		}//END computeGeometricFunctionData

	/*
	 * @brief Computes the restriction of a RationalFunction to a cycle contained in its [[DOMAIN]].
	 */
	template <typename Addition>
		perl::Object restrict(perl::Object function, perl::Object cycle, bool compute_lattice) {
			//First check if the function is given by a polynomial quotient
			if(function.exists("NUMERATOR") && function.exists("DENOMINATOR")) {
				Polynomial<TropicalNumber<Addition> > num = function.give("NUMERATOR");
				Polynomial<TropicalNumber<Addition> > den = function.give("DENOMINATOR");
				perl::Object domain = function.give("DOMAIN");
				RefinementResult ref_cycle = refinement(cycle,domain, false,false,false,true,compute_lattice);
				perl::Object result(perl::ObjectType::construct<Addition>("RationalFunction"));
					result.take("DOMAIN") << ref_cycle.complex;
					result.take("NUMERATOR") << num;
					result.take("DENOMINATOR") << den;
				return result;
			}

			//Otherwise, refine the cycle along the domain and compute representations in terms
			//of the domain
			perl::Object domain = function.give("DOMAIN");
			Vector<Rational> vertex_values = function.give("VERTEX_VALUES");
			Vector<Rational> lineality_values = function.give("LINEALITY_VALUES");
			Vector<Rational> full_values = vertex_values | lineality_values;
			RefinementResult ref_cycle = refinement(cycle, domain, false,true,false,true,compute_lattice);

			Matrix<Rational> vertexRep = ref_cycle.rayRepFromY;
			Matrix<Rational> linealityRep = ref_cycle.linRepFromY;

			Vector<Rational> new_vertex_rep = vertexRep * full_values;
			Vector<Rational> new_lineality_rep = linealityRep * full_values;

			perl::Object result(perl::ObjectType::construct<Addition>("RationalFunction"));
				result.take("DOMAIN") << ref_cycle.complex;
				result.take("VERTEX_VALUES") << new_vertex_rep;
				result.take("LINEALITY_VALUES") << new_lineality_rep;

			return result;
		}

	// PERL WRAPPER //////////////////////////

	FunctionTemplate4perl("computePolynomialDomain<Addition>(Polynomial<TropicalNumber<Addition> >)");
	FunctionTemplate4perl("computeDomain<Addition>(RationalFunction<Addition>) : void");
	FunctionTemplate4perl("computeGeometricFunctionData<Addition>(RationalFunction<Addition>) : void");
	FunctionTemplate4perl("restrict<Addition>(RationalFunction<Addition>, Cycle<Addition>;$=0)");


}}
