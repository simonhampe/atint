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

	Implements the addition of two morphisms.
	*/

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/refine.h"

namespace polymake { namespace tropical {

	using namespace atintlog::donotlog;
	//   using namespace atintlog::dolog;
	//using namespace atintlog::dotrace;
	
	/*
	 * @brief Takes a morphism and a Cycle whose support is equal to f's [[DOMAIN]]
	 * and returns f as a morphism on that cycle.
	 * @param Morphism f
	 * @param Cycle X
	 * @param bool need_to_refine Whether the cycle needs to be refined or it is
	 * already fine in f's [[DOMAIN]].
	 * @tparam Addition Min or Max
	 * @return Morphism The "refined" morphism
	 */
	template <typename Addition>
		perl::Object refine_domain(perl::Object f, perl::Object cycle, bool need_to_refine) {
			perl::Object oldDomain = f.give("DOMAIN");
			bool has_matrix = f.exists("MATRIX");
			RefinementResult r = refinement(cycle, oldDomain, !has_matrix,false,false,need_to_refine,false);

			perl::Object nDomain = r.complex;
			if(has_matrix) {
				Matrix<Rational> matrix = f.give("MATRIX");
				Vector<Rational> translate = f.give("TRANSLATE");
				perl::Object result(perl::ObjectType::construct<Addition>("Morphism"));
					result.take("DOMAIN") << nDomain;
					result.take("MATRIX") << matrix;
					result.take("TRANSLATE") << translate;
				return result;
			}

			Matrix<Rational> rayRep = r.rayRepFromX;
			Matrix<Rational> linRep = r.linRepFromX;
		
			Matrix<Rational> vertex_values = f.give("VERTEX_VALUES");
			Matrix<Rational> lineality_values = f.give("LINEALITY_VALUES");
			Matrix<Rational> total_values = T( vertex_values / lineality_values);
			int target_dim = std::max(vertex_values.cols(), lineality_values.cols());

			Matrix<Rational> ndom_vertices = nDomain.give("SEPARATED_VERTICES");
			Matrix<Rational> ndom_lineality = nDomain.give("LINEALITY_SPACE");

			Matrix<Rational> rValues(0,target_dim);
			Matrix<Rational> lValues(0,target_dim);

			for(int r = 0; r < ndom_vertices.rows(); r++) {
				rValues /= total_values * rayRep.row(r);
			}

			for(int l = 0; l < ndom_lineality.rows(); l++) {
				lValues /= total_values * linRep.row(l);
			}

			perl::Object result(perl::ObjectType::construct<Addition>("Morphism"));
				result.take("DOMAIN") << nDomain;
				result.take("VERTEX_VALUES") << rValues;
				result.take("LINEALITY_VALUES") << lValues;
			return result;
		}

	/**
	  @brief Computes the sum of two morphisms (which should be defined on the same support)
	  @param perl::Object f A Morphism object
	  @param perl::Object g A Morphism object, whose DOMAIN has the same support as f's DOMAIN and whose image has the same ambient dimension as f's image
	  @return perl::Object A Morphism object representing f+g
	  */
	template <typename Addition>
	perl::Object add_morphisms(perl::Object f, perl::Object g) {
		//First we treat the special case where both are global
		bool f_global = f.exists("MATRIX");
		bool g_global = g.exists("MATRIX");
		Matrix<Rational> sum_matrix;
		Vector<Rational> sum_translate;
		if(f_global && g_global) {
			Matrix<Rational> fmatrix = f.give("MATRIX");
			Vector<Rational> ftranslate = f.give("TRANSLATE");
			Matrix<Rational> gmatrix = g.give("MATRIX");
			Vector<Rational> gtranslate = g.give("TRANSLATE");

			sum_matrix = fmatrix + gmatrix;
			sum_translate = ftranslate + gtranslate;
		}

		//dbgtrace << "Homogenizing where necessary" << endl;

		//First we homogenize where necessary
		perl::Object fDomain = f.give("DOMAIN");
		perl::Object gDomain = g.give("DOMAIN");
		//dbgtrace << "Computing refinement " << endl;

		//Then compute the common refinement of the domains
		RefinementResult r = refinement(fDomain,gDomain,false,false,false,true,false);
		perl::Object nDomain = r.complex;
		
		//If the map is given by a matrix, we're done
		if(f_global && g_global) {
			perl::Object result(perl::ObjectType::construct<Addition>("Morphism"));
				result.take("MATRIX") << sum_matrix;
				result.take("TRANSLATE") << sum_translate;
				result.take("DOMAIN") << nDomain;
			return result;
		}

		//Otherwise, compute values on the new domain
		perl::Object f_refined = refine_domain<Addition>(f, nDomain, false);
		perl::Object g_refined = refine_domain<Addition>(g, nDomain, false);

		Matrix<Rational> fref_vert = f_refined.give("VERTEX_VALUES");
		Matrix<Rational> gref_vert = g_refined.give("VERTEX_VALUES");
		Matrix<Rational> fref_lin  = f_refined.give("LINEALITY_VALUES");
		Matrix<Rational> gref_lin  = g_refined.give("LINEALITY_VALUES");
		
		perl::Object result(perl::ObjectType::construct<Addition>("Morphism"));
		result.take("DOMAIN") << nDomain;
		result.take("VERTEX_VALUES") << fref_vert + gref_vert;
		result.take("LINEALITY_VALUES") << fref_lin + gref_lin;

		return result;
	}



	UserFunctionTemplate4perl("# @category Morphisms"
			"# Computes the sum of two morphisms. Both [[DOMAIN]]s should have the same support"
			"# and the target spaces should have the same ambient dimension"
			"# The domain of the result will be the common refinement of the two domains."
			"# @param Morphism f"
			"# @param Morphism g"
			"# @tparam Addition Min or Max, both morphisms should use the same one."
			"# @return Morphism",
			"add_morphisms<Addition>(Morphism<Addition>, Morphism<Addition>)");
}}

