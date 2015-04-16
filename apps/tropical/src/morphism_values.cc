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

	Computes basic properties of Morphism
	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {


	/*
	 * @brief Computes the [[DOMAIN]] as the projective torus of right dimension
	 * from a given [[MATRIX]].
	 */
	template <typename Addition>
		void computeDomainFromMatrix(perl::Object morphism) {
			Matrix<Rational> mat = morphism.give("MATRIX");
			morphism.take("DOMAIN") << projective_torus<Addition>(mat.cols()-1,1);
		}

	/*
	 * @brief Computes [[VERTEX_VALUES]] and [[LINEALITY_VALUES]] from
	 * [[DOMAIN]], [[MATRIX]] and [[TRANSLATE]].
	 */
	void computeValuesFromMatrix(perl::Object morphism) {
		//Extract values
		perl::Object domain = morphism.give("DOMAIN");
		Matrix<Rational> rays = domain.give("VERTICES");
		Matrix<Rational> lineality = domain.give("LINEALITY_SPACE");
		Matrix<Rational> matrix = morphism.give("MATRIX");
		Vector<Rational> translate = morphism.give("TRANSLATE");

		Matrix<Rational> vertex_values = T(matrix * T(rays.minor(All,~scalar2set(0))));
		Matrix<Rational> lineality_values = T(matrix * T(lineality.minor(All,~scalar2set(0))));

		//For each nonfar vertex, we have to add the translate
		for(int r = 0; r < rays.rows(); r++) {
			if(rays(r,0) != 0) 
				vertex_values.row(r) += translate;
		}

		morphism.take("VERTEX_VALUES") << vertex_values;
		morphism.take("LINEALITY_VALUES") << lineality_values;
	}


	FunctionTemplate4perl("computeDomainFromMatrix<Addition>(Morphism<Addition>) : void");
	Function4perl(&computeValuesFromMatrix, "computeValuesFromMatrix(Morphism) : void");

}}
