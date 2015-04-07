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

	This file provides functionality to compute certain special tropical varieties
	*/

#include "polymake/client.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/PowerSet.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Array.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Graph.h"
#include "polymake/linalg.h"

namespace polymake { namespace tropical {

	//using namespace atintlog::donotlog;
	//using namespace atintlog::dolog;
	//using namespace atintlog::dotrace;

	template <typename Addition>
		perl::Object empty_cycle(int ambient_dim) {
			perl::Object cycle(perl::ObjectType::construct<Addition>("Cycle"));
				cycle.take("VERTICES") << Matrix<Rational>(0,ambient_dim);
				cycle.take("MAXIMAL_POLYTOPES") << Array<Set<int> >();
				cycle.take("WEIGHTS") << Vector<Integer>();
				cycle.take("PROJECTIVE_AMBIENT_DIM") << ambient_dim;
				cycle.set_description() << "Empty cycle in dimension " << ambient_dim;

			return cycle;
		}//END empty_cycle

	template <typename Addition>
		perl::Object point_collection(Matrix<Rational> m, Vector<Integer> weights) {
			//Sanity check
			if(m.rows() == 0)
				throw std::runtime_error("No points given.");
			if(m.rows() != weights.dim()) 
				throw std::runtime_error("Number of points does not match number of weights");
			
			//Create vertices
			m = ones_vector<Rational>(m.rows()) | m;

			//Create polytopes
			Array<Set<int> > polytopes(m.rows());
			for(int i = 0; i < polytopes.size(); i++)
				polytopes[i] = scalar2set(i);

			perl::Object cycle(perl::ObjectType::construct<Addition>("Cycle"));
				cycle.take("VERTICES") << m;
				cycle.take("MAXIMAL_POLYTOPES") << polytopes;
				cycle.take("WEIGHTS") << weights;
		

			return cycle;
		}//END point_collection


	template <typename Addition> 
		perl::Object uniform_linear_space(const int n, const int k) {

			//Ensure that dimensions match
			if(k > n) 
				throw std::runtime_error("Cannot create uniform linear space. Fan dimension is larger than ambient dimension.");
			if(k < 0 || n < 0) 
				throw std::runtime_error("Cannot create uniform linear space. Negative dimension provided.");
			if(k == 0) {
				return point_collection<Addition>( Matrix<Rational>(1,n+1), ones_vector<Integer>(1));	
			}

			//Create rays
			Matrix<Rational> vertices(unit_matrix<Rational>(n+1));
			vertices = zero_vector<Rational>(n+1) | vertices;
			vertices *= Addition::orientation();
			vertices = unit_vector<Rational>(n+2,0) / vertices;

			//Create cones
			Array<Set<int> > polytopes = all_subsets_of_k( sequence(1,n+1), k);
			for(int i = 0; i < polytopes.size(); i++) 
				polytopes[i] += 0;

			//Create weights
			Vector<Integer> weights = ones_vector<Integer>(polytopes.size());

			//Create final object
			perl::Object fan(perl::ObjectType::construct<Addition>("Cycle"));
			fan.take("VERTICES") << vertices;
			fan.take("MAXIMAL_POLYTOPES") << polytopes;
			fan.take("WEIGHTS") << weights;
			fan.set_description() << "Uniform linear space of dimension " << k << " in dimension " << n;
			return fan;

		}//END uniform_linear_space

	template <typename Addition>
		perl::Object halfspace_subdivision(Rational a, Vector<Rational> g, Integer weight) {
			//Sanity check
			if(g == zero_vector<Rational>(g.dim())) 
					throw std::runtime_error("Zero vector does not define a hyperplane.");
			if(g*ones_vector<Rational>(g.dim()) != 0) 
					throw std::runtime_error("Normal vector must be homogenous, i.e. sum of entries must be zero");

			//Create vertices
			Matrix<Rational> vertices(0,g.dim());
			vertices /= g;
			vertices /= (-g);
			vertices = zero_vector<Rational>(2) | vertices;

			//Create lineality
			Matrix<Rational> lineality = null_space(g).minor(~scalar2set(0),All);
			lineality = zero_vector<Rational>(lineality.rows()) | lineality;

			//Compute apex
			Rational sum = accumulate(attach_operation(g,operations::square()),operations::add());
			Vector<Rational> apex = (a/sum) *g;
				apex = Rational(1) | apex;
			vertices = apex / vertices;

			Array<Set<int> > polytopes(2);
				polytopes[0] = (sequence(0,3)-1);
				polytopes[1] = (sequence(0,3)-2);

			perl::Object cycle(perl::ObjectType::construct<Addition>("Cycle"));
				cycle.take("VERTICES") << vertices;
				cycle.take("MAXIMAL_POLYTOPES") << polytopes;
				if(lineality.rows() > 0)
					cycle.take("LINEALITY_SPACE") << lineality;
				cycle.take("WEIGHTS") << weight * ones_vector<Integer>(2);

			return cycle;
		}//END halfspace_subdivision

	// PERL WRAPPER -------------------------------------------

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
									"# Creates the empty cycle in a given ambient dimension"
									"# (i.e. it will set the property [[PROJECTIVE_AMBIENT_DIM]]."
									"# @param Int ambient_dim The ambient dimension"
									"# @tparam Addition Max or Min"
									"# @return Cycle The empty cycle",
									"empty_cycle<Addition>($)");
	
	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
							"# Creates a cycle consisting of a collection of points"
							"# with given weights"
							"# @param Matrix<Rational> points The points, in tropical homogeneous coordinates"
							"# (though not with leading ones for vertices)."
							"# @param Vector<Integer> weights The list of weights for the points"
							"# @tparam Addition Max or Min"
							"# @return Cycle The point collection.",
							"point_collection<Addition>(Matrix<Rational>, Vector<Integer>)");

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
			"# Creates the linear space of the uniform matroid of rank k+1 on n+1 variables."
			"# @param Int n The ambient (projective) dimension."
			"# @param Int k The (projective dimension of the fan."
			"# @tparam Addition A The tropical addition (min or max)"
			"# @return Cycle A tropical linear space.",
			"uniform_linear_space<Addition>($,$)");       

	UserFunctionTemplate4perl("# @category Creation functions for specific cycles"
									"# Creates a subdivision of the tropical projective torus"
									"# along an affine hyperplane into two halfspaces."
									"# This hyperplane is defined by an equation gx = a"
									"# @param Rational a The constant coefficient of the equation"
									"# @param Vector<Rational> g The linear coefficients of the equation"
									"# Note that the equation must be homogeneous in the sense that (1,..1)"
									"# is in its kernel, i.e. all entries of g add up to 0."
									"# @param Integer The (constant) weight this cycle should have"
									"# @tparam Addition Max or Min"
									"# @return Cycle The halfspace subdivision",
									"halfspace_subdivision<Addition>($,Vector<Rational>,$)");

}}

