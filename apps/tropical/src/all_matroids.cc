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

	Contains functions to compute the affine transform of a cycle 
	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/PowerSet.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/thomog.h"
#include "polymake/matroid/check_axioms.h"
#include "polymake/Ring.h"
#include "polymake/Polynomial.h"

namespace polymake { namespace tropical {

	/*
	 * Computes all loopfree matroids of given rank on 0 .. n-1
	 */
	perl::ListReturn all_loopfree_of_rank(int n, int r) {
		perl::ListReturn result;

		Vector<Set<int> > r_sets(all_subsets_of_k(sequence(0,n),r));
		Array<Set<int> > all_sets = all_subsets(sequence(0, r_sets.dim()));
		//Iterate all possible subsets of r-sets
		for(Entire<Array< Set<int> > >::iterator as = entire(all_sets); !as.at_end(); as++) {
			Set<int> base_union = accumulate(r_sets.slice( *as), operations::add());
			if(base_union.size() == n) {
				if(matroid::check_basis_exchange_axiom_impl( r_sets.slice( *as))) {
					perl::Object m("matroid::Matroid");
					m.take("N_ELEMENTS") << n;
					m.take("BASES") << r_sets.slice(*as);
					result << m;
				}
			}
		}

		return result;
	}

	/*
	 * Computes all (not only non-isomorphic) matroids on 0 .. n-1
	 */
	perl::ListReturn all_loopfree_matroids_on_n(int n) {
		perl::ListReturn result;
	
		//Iterate ranks

		for(int r = 1; r <= n; r++) {
			cout << "Finding rank " << r << " matroids." << endl;
			Vector<Set<int> > r_sets(all_subsets_of_k(sequence(0,n),r));
			Array<Set<int> > all_sets = all_subsets(sequence(0, r_sets.dim()));
			//Iterate all possible subsets of r-sets
			for(Entire<Array< Set<int> > >::iterator as = entire(all_sets); !as.at_end(); as++) {
				Set<int> base_union = accumulate(r_sets.slice( *as), operations::add());
				if(base_union.size() == n) {
					if(matroid::check_basis_exchange_axiom_impl( r_sets.slice( *as))) {
						perl::Object m("matroid::Matroid");
							m.take("N_ELEMENTS") << n;
							m.take("BASES") << r_sets.slice(*as);
						result << m;
					}
				}
			}
		}
		

		return result;
	}

	//Takes a list of matroids m_1,...,m_k closed under intersection (0 is allowed)
	//and computes the ideal I s.t. Z[m_1,...,m_k] = Z[x_1,...,x_k]/I
	perl::Object matroid_ring_ideal(Array<perl::Object> matroids) {
		Ring<Rational> r(matroids.size());
		Array<Polynomial<Rational> > vars = r.variables();
		Vector<Polynomial<Rational> > vars_vector(vars);

		Map< Array<Set<int> >, int> index_map;
		for(int i = 0; i < matroids.size(); i++) {
			Array<Set<int> > bases = matroids[i].give("BASES");
			index_map[bases] = i;
		}
		index_map[ Array<Set<int> >()] = matroids.size();

		Array<Set<int> > pairs = all_subsets_of_k( sequence(0, matroids.size()), 2);
		Array<Polynomial<Rational> > relations(pairs.size());

		int pindex = 0;
		for(Entire<Array<Set<int> > >::iterator p = entire(pairs); !p.at_end(); p++, pindex++) {
			cout << "Computing " << pindex+1 << " of " << pairs.size() << endl;
			Vector<int> pvector( *p);
			perl::Object fan_1 = CallPolymakeFunction("matroid_fan_max", matroids[ pvector[0]]);
			perl::Object fan_2 = CallPolymakeFunction("matroid_fan_max", matroids[ pvector[1]]);
			perl::Object inter = CallPolymakeFunction("intersect", fan_1, fan_2);
			if(CallPolymakeFunction("is_empty",inter)) {
				relations[pindex] = vars_vector[ pvector[0]] * vars_vector[ pvector[1]];
			}
			else {
				perl::Object imatroid = CallPolymakeFunction("matroid_from_fan",inter);
				Array<Set<int> > ibases = imatroid.give("BASES");
				int iindex = index_map[ibases];
				relations[pindex] = vars_vector[ pvector[0]] * vars_vector[ pvector[1]] - vars_vector[iindex];
			}
		}

		perl::Object ideal("ideal::Ideal");
			ideal.take("GENERATORS") << relations;
		return ideal;
	}

	UserFunction4perl("", &all_loopfree_of_rank,"all_loopfree_of_rank($,$) : returns(@)");

	UserFunction4perl("", &all_loopfree_matroids_on_n, "all_loopfree_matroids_on_n($) : returns(@)");

	UserFunction4perl("", &matroid_ring_ideal, "matroid_ring_ideal(matroid::Matroid+)");

}}
