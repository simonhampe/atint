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

	//Takes a list of matroid_fans and computes all possible intersections and the relations between them.
	perl::Object matroid_intersection_ideal(Array<perl::Object> matroid_fans) {
		Ring<Rational> r(matroid_fans.size());
		Array<Polynomial<Rational> > vars = r.variables();
		Vector<Polynomial<Rational> > vars_vector(vars);
		Vector<Polynomial<Rational> > relations;

		perl::Object first_matroid = matroid_fans[0];
		int n = first_matroid.give("PROJECTIVE_AMBIENT_DIM");

		//Maps lists of indices to the corresponding intersection product
		Vector<Vector<int> > current_tuples;
		Map< Vector<int>, perl::Object> intersection_cycles;
		Map< int, Polynomial<Rational> > var_map;
		for(int i = 0; i < matroid_fans.size(); i++) {
			Vector<int> atom(scalar2set(i));
			current_tuples |= atom;
			intersection_cycles[atom] = matroid_fans[i]; 
			var_map[i] = vars_vector[i];
		}

		cout << current_tuples << endl;

		for(int j = 2; j <= n; j++) {
			cout << "Codimension " << j-1 << endl;
			Vector<Vector<int> > next_tuples;
			for(int ct = 0; ct < current_tuples.dim(); ct++) {
				for(int next_index = current_tuples[ct][current_tuples[ct].dim()-1]; next_index < matroid_fans.size(); next_index++) {
					next_tuples |= (current_tuples[ct] | next_index);
					perl::Object isection = CallPolymakeFunction("intersect", intersection_cycles[ current_tuples[ct]], matroid_fans[next_index]);
					intersection_cycles[ next_tuples[next_tuples.dim()-1] ] = isection;
				}
			}
			//Check equality of intersection products
			for(int p1 = 0; p1 < next_tuples.dim(); p1++) {
				if(CallPolymakeFunction("is_empty", intersection_cycles[next_tuples[p1]])) {
					relations  |=  accumulate(  attach_operation( next_tuples[p1],pm::operations::associative_access<Map<int, Polynomial<Rational> >, int >(&var_map)), operations::mul()  );
				}
				for(int p2 = p1+1; p2 < next_tuples.dim(); p2++) {
					if(CallPolymakeFunction("check_cycle_equality", intersection_cycles[next_tuples[p1]], intersection_cycles[next_tuples[p2]])) {
						relations |= ( accumulate(  attach_operation( next_tuples[p1],pm::operations::associative_access<Map<int, Polynomial<Rational> >, int >(&var_map)), operations::mul()  ) -
											accumulate(  attach_operation( next_tuples[p2],pm::operations::associative_access<Map<int, Polynomial<Rational> >, int >(&var_map)), operations::mul()  ));
					}
				}
			}
			cout << next_tuples << endl;
			current_tuples = next_tuples;
		}


		perl::Object result("ideal::Ideal");
			result.take("GENERATORS") << relations;
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

	perl::Object matroid_intersection(perl::Object m1, perl::Object m2) {
		Array<Set<int> > b1 = m1.give("BASES");
		Array<Set<int> > b2 = m2.give("BASES");
		int r1 = m1.give("RANK");
		int r2 = m2.give("RANK");
		int n = m1.give("N_ELEMENTS");

		int inter_rank = r1 + r2 - n;

		Set<Set<int> > newbases;
		for(int i1 = 0; i1 < b1.size(); i1++) {
			for(int i2 = 0; i2 < b2.size(); i2++) {
				Set<int> c = b1[i1] * b2[i2];
				if(c.size() == inter_rank) newbases += c;
			}
		}

		perl::Object result("matroid::Matroid");
			result.take("N_ELEMENTS") << n;
			result.take("BASES") << Array<Set<int> >(newbases);
		return result;
	}

	UserFunction4perl("", &all_loopfree_of_rank,"all_loopfree_of_rank($,$) : returns(@)");

	UserFunction4perl("", &all_loopfree_matroids_on_n, "all_loopfree_matroids_on_n($) : returns(@)");

	UserFunction4perl("", &matroid_ring_ideal, "matroid_ring_ideal(matroid::Matroid+)");

	UserFunction4perl("", &matroid_intersection_ideal, "matroid_intersection_ideal(Cycle+)");

	UserFunction4perl("", &matroid_intersection, "matroid_intersection(matroid::Matroid, matroid::Matroid)");

}}
