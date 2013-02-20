/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 * 
 * ---
 * Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 * 
 * 
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/fan_piecewise_divisor.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
//   //using namespace atintlog::dolog;
// //   using namespace atintlog::dotrace;
//   
//   //Returns all possible vectors of length "length", filled with elements in "content"
//   Matrix<int> fill_value_vector(int length, Vector<int> content) {
//     Matrix<int> result(0,length);
//     if(content.dim() == 1) {
//       result /= (content[0] * ones_vector<int>(length));
//       return result;
//     }
//     //Go through all subsets
//     Array<Set<int> > subsets = all_subsets(sequence(0,length));
//     for(int s = 0; s < subsets.size(); s++) {
//       //We want each entry at least once, so we leave enough space
// //       if(length - subsets[s].size() >= content.dim() -1 && subsets[s].size() > 0) {
// 	//Fill subset s with first entry
// 	Vector<int> v(length);
// 	v.slice(subsets[s]) = content[0] * ones_vector<int>(subsets[s].size());
// 	//Continue recursively:
// 
// 	Matrix<int> recursive = 
// 	  fill_value_vector(length - subsets[s].size(), content.slice(~scalar2set(0)));
// 	for(int r = 0; r < recursive.rows(); r++) {
// 	  v.slice(~subsets[s]) = recursive.row(r);
// 	  result /= v;
// 	}
//       
// //       }
//     }
//     return result;
//   }
// 
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   Matrix<Rational> eq_matrix(perl::Object div) {
//     
//     Set<int> first_half = sequence(0,3);
//     Set<int> second_half = sequence(3,3);
//     //Extract values
//       Matrix<Rational> rays = div.give("RAYS");
//       IncidenceMatrix<> codim = div.give("CODIM_1_FACES");
//       IncidenceMatrix<> codimInMax = div.give("CODIM_1_IN_MAXIMAL_CONES");
//       Map<int, Map<int, Vector<Rational> > > lnFunctionVector = div.give("LATTICE_NORMAL_FCT_VECTOR");
//       Matrix<Rational> lsumFunctionVector = div.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
//       Vector<Integer> weights = div.give("TROPICAL_WEIGHTS");
//       
//       //Find the four diagonal rays (if they exist)
//       Set<int> diag_rays;
//       for(int r = 0; r < rays.rows(); r++) {
// 	if(rays.row(r).slice(first_half) == rays.row(r).slice(second_half)) diag_rays += r;
//       }
//       
//       //Find the six cones that consist of diagonal rays
//       Set<int> diag_cones;
//       for(int c = 0; c < codim.rows(); c++) {
// 	if((codim.row(c) * diag_rays).size() == 2) diag_cones += c;
//       }
//       
// //       dbgtrace << "Diag rays " << diag_rays << ", " << "Diag cones " << diag_cones << endl;
//       
//       //Now check if we can cut out the diagonal with the appropriate linear equation system
//       Matrix<Rational> eq(0,12);
//       for(int co = 0; co < codim.rows(); co++) {
// 	Vector<Rational> v = lsumFunctionVector.row(co);
// 	Set<int> adjacent = codimInMax.row(co);
// 	for(Entire<Set<int> >::iterator mc = entire(adjacent); !mc.at_end(); mc++) {
// 	    v += weights[*mc] * (lnFunctionVector[co])[*mc];
// 	}
// 	eq /= v;
//       }
// //       dbgtrace << "Matrix " << eq << endl;
//       Vector<Rational> desired = zero_vector<Rational>(eq.rows());
// 	desired.slice(diag_cones) = ones_vector<Rational>(6);
// 	
// 	return (eq | desired);
//   }
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   Matrix<Rational> test_diagonal_combinations(perl::Object D) {
//     
// //     Vector<Rational> first_index;
// //     Vector<Rational> second_index;
// //     Vector<Rational> third_index;
// //     for(int i = 0; i < 4; i++) {
// //       for(int j = i+1; j < 4; j++) {
// // 	for(int k = j+1; k < 4; k++) {
// // 	  first_index |= i;
// // 	  second_index |= j;
// // 	  third_index |= k;
// // 	}
// //       }
// //     }
//     
//     int translate = 4;
//     
//     pm::cout << "Creating possibilities" << endl;
//     Matrix<int> possibilities = fill_value_vector(12, Vector<int>(sequence(0,3)));
//     
//     
//     Matrix<Rational> result(0,12);
//     for(int p = 0; p < possibilities.rows(); p++) {
//       pm::cout << "Testing " << p << " of " << possibilities.rows() << ", " << result.rows() << " found " << endl;
//       
//       int i = 0;
//       Vector<Set<int> > piecewise;
//       Vector<Integer> signs;
//       for(int a = 0; a < 4; a++) {
// 	for(int b = a+1; b < 4; b++) {
// 	   for(int c = b+1; c < 4; c++) {
// 	    Set<int> cone;
// 	    cone += (a + translate * possibilities(p,i));
// 	    cone += (b + translate * possibilities(p,i+1));
// 	    cone += (c + translate * possibilities(p,i+2));
// 	    piecewise |= cone;
// 	    signs |= (possibilities(p,i) == 2? -1 : 1) * (possibilities(p,i+1) == 2? -1 : 1) * (possibilities(p,i+2) == 2? -1 : 1);
// 	    i+=3;
// 	   }
// 	}
//       }
//       
//       perl::Object div  = piecewise_divisor(D, piecewise, signs);
//       
//       //Extract values
//       Vector<Integer> weights = div.give("TROPICAL_WEIGHTS");
//       Matrix<Rational> rays = div.give("RAYS");
//       if(weights == ones_vector<Integer>(4)) {
// 	if(rays.minor(All,sequence(0,3)) == rays.minor(All,sequence(3,3))) {
// 	    result /= possibilities.row(p);
// 	    pm::cout << "WORKS" << endl;
// 	}
//       }
//       
//       
//       
//       
//     }//END go through possibilities
//     
//     
//     return result;
//   }
//   
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   Matrix<Rational> test_sign_combinations(perl::Object d, IncidenceMatrix<> cones) {
//     
//     Vector<int> signs; signs |= 1; signs |= -1;
//     Matrix<int> possibilities = fill_value_vector(10,signs);
//     
//     Matrix<Rational> result(0,10);
//     
//     for(int s = 0; s < possibilities.rows(); s++) {
//       pm::cout << "Testing " << s << " of " << possibilities.rows() << endl;
//       
//       perl::Object div  = piecewise_divisor(d, cones, Vector<Integer>(possibilities.row(s)));
//       
//       //Extract values
//       Vector<Integer> weights = div.give("TROPICAL_WEIGHTS");
//       Matrix<Rational> rays = div.give("RAYS");
//       if(weights == ones_vector<Integer>(6)) {
// 	if(rays.minor(All,sequence(0,3)) == rays.minor(All,sequence(3,3))) {
// 	    result /= possibilities.row(s);
// 	    pm::cout << "WORKS" << endl;
// 	}
//       }
//       
//     }
//     
//     return result;
//     
//   }

  
//   
  // ------------------------- PERL WRAPPERS ---------------------------------------------------

 
//   Function4perl(&eq_matrix, "l32em(WeightedComplex)");
//   Function4perl(&test_diagonal_combinations, "l32td(WeightedComplex)");
//   Function4perl(&test_sign_combinations, "l43ts(WeightedComplex,IncidenceMatrix)");
  
}}