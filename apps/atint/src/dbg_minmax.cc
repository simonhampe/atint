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
#include "polymake/IncidenceMatrix.h"
#include "polymake/Graph.h"
#include "polymake/linalg.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
//   //using namespace atintlog::dolog;
//   //using namespace atintlog::dotrace;
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   perl::Object computeMinMaxSpace(perl::Object function) {
//     perl::Object domain = function.give("DOMAIN");
//       Matrix<Rational> rays = domain.give("RAYS");
//       IncidenceMatrix<> cones = domain.give("MAXIMAL_CONES");
//       
//     Vector<Rational> values = function.give("RAY_VALUES");
//       
//       
//     int number_of_vars = 5*cones.rows();
//     Matrix<Rational> ineqs(0,number_of_vars+1);
//     Matrix<Rational> eqs(0, number_of_vars+1);
//     
//     //Create equalities:
//     for(int sigma = 0; sigma < cones.rows(); sigma++) {
//       Set<int> rays_in_sigma = cones.row(sigma);
//       for(Entire<Set<int> >::iterator r = entire(rays_in_sigma); !r.at_end(); r++) {
// 	Vector<Rational> r_sigma_eq(number_of_vars+1);
// 	r_sigma_eq[0] = values[*r];
// 	r_sigma_eq.slice(sequence(5*sigma +1,5)) = -rays.row(*r);
// 	eqs /= r_sigma_eq;
//       }
//     }
//     
//     //Create inequalities
//     IncidenceMatrix<> cones_of_rays = T(cones);
//     for(int r = 0; r < rays.rows(); r++) {
//       Set<int> non_adjacent = sequence(0,cones.rows()) - cones_of_rays.row(r);
//       for(Entire<Set<int> >::iterator sigma = entire(non_adjacent); !sigma.at_end(); sigma++) {
// 	Vector<Rational> r_sigma_ineq(number_of_vars+1);
// 	r_sigma_ineq[0] = values[r];
// 	r_sigma_ineq.slice(sequence(5* (*sigma) + 1,5)) = - rays.row(r);
// 	ineqs /= r_sigma_ineq;
//       }
//     }
//     
//     
//     
//     perl::Object ptope("polytope::Polytope");
//       ptope.take("INEQUALITIES") << ineqs;
//       ptope.take("EQUATIONS") << eqs;
//       
//     return ptope;
//   }
//   

     IncidenceMatrix<> doflip(IncidenceMatrix<> cells, Set<int> tcone1, Set<int> tcone2) {
	
	//Compute complementary subdivision
	Set<int> inter = tcone1 * tcone2;
	Set<int> symdif = (tcone1 + tcone2) - inter;
	Vector<Set<int> > triang_simpl;
	    for(Entire<Set<int> >::iterator i = entire(inter);!i.at_end(); i++) {
		triang_simpl |= (symdif + *i);
	    }
	  
	//Go through cones, find all the ones containing a cone, remember links
	Vector<Set<int> > links;
	Set<int> containing_cones;
	for(int cone = 0; cone < cells.rows(); cone++) {
	    Set<int> rest;
	    if( (cells.row(cone) * tcone1).size() == tcone1.size() ) {
	      rest = cells.row(cone) - tcone1;
	    }
	    if( (cells.row(cone) * tcone2).size() == tcone2.size() ) {
	      rest = cells.row(cone) - tcone2;
	    }
	    if(rest.size() > 0) {
		containing_cones += cone;
		bool found_one = false;
		for(int l = 0; l < links.dim(); l++) {
		    if( (links[l] * rest).size() == rest.size() ) {
		      found_one = true; break;
		    }
		}
		if(!found_one) {
		  links |= rest;
		}
	    }	    
	}
	
	//Replace cones
	Vector<Set<int> > result;
	for(int c = 0; c < cells.rows(); c++) {
	    if(!containing_cones.contains(c)) result |= cells.row(c);
	}
	for(int l = 0; l < links.dim(); l++) {
	    result |= (links[l] + triang_simpl[0]);
	    result |= (links[l] + triang_simpl[1]);
	}
	
	return IncidenceMatrix<>(result);
       
       
     }//END doflip

//   // ------------------------- PERL WRAPPERS ---------------------------------------------------
//   
//   Function4perl(&computeMinMaxSpace,"cmms(RationalFunction)");

    Function4perl(&doflip,"doflip(IncidenceMatrix,Set<Int>,Set<Int>)");
  
}}