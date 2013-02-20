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
 * This file contains functions to check linear equivalence of cycles
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/linalg.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  bool check_cycle_equality(perl::Object X, perl::Object Y, bool check_weights = true) {
    
    //dbgtrace << "Extracting values " << endl;
    //Extract values
    Matrix<Rational> xrays = X.give("RAYS");
    IncidenceMatrix<> xcones = X.give("MAXIMAL_CONES");
    Matrix<Rational> xlin = X.give("LINEALITY_SPACE");
    int xambi = X.give("CMPLX_AMBIENT_DIM");
    Vector<Integer> xweights;
    if(X.exists("TROPICAL_WEIGHTS")) {
      X.give("TROPICAL_WEIGHTS") >> xweights;
    }
    else check_weights = false;
    
    Matrix<Rational> yrays = Y.give("RAYS");
    IncidenceMatrix<> ycones = Y.give("MAXIMAL_CONES");
    Matrix<Rational> ylin = Y.give("LINEALITY_SPACE");
    int yambi = Y.give("CMPLX_AMBIENT_DIM");
    Vector<Integer> yweights;
    if(Y.exists("TROPICAL_WEIGHTS")) {
      Y.give("TROPICAL_WEIGHTS") >> yweights;
    }
    else check_weights = false;
    
    //dbgtrace << "Checking equality of dimensions " << endl;
    
    //Check dimensional equality
    if(xambi != yambi) return false;
    
    //dbgtrace << "Checkking equality of lineality spaces " << endl;
    
    //Check equality of lineality spaces
    if(rank(xlin) == rank(ylin)) {
      if(rank(xlin / ylin) > rank(xlin)) return false;
    }
    else return false;
    
    //dbgtrace << "Finding ray permutation " << endl;
    
    //Find ray permutation
    if(xrays.rows() != yrays.rows()) return false;
    Map<int,int> permutation;
    for(int x = 0; x < xrays.rows(); x++) {
      for(int y = 0; y < yrays.rows(); y++) {
	//Check if yray = xray modulo lineality
	Matrix<Rational> xray_matrix = xrays.row(x) / xlin;
	if(rank(xray_matrix) == rank(xray_matrix / yrays.row(y))) {
	    permutation[x] = y;
	    break;
	}
	//If we arrive here, there is no ray matching x
	if(y == yrays.rows()-1) return false;
      }
    }//END compute ray permutation
    
    //dbgtrace << "Ray permutation is " << permutation << endl;
    
    //dbgtrace << "Matching cones " << endl;
    
    //Now check if all cones are equal
    Set<int> matched_cones;
    for(int xc = 0; xc < xcones.rows(); xc++) {
      //dbgtrace << "Matching cone " << xcones.row(xc) << endl;
      //Compute permuted cone
      Set<int> perm_cone = attach_operation(xcones.row(xc), pm::operations::associative_access<Map<int,int>, int>(&permutation));
      //dbgtrace << "Permuted cone is " << perm_cone << endl;
      //Find this cone in Y
      for(int yc = 0; yc < ycones.rows(); yc++) {
	if(!matched_cones.contains(yc)) {
	  if(ycones.row(yc).size() == perm_cone.size()) {
	    //dbgtrace << "Comparing with " << ycones.row(yc) << endl;
	      if((ycones.row(yc) * perm_cone).size() == perm_cone.size()) {
		matched_cones += yc;
		//Check equality of weights, if necessary
		if(check_weights) {
		  //dbgtrace << "Checking weight equality" << endl;
		  if(xweights[xc] != yweights[yc]) return false;
		}
		break;
	      }
	  }//END check cone equality
	}
	//If we arrive here, there is no match:
	if(yc == ycones.rows()-1) return false;
	
      }//END iterate Y cones
    }//END iterate X cones
    
    //dbgtrace << "Checking if all cones were matched " << endl;
    
    //Check if we actually matched ALL Y cones
    if(matched_cones.size() < ycones.rows()) return false;
    
    return true;
    
  }//END function check_cycle_equality
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Basic polyhedral properties"
		    "# This takes two pure-dimensional polyhedral complexes and checks if they are equal"
		    "# i.e. if they have the same lineality space, the same rays (modulo lineality space)"
		    "# and the same cones. Optionally, it can also check if the weights are equal"
		    "# @param WeightedComplex X A weighted complex"
		    "# @param WeightedComplex Y A weighted complex"
		    "# @param Bool check_weights Whether the algorithm should check for equality of weights. "
		    "# This parameter is optional and true by default"
		    "# @return Bool Whether the cycles are equal",
		    &check_cycle_equality, "check_cycle_equality(WeightedComplex,WeightedComplex;$=1)");
  
}}