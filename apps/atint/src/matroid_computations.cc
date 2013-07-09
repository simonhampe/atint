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
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 This file provides convenience methods for creating locally restricted
 tropical varieties from given varieties
 */


#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  perl::Object matroidFromFan(perl::Object fan) {
    
    //Dehomogenize, if necessary.
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    if(uses_homog) {
	fan = fan.CallPolymakeMethod("dehomogenize");
    }
    
    //Extract values
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    IncidenceMatrix<> cones = fan.give("MAXIMAL_CONES");
    int rank = fan.give("CMPLX_DIM");
    int no_of_el = fan.give("CMPLX_AMBIENT_DIM");
    
    //Compute all potential bases
    Set<Set<int> > bases(pm::all_subsets_of_k(sequence(0,no_of_el),rank));
    
    //Now go through all potential flat vectors of cardinality >= rank. The corr. set is a flat, iff
    // the fan contains that vector. If it is a flat, remove all bases contained in it.
    for(int card = rank; card < no_of_el; card++) {
	Array<Set<int> > potential_flats = pm::all_subsets_of_k(sequence(0,no_of_el),card);
	for(int flat = 0; flat < potential_flats.size(); flat++) {
	    Vector<Rational> inc_vec(no_of_el);
	      inc_vec.slice(potential_flats[flat]) = -ones_vector<Rational>(card);
	    bool is_in_fan = false;
	    for(int mc = 0; mc < cones.rows(); mc++) {
	      if(is_ray_in_cone(rays.minor(cones.row(mc),All),linspace,inc_vec)) {
		  is_in_fan = true; break;
	      }
	    }//END iterate maximal cones
	    if(is_in_fan) {
	      //Compute all subsets of cardinality the rank and remove them from the bases
	      bases -= Set<Set<int> >(pm::all_subsets_of_k(potential_flats[flat],rank));	      
	    }
	}//END iterate flats of given cardinality
	
    }//END iterate flat cardinalities
    
    Vector<Set<int> > final_bases(bases);
    
    perl::Object result("matroid::Matroid");
      result.take("N_ELEMENTS") << no_of_el;
      result.take("BASES") << final_bases;
      
    return result;
    
  }//END matroidFromFan
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Matroids"
		    "# Takes a tropical fan (that should have degree 1) and computes the matroid that has this fan"
		    "# as bergman fan"
		    "# @param WeightedComplex fan A tropical fan that should be the Bergman fan of a matroid (NOT"
		    "# modulo lineality space). Be aware that this implementation uses (-1/0)-incidence vectors"
		    "# for flats, not (1/0)."
		    "# @return matroid::Matroid The corresponding matroid (given in terms of bases)",
		    &matroidFromFan,"matroid_from_fan(WeightedComplex)");
  
}}