/*
 This *program is free software; you can redistribute it and/or
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
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/WeightedComplexRules.h"

namespace polymake { namespace atint { 
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  perl::Object local_restrict(perl::Object complex, IncidenceMatrix<> cones) {
    //Extract values
    Vector<Set<int> > maximalCones = complex.give("MAXIMAL_CONES");
    Matrix<Rational> rays = complex.give("RAYS");
    Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
    
    //Find out which cones are no longe compatible
    Set<int> remainingCones;
    for(int c = 0; c < maximalCones.dim(); c++) {
      if(is_coneset_compatible(maximalCones[c], cones)) {
	remainingCones += c;
      }
    }
    
    //Adapt cone description and ray indices
    maximalCones = maximalCones.slice(remainingCones);
    weights = weights.slice(remainingCones);
    Set<int> usedRays = accumulate(maximalCones,operations::add());
    rays = rays.minor(usedRays,All);
    IncidenceMatrix<> newMaximalCones(maximalCones);
      newMaximalCones = newMaximalCones.minor(All,usedRays);
    cones = cones.minor(All,usedRays);
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << newMaximalCones;
      result.take("LINEALITY_SPACE") << linspace;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("LOCAL_RESTRICTION") << cones;
    
    return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry /Local geometry"
		    "# This takes a tropical variety and an IncidenceMatrix describing a set"
		    "# of cones (not necessarily maximal ones) of this variety. It will then"
		    "# create a variety that contains all compatible maximal cones and is"
		    "# locally restricted to the given cone set."
		    "# @param WeightedComplex complex An arbitrary weighted complex"
		    "# @param IncidenceMatrix cones A set of cones, indices refer to RAYS"
		    "# @return WeightedComplex The same complex, locally restricted to the given"
		    "# cones",
		    &local_restrict, "local_restrict(WeightedComplex,IncidenceMatrix)");
  
}}


