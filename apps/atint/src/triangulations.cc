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
 * This file contains functions to compute triangulations (of fans, complexes,...)
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/polytope/beneath_beyond.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/WeightedComplexRules.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object triangulateFan(perl::Object fan) {
    Matrix<Rational> rays = fan.give("RAYS");
    IncidenceMatrix<> cones = fan.give("MAXIMAL_CONES");
    Matrix<Rational> lin_space = fan.give("LINEALITY_SPACE");
    
    int fan_dim = fan.give("FAN_DIM");
    
    bool weights_exist = fan.exists("TROPICAL_WEIGHTS");
    Vector<Integer> weights;
    if(weights_exist) {
      fan.give("TROPICAL_WEIGHTS") >> weights;
    }
    
    int no_of_rays = fan_dim - lin_space.rows();
    
    //Will contain the triangulating cones and their weights
    Vector<Set<int> > triangleCones;
    Vector<Integer> triangleWeights;
    
    //dbgtrace << "Extracted values" << endl;
    
    //Go through all cones
    for(int sigma = 0; sigma < cones.rows(); sigma++) {
	//dbgtrace << "Subdividing cone " << sigma << " = " << cones.row(sigma) << endl;
	if(cones.row(sigma).size() > no_of_rays) {
	    //dbgtrace << "Setting up beneath and beyond algorithm" << endl;
	    Matrix<Rational> sigmarays = rays.minor(cones.row(sigma),All);
	    polytope::beneath_beyond_algo<Rational> algo( sigmarays,false);
	    //dbgtrace << "Computing data" << endl;
	    algo.compute(entire(sequence(0,sigmarays.rows())));
	    Array<Set<int> > triang_seq = algo.getTriangulation();
	    //dbgtrace << "Triangulation reads " << triang_seq << endl;
	    Vector<int> rays_as_list(cones.row(sigma));
	    for(int tr = 0; tr < triang_seq.size(); tr++) {
	      triangleCones |= Set<int>(rays_as_list.slice(triang_seq[tr]));
	      if(weights_exist) triangleWeights |= weights[sigma];
	    }
	}
	else {
	    //dbgtrace << "Is already simplicial\n" << endl;
	    triangleCones |= cones.row(sigma);
	    if(weights_exist) {
	      triangleWeights |= weights[sigma];
	    }
	}	
	//dbgtrace << "Weights are " << triangleWeights << endl;
    }//END iterate all cones
    
    //dbgtrace << "Triangulation complete " << endl;
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << triangleCones;
      result.take("LINEALITY_SPACE") << lin_space;
      if(weights_exist) result.take("TROPICAL_WEIGHTS") << triangleWeights;
    
    return result;
  }//END triangulateFan
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Triangulations"
		    "# Takes a polyhedral fan and computes a triangulation"
		    "# @param WeightedComplex F A polyhedral fan (i.e. USES_HOMOGENEOUS_C = FALSE),"
		    "# possibly with weights"
		    "# @return WeightedComplex A simplicial refinement of F",
		    &triangulateFan,"triangulateFan(WeightedComplex)");
  
}}