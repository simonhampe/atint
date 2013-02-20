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
#include "polymake/linalg.h"
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object insert_rays(perl::Object fan, Matrix<Rational> add_rays) {
    
    //dbgtrace << "Triangulating fan " << endl;
    
    //Triangulate fan first
    fan = triangulateFan(fan);
    
    //Extract values
    Matrix<Rational> rays = fan.give("RAYS");
    Vector<Set<int> > cones = fan.give("MAXIMAL_CONES");
    bool weights_exist = fan.exists("TROPICAL_WEIGHTS");
    Vector<Integer> weights;
    if(weights_exist) {
      fan.give("TROPICAL_WEIGHTS") >> weights;
    }
    
    //First we check if any of the additional rays is already a ray of the fan
    
    //dbgtrace << "Checking if rays exist " << endl;
    
    Set<int> rays_to_check;
    for(int r = 0; r < add_rays.rows(); r++) {
      bool is_contained = false;
      for(int oray = 0; oray < rays.rows(); oray++) {
	if(rays.row(oray) == add_rays.row(r)) {
	    is_contained = true; break;
	}
      }
      if(!is_contained) rays_to_check += r;
    }//END check if rays are already in fan
    
    add_rays = add_rays.minor(rays_to_check,All);
    
    //dbgtrace << "Rays remaining: " << rays_to_check << endl;
    
    //dbgtrace << "Inserting remaining rays " << endl;
    
    //Now iterate over the remaining rays
    for(int nr = 0; nr < add_rays.rows(); nr++) {
      //dbgtrace << "Inserting ray nr. " << nr << endl;
      
      Vector<Set<int> > nr_cones;
      Vector<Integer> nr_weights;
      
      rays /= add_rays.row(nr);
      int nr_ray_index = rays.rows()-1;
      
      //Go through all cones and check if they contain this ray
      for(int mc = 0; mc < cones.dim(); mc++) {
	 //dbgtrace << "Checking cone nr. " << mc << endl;
	 bool contains_ray = false;
	 Matrix<Rational> relations = 
	  null_space( T(rays.minor(cones[mc],All) / add_rays.row(nr)));
	 //dbgtrace << "Relation is " << relations << endl;
	 //Since fan is simplicial, this matrix can have at most one relation
	 if(relations.rows() > 0) {
	    //Check if - assuming the last coefficient is < 0 - all other coeffs are >= 0
	    //Remember those that are > 0
	    if(relations(0,relations.cols()-1) > 0) relations *= -1;
	    contains_ray = true;
	    Set<int> greater_zero;
	    for(int c = 0; c < relations.cols()-1; c++) {
	      if(relations(0,c) < 0) {
		contains_ray = false;
	      }
	      if(relations(0,c) > 0) {
		greater_zero += c;
	      }
	    }//END check coefficients
	    
	    //If it is contained, subdivide the cone accordingly:
	    //For each non-zero-coefficient: take the codimension one face obtained
	    //by removing the corresponding ray and add the new ray
	    if(contains_ray) {
	      //dbgtrace << "Non-zero coeffs: " << greater_zero << endl;
	      Vector<int> rays_as_list(cones[mc]);
	      for(Entire<Set<int> >::iterator gzrays = entire(greater_zero); !gzrays.at_end(); gzrays++) {
		Set<int> nr_cone(rays_as_list.slice(~scalar2set(*gzrays)));
		nr_cone += nr_ray_index;
		nr_cones |= nr_cone;
		if(weights_exist) nr_weights |= weights[mc];
	      }//END iterate rays with coeff > 0
	    }//END subdivide cone
	    
	 }//END if there is a relation
	 
	 //If it does not contain the ray, just copy it
	 if(!contains_ray) {
	   //dbgtrace << "Does not contain, copying..." << endl;
	   nr_cones |= cones[mc];
	   if(weights_exist) nr_weights |= weights[mc];
	 }
      }//END iterate maximal cones
      
      //Copy new cones and weights
      cones = nr_cones;
      weights = nr_weights;
      
      //dbgtrace << "Cones are " << cones << endl;
      
    }//END iterate new rays
    
    //Create result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones;
      if(weights_exist) result.take("TROPICAL_WEIGHTS") << weights;
      
    return result;
    
    
  }//END insert_rays
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Basic polyhedral operations"
		    "# Takes a polyhedral fan and computes a triangulation"
		    "# @param WeightedComplex F A polyhedral fan (i.e. USES_HOMOGENEOUS_C = FALSE),"
		    "# possibly with weights"
		    "# @return WeightedComplex A simplicial refinement of F",
		    &triangulateFan,"triangulateFan(WeightedComplex)");
  
  UserFunction4perl("# @category Basic polyhedral operations"
		    "# Takes a polyhedral fan wihtout lineality and a list of rays and triangulates the fan"
		    "# such that it contains these rays"
		    "# @param WeightedComplex F A polyhedral fan (i.e. USES_HOMOGENEOUS_C = FALSE),"
		    "# possibly with weights but with no lineality space"
		    "# @param Matrix<Rational> R A list of normalized rays (as row vectors), which will "
		    "# be added to the fan"
		    "# @return WeightedComplex A triangulation of F that contains all the "
		    "# original rays of F plus the ones in R",
		    &insert_rays, "insert_rays(WeightedComplex,Matrix<Rational>)");
  
}}