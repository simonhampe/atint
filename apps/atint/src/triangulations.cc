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
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/WeightedComplexRules.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object triangulateFan(perl::Object fan) {
    //Extract values
    Matrix<Rational> rays = fan.give("RAYS");
    Vector<Set<int> > cones = fan.give("MAXIMAL_CONES");
    Vector<Set<int> > codim = fan.give("CODIM_1_FACES");
    IncidenceMatrix<> codimInMax = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> codimInMaxTranspose = T(codimInMax);
    //Convert transpose to vector of sets
    Vector<Set<int> > maxHasCodim;
      for(int r = 0; r < codimInMaxTranspose.rows(); r++) {
	maxHasCodim |= codimInMaxTranspose.row(r);
      }
    
    Matrix<Rational> lin_space = fan.give("LINEALITY_SPACE");
    int lineality_dim = fan.give("LINEALITY_DIM");
    int fan_dim = fan.give("FAN_DIM");
    
    bool weights_exist = fan.exists("TROPICAL_WEIGHTS");
    Vector<Integer> weights;
    if(weights_exist) {
      fan.give("TROPICAL_WEIGHTS") >> weights;
    }
    
    int no_of_rays = fan_dim - lineality_dim;
    
    
    //This will contain the indices of cones we keep at the end
    Set<int> listOfCones = sequence(0,cones.dim());
    
    std::list<int> queue;
      for(int i = 0; i < cones.dim(); i++) { queue.push_back(i);}
    
    //We do an inductive barycentric triangulation of each cone: If it has too many rays,
    //add an interior ray, replace the cone by all cones spanned by a codimension one face and this
    //additional ray. Then repeat the process for the newly added cones.
    while(queue.size() > 0) {
      int sigma = queue.front();
	queue.pop_front();
	
      //If the cone has too many rays, we have to triangulate
      if(cones[sigma].size() > no_of_rays) {
	//Remove cone from used Cones
	listOfCones -= sigma;
	
	//Create new ray as sum of old ones and normalize
	Vector<Rational> nray = accumulate(rows(rays.minor(cones[sigma],All)), operations::add());
	for(int c = 0; c < nray.dim(); c++) {
	  if(nray[c] != 0) {
	    nray /= abs(nray[c]);
	    break;
	  }
	}
	
	//Add ray to list of rays
	rays /= nray;
	int nray_index = rays.rows()-1;
	
	//Add a new cone for each codimension one face
	Set<int> sigmafaces = maxHasCodim[sigma];
	for(Entire<Set<int> >::iterator tau = entire(sigmafaces); !tau.at_end(); tau++) {
	    //Add new maximal cone
	    Set<int> new_maximal = codim[*tau] + nray_index;
	    Vector<Set<int> > new_maximal_list; new_maximal_list |= new_maximal;
	    cones |= new_maximal;
	    int ncone_index = cones.dim()-1;
	    maxHasCodim |= Set<int>();
	    queue.push_back(ncone_index);
	    listOfCones += ncone_index;
	    if(weights_exist) {
	      weights |= weights[sigma];
	    }
	    
	    //Compute the new codimension one faces
	    CodimensionOneResult coresult = 
	      calculateCodimOneData(rays, new_maximal_list, false, lin_space, IncidenceMatrix<>());
	    for(int r = 0; r < coresult.codimOneCones.rows(); r++) {
	      codim |= coresult.codimOneCones.row(r);
	      maxHasCodim[ncone_index] += (codim.dim()-1);
	    }
	}//END iterate codim one faces
      }//END if not simplicial
    }//END iterate queue
    
    //Create result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones.slice(listOfCones);
      result.take("LINEALITY_SPACE") << lin_space;
      if(weights_exist) {
	result.take("TROPICAL_WEIGHTS") << weights.slice(listOfCones);
      }
      
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