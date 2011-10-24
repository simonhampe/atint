/*
 T his program is free software; you can redistribute it and/or                             *
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
 Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 This file contains the implementation of the generalized refinement function
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Integer.h"
#include "polymake/linalg.h"

namespace polymake { namespace atint { 
    
  using polymake::polytope::cdd_interface::solver;

  //using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;

  //Documentation see header
  RefinementResult refinement(perl::Object X, perl::Object Y, bool repFromX, bool repFromY,bool computeAssoc,bool refine) {
    solver<Rational> sv;
    
    //Extract values of the variety
    bool x_uses_homog = X.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> x_rays = X.give("RAYS");
    IncidenceMatrix<> x_cones = X.give("MAXIMAL_CONES");
    Matrix<Rational> x_lineality = X.give("LINEALITY_SPACE");
    int ambient_dim = x_rays.cols() < x_lineality.cols() ? x_lineality.cols() : x_rays.cols();
    int x_dimension = X.give("CMPLX_DIM");	
    Array<Integer> weights; bool weightsExist = false;
    if(X.exists("TROPICAL_WEIGHTS")) {
      weights = X.give("TROPICAL_WEIGHTS");
      weightsExist = true;	
    }
    
    dbgtrace << "Extracted X-values" << endl;
    
    //Extract values of the container
    bool y_uses_homog = Y.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> y_rays = Y.give("RAYS");
    IncidenceMatrix<> y_cmplx_cones = Y.give("CMPLX_MAXIMAL_CONES");
    IncidenceMatrix<> y_cones = Y.give("MAXIMAL_CONES");
    Matrix<Rational> y_lineality = Y.give("LINEALITY_SPACE");
      
    dbgtrace << "Extracted Y-values" << endl;
    
    //Prepare result variables
    perl::Object complex("WeightedComplex");
      Matrix<Rational> c_rays(0,ambient_dim);
      Matrix<Rational> c_lineality;
      int c_lineality_dim;
      Vector<Set<int> > c_cones;
      Vector<Integer> c_weights;
    Matrix<Rational> rayRepFromX;
    Matrix<Rational> rayRepFromY;
    Matrix<Rational> linRepFromX;
    Matrix<Rational> linRepFromY;
    Vector<int> associatedRep;
    
    dbgtrace << "Prepared result variables" << endl;
    
    //Step 1: Compute the lineality space
    if(refine) {
      //Compute the intersection of the two spaces
      Matrix<Rational> i_lineality = T(x_lineality  / (-y_lineality));
	dbgtrace << "Computing kernel of " << i_lineality << endl;
      Matrix<Rational> dependence =  null_space(i_lineality);
      c_lineality = dependence.minor(All,sequence(0,x_lineality.rows())) * x_lineality;
	dbgtrace << "Result: " << c_lineality << endl;
      c_lineality_dim = rank(c_lineality);
      //Compute X-rep if necessary
      if(repFromX) {
	linRepFromX = dependence.minor(All,sequence(0,x_lineality.rows()));
      }
      if(repFromY) {
	linRepFromY = dependence.minor(All,sequence(x_lineality.rows(),y_lineality.rows()));
      }
    }
    else {
      c_lineality = x_lineality;
      c_lineality_dim = X.give("LINEALITY_DIM");
      if(repFromX) {
	linRepFromX = unit_matrix<Rational>(x_lineality.rows());
      }
      if(repFromY) {
	linRepFromY = null_space(T(x_lineality / (-y_lineality))).minor(All,sequence(x_lineality.rows(),y_lineality.rows()));
      }
    }
    
    dbgtrace << "Computed lineality space" << endl;
    
    //Step 2: Compute cone refinement and ray representations.
    
    
    
//     
//     
//     //Prepare result variables
//     Matrix<Rational> newrays(0,ambient_dim);
//     Matrix<Rational> newlineality(0,ambient_dim);
//     Vector<Set<int> > newcones;
//     Vector<Integer> newweights;
//     Vector<Rational> newrayvalues;
//     Vector<Rational> newlinvalues;
//     //This variable maps indices of cones of the variety to ray index sets occuring in newcones.
//     //More precisely, for cone i, refinementSet[i] is a set of all new cones
//     //obtained by intersecting cone i.
//     Map<int,Set<Set<int> > > refinementSet;
//     
//     //Compute lineality space
//     Matrix<Rational> interlinmatrix = 
//       lineality_space | zero_matrix<Rational>(lineality_space.rows(),re_lineality_space.cols());    
//     interlinmatrix /= 
//       (zero_matrix<Rational>(re_lineality_space.rows(),lineality_space.cols()) | re_lineality_space);
//     newlineality = null_space(interlinmatrix);
//     newlineality = newlineality.minor(All,sequence(0,lineality_space.cols()));
//     int newlineality_dim = rank(newlineality);
//     
//     //If the variety or container is only a lineality space, we have to treat it as a cone for intersection
//     bool fanHasOnlyLineality = cones.rows() == 0;
//     bool containerOnlyHasLineality = re_cones.rows() == 0;
//     
//     //Before we start intersecting, we compute the H-reps of all container cones, so we don't 
//     //have to do it several times
//     Vector<Matrix<Rational> > inequalities;
//     Vector<Matrix<Rational> > equalities;
//     for(int c = 0; c < re_cones.cols(); c++) {
//       std::pair<Matrix<Rational>, Matrix<Rational> > facets = sv.enumerate_facets(zero_vector<Rational>() | re_rays.minor(re_cones.row(c),All),zero_vector<Rational>() | re_lineality_space, true,false);
//       inequalities |= facets.first;
//       equalities |= facets.second;
//     }
//     if(containerOnlyHasLineality) {
//       equalities |= null_space(re_lineality_space);
//       inequalities |= Matrix<Rational>(0,equalities[0].cols());
//     }
//     
//     //Iterate all maximal cones of variety
//     for(int vc = 0; vc < cones.rows() + (fanHasOnlyLineality? 1 : 0); vc++) {
//       //We compute the H-representation of the variety cone
//       std::pair<Matrix<Rational>,Matrix<Rational> > vfacet = sv.enumerate_facets(
// 	fanHasOnlyLineality? Matrix<Rational>(0,ambient_dim) : 
// 		      zero_vector<Rational>() | rays.minor(cones.row(vc),All),
// 		      zero_vector<Rational>() | lineality_space,true,false);
//       //Iterate all maximal cones of the container 
//       for(int co = 0; co < re_cones.rows() + (containerOnlyHasLineality? 1 : 0); co++) {
// 	//Compute the intersection rays
// 	Matrix<Rational> inter = sv.enumerate_vertices(
// 	    vfacet.first / inequalities[co], vfacet.second / equalities[co],true,true).first;
// 	inter = inter.minor(All,~scalar2set(0));
// 	
// 	//Check if the intersection has the correct dimension
// 	if(rank(inter) + newlineality_dim == dimension + (uses_homog? 1 : 0)) {
// 	  //Now we canonicalize the rays
// 	  for(int rw = 0; rw < inter.rows(); rw++) {
// 	    //Vertices start with a 1
// 	    if(uses_homog && inter(rw,0) != 0) {
// 	      inter.row(rw) /= inter(rw,0);
// 	    }
// 	    //The first non-zero entry in a directional ray is +-1
// 	    else {
// 	      for(int cl = 0; cl < inter.cols();cl++) {
// 		if(inter(rw,cl) != 0) {
// 		  inter.row(rw) /= abs(inter(rw,cl));
// 		  break;
// 		}
// 	      }
// 	    }
// 	  } //END canonicalize rays
// 	  
// 	  //Now we assign indices to the rays
// 	  int noOfRays = newrays.rows();
// 	  Set<int> rayIndices;
// 	  Set<int> newRayIndices;
// 	  for(int nr = 0; nr < inter.rows(); nr++) {
// 	    for(int r = 0; r < noOfRays; r++) {
// 	      if(newrays.row(r) == inter.row(nr)) {
// 		  rayIndices += r;
// 		  break;
// 	      }
// 	      if(r == noOfRays-1) {
// 		  newrays /= inter.row(nr);
// 		  rayIndices += (newrays.rows()-1);
// 		  newRayIndices += (newrays.rows()-1);
// 	      }
// 	    }
// 	  } //END assign ray indices
// 	  
// 	  //If there are no new rays, we have to check, whether this cone already
// 	  //exists (in this case we simply jump to the next case)
// 	  if(newRayIndices.size() == 0) {
// 	    if(refinementSet[vc].contains(rayIndices)) continue;
// 	  }
// 	  
// 	  //Now we add the new cone
// 	  newcones |= rayIndices;
// 	  refinementSet[vc] += rayIndices;
// 	  if(weightsExist) newweights |= weights[vc];
// 	  
// 	  //Now we compute the function values, if necessary
// 	  if(isFunction) {
// 	    //Iterate all new rays
// 	    for(Entire<Set<int> >::iterator nr = entire(newRayIndices); !nr.at_end(); nr++) {
// 	      //If the function is a general function, we just compute a representation of each
// 	      //new vector in the old rays and use this to compute its value
// 	      if(!isMinMaxFunction) {
// 		Vector<Rational> repv = functionRepresentationVector(re_cmplx_cones.row(co),newrays.row(*nr),							     ambient_dim,uses_homog,re_rays,re_lineality_space,
// 							    re_lineality_dim);
// 		newrayvalues |= (repv * values);
// 	      }
// 	      //Otherwise compute values for 
// 	      else {
// 		
// 	      }
// 	    }//END iterate all new rays
// 	  }//END compute function values
// 	  
// 	}//END if full-dimensional
// 	
// 	
//       }//END iterate container cones
//     }//END iterate variety cones
    
    //Insert values
    
    //Copy return values
    if(refine) {
      complex.take("RAYS") << c_rays;
      complex.take("MAXIMAL_CONES") << c_cones;
      complex.take("LINEALITY_SPACE") << c_lineality;
      complex.take("USES_HOMOGENEOUS_C") << x_uses_homog;
      if(weightsExist) complex.take("TROPICAL_WEIGHTS") << c_weights;
    }
    else {
      complex = X;
    }
    RefinementResult result;
      result.complex = complex;
      result.rayRepFromX = rayRepFromX;
      result.rayRepFromY = rayRepFromY;
      result.linRepFromX = linRepFromX;
      result.linRepFromY = linRepFromY;
      result.associatedRep = associatedRep;
    return result;
    
  }//END function refine
  
  //TODO: DEBUG. REMOVE
  perl::Object reftest(perl::Object X, perl::Object Y, bool repFromX, bool repFromY,bool computeAssoc,bool refine) {
    RefinementResult r;
    r = refinement(X, Y, repFromX, repFromY,computeAssoc,refine);
    pm::cout << "Xrayrep: " << r.rayRepFromX << endl;
    pm::cout << "Xlinrep: " << r.linRepFromX << endl;
    pm::cout << "Yrayrep: " << r.rayRepFromY << endl;
    pm::cout << "Ylinrep: " << r.linRepFromY << endl;
    pm::cout << "assoc rep: " << r.associatedRep << endl;
    return r.complex;
  }

// ------------------------- PERL WRAPPERS ---------------------------------------------------

Function4perl(&reftest,"reftest(WeightedComplex,WeightedComplex,$,$,$,$)");

}}

