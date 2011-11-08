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
#include "polymake/atint/normalvector.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Integer.h"
#include "polymake/linalg.h"

namespace polymake { namespace atint { 
    
  using polymake::polytope::cdd_interface::solver;

  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;

  //Documentation see header
  RefinementResult refinement(perl::Object X, perl::Object Y, bool repFromX, bool repFromY,bool computeAssoc,bool refine) {
    solver<Rational> sv;
    
    //Extract values of the variety
    bool x_uses_homog = X.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> x_rays = X.give("RAYS");
    IncidenceMatrix<> x_cones = X.give("MAXIMAL_CONES");
    Matrix<Rational> x_cmplx_rays = repFromX? X.give("CMPLX_RAYS") : Matrix<Rational>();
    IncidenceMatrix<> x_cmplx_cones = repFromX? X.give("CMPLX_MAXIMAL_CONES") : IncidenceMatrix<>();
    Matrix<Rational> x_lineality = X.give("LINEALITY_SPACE");
    int x_lineality_dim = X.give("LINEALITY_DIM");
    int ambient_dim = x_rays.cols() < x_lineality.cols() ? x_lineality.cols() : x_rays.cols();
    int x_dimension = X.give("CMPLX_DIM");	
    Array<Integer> weights; bool weightsExist = false;
    if(X.exists("TROPICAL_WEIGHTS")) {
      weights = X.give("TROPICAL_WEIGHTS");
      weightsExist = true;	
    }
    
    dbgtrace << "Extracted X-values" << endl;
    
    //Extract values of the container
    Matrix<Rational> y_rays = Y.give("RAYS");
    Matrix<Rational> y_cmplx_rays = repFromY? Y.give("CMPLX_RAYS") : Matrix<Rational>();
    IncidenceMatrix<> y_cmplx_cones = repFromY? Y.give("CMPLX_MAXIMAL_CONES") : IncidenceMatrix<>();
    IncidenceMatrix<> y_cones = Y.give("MAXIMAL_CONES");
    Matrix<Rational> y_lineality = Y.give("LINEALITY_SPACE");
    int y_lineality_dim = Y.give("LINEALITY_DIM");
      
    dbgtrace << "Extracted Y-values" << endl;
    
    //Prepare result variables
    perl::Object complex("WeightedComplex");
      Matrix<Rational> c_rays(0,ambient_dim);
      Matrix<Rational> c_lineality(0,ambient_dim);
      int c_lineality_dim = 0;
      Vector<Set<int> > c_cones;
      Vector<Integer> c_weights;
    Matrix<Rational> rayRepFromX;
    Matrix<Rational> rayRepFromY;
    Matrix<Rational> linRepFromX;
    Matrix<Rational> linRepFromY;
    Vector<int> associatedRep;
    
    //If we don't refine, we already know most of the results concerning X
    if(!refine) {
      c_rays = x_rays;
      for(int xr = 0; xr < x_cones.rows(); xr++) c_cones |= x_cones.row(xr);
      c_weights = weights;
      //if(computeAssoc) associatedRep = Vector<int>(x_cones.rows());
      if(repFromX) rayRepFromX = unit_matrix<Rational>(x_cmplx_cones.rows()) |
		    zero_matrix<Rational>(x_cmplx_cones.rows(),x_lineality.rows());
    }
    
    dbgtrace << "Prepared result variables" << endl;
    
    //Step 1: Compute the lineality space ----------------------------------
    if(x_lineality.rows() != 0 && y_lineality.rows() != 0) {
      if(refine) {
	//Compute the intersection of the two spaces
	//We compute the kernel of (x_lineality | -y_lineality)
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
	  //We compute the kernel of (x_lineality | -y_lineality) and invert the first
	  //half to obtain a description of x_lineality in terms of y_lineality
	  linRepFromY = null_space(T(x_lineality / (-y_lineality)));
	  linRepFromY = (inv(linRepFromY.minor(All,sequence(0,x_lineality.rows()))) * linRepFromY).minor(All,sequence(x_lineality.rows(),y_lineality.rows()));
	}
      }
    }
    
    dbgtrace << "Computed lineality space" << endl;
    
    //Step 2: Compute cone refinement and ray representations. -----------------
    
    //If any of the two is just a lineality space, we still have to consider this as a single
    //maximal cone for refinement / representation
    bool x_onlylineality = x_cones.rows() == 0; 
    bool y_onlylineality = y_cones.rows() == 0; 
    
    //Compute a facet representation for the Y-cones
    Vector<std::pair<Matrix<Rational>,Matrix<Rational> > > y_equations;
    for(int yc = 0; yc < (y_onlylineality? 1 : y_cones.rows()); yc++) {
      y_equations |= 
	sv.enumerate_facets( zero_vector<Rational>() | (
	   y_onlylineality? Matrix<Rational>(0,y_lineality.cols()) : y_rays.minor(y_cones.row(yc),All)),
	   zero_vector<Rational>() | y_lineality, true,false);
    }
    
    //This saves for each x-cone (index in x_cones) the cones that have been created as refinements.
    //Since doubles can only occur when the same x-cone is contained in the intersection of several y-cones,
    //we only have to check the cones in xrefinements[xc]
    Vector< Set<Set<int> > > xrefinements(x_onlylineality? 1 :  x_cones.rows());
    
    //These variables save for each cone of the intersection the indices of the cones in X and Y containing it
    Vector<int> xcontainers;
    Vector<int> ycontainers;
    
    if(refine || repFromX || repFromY) {
      //Iterate all cones of X
      for(int xc = 0; xc < (x_onlylineality? 1 : x_cones.rows()); xc++)  {
	//Initalize refinement cone set
	xrefinements[xc] = Set<Set<int> >();
	//Compute a facet representation for the X-cone
	std::pair<Matrix<Rational>, Matrix<Rational> > x_equations = 
	  sv.enumerate_facets( zero_vector<Rational>() | (
	    x_onlylineality? Matrix<Rational>(0,x_lineality.cols()) : x_rays.minor(x_cones.row(xc),All)),
	    zero_vector<Rational>() | x_lineality, true,false);
	  
	//Iterate all cones of Y
	for(int yc = 0; yc < (y_onlylineality? 1 : y_cones.rows()); yc++) {
	  //Compute a V-representation of the intersection
	  Matrix<Rational> interrays = sv.enumerate_vertices(x_equations.first / y_equations[yc].first,
		      x_equations.second / y_equations[yc].second,true,true).first;
	  interrays = interrays.minor(All,~scalar2set(0));
	  
	  //Check if it is full-dimensional (and has at least one ray - lin.spaces are not interesting)
	  if(interrays.rows() > 0 && rank(interrays) + c_lineality_dim - (x_uses_homog? 1 : 0) == x_dimension) {
	    //If we refine, add the cone. Otherwise just remember the indices
	    dbgtrace << "Inter rays: " << interrays << endl;
	    Set<int> interIndices;
	    if(!refine) {
	      //Copy indices
	      interIndices = x_cones.row(xc);
	      //Compute assoc-rep. and identify new rays
  // 	    int vertex = -1; Set<int> newDirectionalRays;
  // 	    for(Entire<Set<int> >::iterator ir = entire(interIndices); !ir.at_end(); ir++) {
  // 	      if(computeAssoc) {
  // 		if(c_rays(*ir,0) == 1) {
  // 		    if(vertex == -1) vertex = *ir;
  // 		    associatedRep[*ir] = *ir;
  // 		}
  // 		else newDirectionalRays += (*ir);
  // 	      }
  // 	    }
  // 	    if(computeAssoc) {
  // 	      for(Entire<Set<int> >::iterator dr = entire(newDirectionalRays); !dr.at_end(); dr++) {
  // 		associatedRep[*dr] = vertex;
  // 	      }
  // 	    }	    
	    }
	    else {
	      //Now we canonicalize the rays and assign ray indices
	      //We also use this opportunity to sort the rays 
	      Set<int> newRays;
  // 	    , newDirectionalRays;
  // 	    int associatedVertex = -1;
	      for(int rw = 0; rw < interrays.rows(); rw++) {
		//bool isVertex = false;
		//Vertices start with a 1
		if(x_uses_homog && interrays(rw,0) != 0) {
		  interrays.row(rw) /= interrays(rw,0);
		  //isVertex = true;
		}
		//The first non-zero entry in a directional ray is +-1
		else {
		  for(int cl = 0; cl < interrays.cols();cl++) {
		    if(interrays(rw,cl) != 0) {
		      interrays.row(rw) /= abs(interrays(rw,cl));
		      break;
		    }
		  }
		}
		
		dbgtrace << "Considering row " << interrays.row(rw) << endl;
		
		//Go through the existing rays and compare
		int nrays = c_rays.rows();
		int newrayindex = -1;
		for(int oray = 0; oray < nrays; oray++) {
		  if(interrays.row(rw) == c_rays.row(oray)) {
		      newrayindex = oray;
		      break;
		  }
		}
		if(newrayindex == -1) {
		  c_rays /= interrays.row(rw);
		  newrayindex = c_rays.rows()-1;
		  newRays += newrayindex;
  // 		if(!isVertex) { newDirectionalRays += newrayindex;}
  // 		//Also add rows/entries in the assoc representation variable
  // 		if(computeAssoc) associatedRep |= 0;
		}
		interIndices += newrayindex;
		//If this is the first vertex, save the index separately
  // 	      if(isVertex && computeAssoc) {
  // 		associatedRep[newrayindex] = newrayindex;
  // 		if(associatedVertex == -1) associatedVertex = newrayindex;	      		
  // 	      }
	      } //END canonicalize rays
	      
	      dbgtrace << "Ray indices " << interIndices << endl;
	      dbgtrace << "new rays: " << newRays << endl;
  // 	    dbgtrace << "directional: " << newDirectionalRays << endl;
  // 	    dbgtrace << "Associated vertex: " << associatedVertex << endl;
	      
	      //Check if the cone exists - if there are new rays, then the cone must be new as well
	      bool addCone = newRays.size() > 0;
	      if(!addCone) addCone = !xrefinements[xc].contains(interIndices);
	      //If the cone is new, add it
	      if(addCone) {
		dbgtrace << "Adding new cone" << endl;
		c_cones |= interIndices;
		if(weightsExist) c_weights |= weights[xc];
		xrefinements[xc] += interIndices;
		xcontainers |= xc;
		ycontainers |= yc;
	      }
	      
	    /*  //Compute assocRep
	      if(computeAssoc) {
		for(Entire<Set<int> >::iterator dr = entire(newDirectionalRays); !dr.at_end(); dr++) {
		  associatedRep[*dr] = associatedVertex;
		}
	      }*/	    
	    } //END canonicalize intersection cone and add it
		    
	    //If we do not refine, we only need to find one y-cone containing the x-cone
	    if(!refine) break;	  
	  } //END if full-dimensional
	}//END iterate y-cones
      }//END iterate x-cones
    } //END if intersection is necessary?
    
    //Copy return values into the fan
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
    
    //To compute representations of CMPLX_RAYS, we naturally have to compute the CMPLX_RAYS first
    if((repFromX && refine) || repFromY || computeAssoc) {
      dbgtrace << "Computing representations" << endl;
      Matrix<Rational> c_cmplx_rays = complex.give("CMPLX_RAYS");
      IncidenceMatrix<> c_cmplx_cones = complex.give("CMPLX_MAXIMAL_CONES");
      //Initialize rep matrices to proper size
      rayRepFromX = Matrix<Rational>(c_cmplx_rays.rows(),x_cmplx_rays.rows() + x_lineality.rows());
      rayRepFromY = Matrix<Rational>(c_cmplx_rays.rows(),y_cmplx_rays.rows() + y_lineality.rows());
      //Compute representations for X (mode 0) and/or Y (mode 1)
      for(int mode = 0; mode <= 1; mode++) {
	dbgtrace << "Computing in mode " << mode << endl;
	if((mode == 0 && repFromX) || (mode == 1 && repFromY)) {
	    //Recalls for which ray we already computed a representation
	    Vector<bool> repComputed(c_cmplx_rays.rows());
	    Matrix<Rational> raysForComputation = (mode == 0? x_cmplx_rays : y_cmplx_rays);
	    Matrix<Rational> linForComputation = (mode == 0? x_lineality : y_lineality);
	    int dimForComputation = (mode == 0? x_lineality_dim : y_lineality_dim);
	    //Go through all complex cones
	    for(int cone = 0; cone < c_cmplx_cones.rows(); cone++) {
	      dbgtrace << "Computing rep in cone " << cone << endl;
	      //Go through all rays for which we have not yet computed a representation
	      Set<int> raysOfCone = c_cmplx_cones.row(cone);
	      for(Entire<Set<int> >::iterator r = entire(raysOfCone); !r.at_end(); r++) {
		if(!repComputed[*r]) {
		  repComputed[*r] = true;
		  //The ray used for computing the representation are the rays of the containing
		  //cone (or none, if the corr. fan is only a lineality space)
		  Set<int> rfc;
		  if(mode == 0 && !x_onlylineality) 
		      rfc = x_cmplx_cones.row(xcontainers[cone]);
		  if(mode == 1 && !y_onlylineality)
		      rfc = y_cmplx_cones.row(ycontainers[cone]);
		  //Compute representation vector
		  (mode == 0? rayRepFromX : rayRepFromY).row(*r) =
		      functionRepresentationVector(
			rfc,
			c_cmplx_rays.row(*r),
			ambient_dim,
			false,
			raysForComputation,
			linForComputation, dimForComputation); 				   
		}
	      }
	    }//END iterate all cones
	}//END if mode
      }//END go through X- and Y-mode
      
      if(computeAssoc) {
	//For each cmplx_ray, in which cone does it lie?
	IncidenceMatrix<> c_cmplx_cones_t = T(c_cmplx_cones);
	for(int cr = 0; cr < c_cmplx_rays.rows(); cr++) {
	    if(c_cmplx_rays(cr,0) == 1) {
	      associatedRep |= cr;
	    }
	    else {
	      //Take a cone containing ray cr
	      int mc = *(c_cmplx_cones_t.row(cr).begin());
	      //Find an affine ray in this cone
	      Set<int> mcrays = c_cmplx_cones.row(mc);
	      for(Entire<Set<int> >::iterator mcr = entire(mcrays); !mcr.at_end(); mcr++) {
		if(c_cmplx_rays(*mcr,0) == 1) {
		    associatedRep |= *mcr;
		    break;
		}
	      }
	    }
	}
      }//END if computeAssoc      
    }//END compute nontrivial representations
        
    
    //Insert values
    
    RefinementResult result;
      result.complex = complex;
      result.rayRepFromX = rayRepFromX;
      result.rayRepFromY = rayRepFromY;
      result.linRepFromX = linRepFromX;
      result.linRepFromY = linRepFromY;
      result.associatedRep = associatedRep;
    return result;
    
  }//END function refine
  
  
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

