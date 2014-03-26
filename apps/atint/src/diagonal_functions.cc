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
 * Computes functions that cut out the diagonal on certain tropical varieties
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/linalg.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/minmax_functions.h"
#include "polymake/atint/diagonal_functions.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::ListReturn diagonal_unk(int n, int k) {
    perl::ListReturn p;
    
    //dbgtrace << "Creating flat functions... " << endl;
    
    //Create flat functions
    Array<perl::Object> flat_functions(k);
    for(int j = 1; j < k; j++) {
      Array<Set<int> > jsubsets = pm::all_subsets_of_k(sequence(0,n),j);
      for(int s = 0; s < jsubsets.size(); s++) {
	//Create max function for j-subset of [1,..,n]
	perl::Object f("MinMaxFunction");
	Matrix<Rational> xmatrix(2*j, n);
	Matrix<Rational> ymatrix(2*j, n);
	xmatrix.minor(sequence(0,j),jsubsets[s]) = unit_matrix<Rational>(j);
	ymatrix.minor(sequence(j,j), jsubsets[s]) = unit_matrix<Rational>(j);
	Matrix<Rational> fmatrix = xmatrix | ymatrix | zero_vector<Rational>();
	f.take("FUNCTION_MATRIX") << fmatrix;
	f.take("USES_MIN") << false;
	//Add it
	if(s == 0) flat_functions[j-1] = f;
	else flat_functions[j-1] = add_minmax_functions(flat_functions[j-1],f,true);
      }//END iterate j-subsets
    }//END iterate j-Flats
    //The last flat function is just the maximum over all coordinates
    Matrix<Rational> last_matrix = (unit_matrix<Rational>(2*n) | zero_vector<Rational>(2*n));
    perl::Object last_flat("MinMaxFunction");
      last_flat.take("FUNCTION_MATRIX") << last_matrix;
      last_flat.take("USES_MIN") << false;
    flat_functions[k-1] = last_flat;
          
    //dbgtrace << "Creating diagonal functions " << endl;
    
    //Now we create the diagonal functions
    for(int j = k; j >= 1; j--) {
      
      //dbgtrace << "Creating coefficients for diagonal function " << j << endl;
      
      //Create coefficient vector
      Vector<Integer> lambda = ones_vector<Integer>(k-j+1);
      for(int l = 1; l < lambda.dim()-1; l++) {
	//dbgtrace << "Computing coefficient " << l << endl;
	for(int b = j; b <= j+l-1; b++) {
// 	  dbgtrace << "Substracting " << (Integer::binom(j+l+1,j+l+1-b) * lambda[b-j]) << endl;
	  lambda[l] -= (Integer::binom(j+l,j+l-b) * lambda[b-j]);
	}
      }//END iterate over lambda
      //Create the last coefficient
      if( j < k) {
	for(int b = j; b <= k-1; b++) {
	  lambda[k-j] -= (Integer::binom(n,n-b) * lambda[b-j]);
	}
      }
      
      //dbgtrace << "Coefficients are " << lambda << endl;
      //dbgtrace << "Adding up positive and negative part " << endl;
      
      //Now add up the linear multiples of the MinMaxFunctions 
      //We add up positive and negative multiples separately to avoid intersecting domains as long as possible
      perl::Object positive_part("MinMaxFunction");
	bool pos_has_content = false;
      perl::Object negative_part("MinMaxFunction");
	bool neg_has_content = false;
      for(int l = 0; l < lambda.dim(); l++) {
	//dbgtrace << "Adding function " << l << endl;
	perl::Object scaled_function = scale_minmax_function(flat_functions[j+l-1],abs(lambda[l]));
	if(lambda[l] > 0) {
	  //dbgtrace << "... is positive " << endl;
	  positive_part = pos_has_content? add_minmax_functions(positive_part,scaled_function) : 
					   scaled_function;
	  pos_has_content = true;
	}
	if(lambda[l] < 0) {
	  //dbgtrace << "... is negative " << endl;
	  negative_part = neg_has_content? add_minmax_functions(negative_part, scaled_function) : 
					   scaled_function;
	  neg_has_content = true;
	}
      }
      //Now add positive and negative part
      perl::Object diag_function("MinMaxFunction");
      if(neg_has_content) {
	Matrix<Rational> neg_matrix = negative_part.give("FUNCTION_MATRIX");
	perl::Object neg_neg_part("MinMaxFunction");
	  neg_neg_part.take("FUNCTION_MATRIX") << (-neg_matrix);
	  neg_neg_part.take("USES_MIN") << true;
	diag_function = add_rational_functions(positive_part, neg_neg_part);
      }
      else {
	diag_function = positive_part;
      }
      p << diag_function;
      
    }//END create diagonal functions
    
    
    return p;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  //Documentation see header
  perl::Object lift_variety(perl::Object X) {
    Matrix<Rational> xrays = X.give("RAYS");
    Matrix<Rational> xlin = X.give("LINEALITY_SPACE");
    IncidenceMatrix<> cones = X.give("MAXIMAL_CONES");
    Vector<Integer> weights = X.give("TROPICAL_WEIGHTS");
    bool uses_homog = X.give("USES_HOMOGENEOUS_C");
    int ambdim = X.give("CMPLX_AMBIENT_DIM");
    
    //Lift rays
    xrays = xrays | zero_vector<Rational>();
    //Lift lineality
    if(xlin.rows() == 0) xlin = Matrix<Rational>(0,xrays.cols());
    else xlin = xlin | zero_vector<Rational>();
    Vector<Rational> add_generator = ones_vector<Rational>(ambdim+1);
      if(uses_homog) add_generator = Rational(0) | add_generator;
      
    xlin /= add_generator;
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << xrays; 
      result.take("MAXIMAL_CONES") << cones;
      result.take("LINEALITY_SPACE") << xlin;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
      
    return result;
  }
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  //Documentation see header
  perl::Object project_variety(perl::Object X) {
    //Extract values
    Matrix<Rational> xrays = X.give("RAYS");
    Matrix<Rational> xlin = X.give("LINEALITY_SPACE");
    IncidenceMatrix<> cones = X.give("MAXIMAL_CONES");
    Vector<Integer> weights = X.give("TROPICAL_WEIGHTS");
    bool uses_homog = X.give("USES_HOMOGENEOUS_C");
    
    //Truncate matrices in homogeneous case
    Vector<Rational> vertex_indicator;
    if(uses_homog) {
      if(xrays.rows() > 0) {
	vertex_indicator = xrays.col(0);
	xrays = xrays.minor(All,~scalar2set(0));
      }
      if(xlin.rows() > 0) {
	xlin = xlin.minor(All,~scalar2set(0));
      }
    }
    
    //Create projection matrix
    int cols = xlin.cols() > 0? xlin.cols() : xrays.cols();
    Matrix<Rational> projectionMatrix = unit_matrix<Rational>(cols-1);
    projectionMatrix /= (-ones_vector<Rational>(projectionMatrix.cols()));
    
    
    //Project
    if(xrays.rows() > 0) xrays = xrays * projectionMatrix;
    if(xlin.rows() > 0) xlin = xlin * projectionMatrix;
    
    Set<int> lbasis = basis_rows(xlin);
    xlin = xlin.minor(lbasis,All);
    
    //Put back into homogeneous coordinates if necessary
    if(uses_homog) {
      if(xrays.rows() > 0) xrays = vertex_indicator | xrays;
      if(xlin.rows() > 0) xlin = zero_vector<Rational>() | xlin ;
    }
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << xrays; 
      result.take("MAXIMAL_CONES") << cones;
      result.take("LINEALITY_SPACE") << xlin;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
    return result;
  }
//   //Finally, we mod out the lineality space (1..1) if necessary
//     if(modOutLineality) {
// 	int cols = bergman_lineality.cols();
// 	if(projCoordinate >= cols) projCoordinate = cols-1;
// 	//Create the projection matrix
// 	Matrix<Rational> unitMatrix = unit_matrix<Rational>(cols-1);
// 	Matrix<Rational> projectionMatrix(0,unitMatrix.cols());
// 	
// 	//Insert a -1's- vector at the right position
// 	if(projCoordinate > 0) {
// 	    projectionMatrix /= unitMatrix.minor(sequence(0,projCoordinate),All);
// 	}
// 	projectionMatrix /= - ones_vector<Rational>(unitMatrix.cols());
// 	if(projCoordinate < unitMatrix.rows()) {						
// 	      projectionMatrix /= unitMatrix.minor(sequence(projCoordinate,unitMatrix.rows() -	projCoordinate),All);
// 	}
// 	
// 	//dbgtrace << "Projection matrix is " << projectionMatrix << endl;
// 			   
// 	if(bergman_rays.rows() > 0) bergman_rays = bergman_rays * projectionMatrix;
// 	
// 	//Apply projection to the lineality space, but make sure the remaining rows are a basis
// 	bergman_lineality = bergman_lineality * projectionMatrix;
//     }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Intersection products"
  "# Computes a list of functions that cut out the diagonal on B(U<sup>n</sup><sub>k</sub>) x B(U<sup>n</sup><sub>k</sub>) "
		    "# @param Int n The ambient dimension of B(U<sup>n</sup><sub>k</sub>) "
		    "# @param Int k The dimension of B(U<sup>n</sup><sub>k</sub>) "
		    "# @return RationalFunction An array of functions, cutting out the diagonal",
		    &diagonal_unk, "diagonal_unk($,$)");
  
  Function4perl(&lift_variety,"lift_variety(WeightedComplex)");
  Function4perl(&project_variety,"project_variety(WeightedComplex)");
  
}}