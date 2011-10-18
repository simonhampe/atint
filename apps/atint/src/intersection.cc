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
Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/atint/basicoperations.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/refinement.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/linalg.h"

namespace polymake { namespace atint {

    using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    /**
     @brief Compute the intersection product of two tropical cycles
     @param perl::Object X A WeightedComplex, assumed to be balanced
     @param perl::Object Y A balances WeightedComplex with the same ambient dimension as X
     */
    perl::Object cycle_intersection(perl::Object X, perl::Object Y) {
      dbgtrace << "Extracting values" << endl;
      //Extract values
      int Xcodim = X.give("CMPLX_CODIMENSION");
      int Ycodim = Y.give("CMPLX_CODIMENSION");
      int Xambi  = X.give("CMPLX_AMBIENT_DIM");
      
      dbgtrace << "Checking codimension" << endl;
      
      //If the codimensions of the varieties add up to something larger then CMPLX_AMBIENT_DIM, return the 0-cycle 
      if(Xcodim + Ycodim > Xambi) {
	Matrix<Rational> norays(0,Xambi);
	Vector<Set<int> > nofacets;
	perl::Object zerocycle("WeightedComplex");
	  zerocycle.take("INPUT_RAYS") << norays;
	  zerocycle.take("INPUT_CONES") << nofacets;
	return zerocycle;
      }
      
      dbgtrace << "Computing product" << endl;
      
      //Compute the cross product
      std::vector<perl::Object> XandY;
	XandY.push_back(X); XandY.push_back(Y);
	dbgtrace << "Computed list for product..." << endl;
      perl::Object Z = compute_product_complex(XandY);
      
      dbgtrace << "Computing functions" << endl;
      
      //Compute the diagonal functions
      perl::ListResult psi = ListCallPolymakeFunction("atint::diagonal_functions",Xambi);
      
      //Now intersect with the product
      for(int i = 0; i < Xambi; i++) {
	dbgtrace << "Intersecting function" << i+1 << endl;
	Z = divisorByPLF(Z,psi[i]);  
      }
      
      dbgtrace << "Projecting" << endl;
      
      //Finally project
      bool uses_homog = Z.give("USES_HOMOGENEOUS_C");
      Matrix<Rational> raymatrix = Z.give("RAYS");
	raymatrix = raymatrix.minor(All,sequence(0,uses_homog? Xambi+1 : Xambi));
      Matrix<Rational> linmatrix = Z.give("LINEALITY_SPACE");
	if(linmatrix.rows() > 0) {
	  linmatrix = linmatrix.minor(All,sequence(0,uses_homog? Xambi+1 : Xambi));
	}
      IncidenceMatrix<> cones = Z.give("MAXIMAL_CONES");
      Vector<Integer> weights = Z.give("TROPICAL_WEIGHTS");
      perl::Object result("WeightedComplex");
	result.take("RAYS") << raymatrix;
	result.take("LINEALITY_SPACE") << linmatrix,
	result.take("MAXIMAL_CONES") << cones;
	result.take("TROPICAL_WEIGHTS") << weights;
	result.take("USES_HOMOGENEOUS_C") << uses_homog;
      return result;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see perl wrapper
    perl::Object recession_fan(perl::Object complex) {
      bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
	if(!uses_homog) return complex;
      
      //Extract values
      Matrix<Rational> rays = complex.give("RAYS");
      Set<int> directional = complex.give("DIRECTIONAL_RAYS");
      IncidenceMatrix<> cones = complex.give("MAXIMAL_CONES");
      Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
      Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
      int dim = complex.give("CMPLX_DIM");
      std::string desc = complex.give("DESCRIPTION");
      
      Matrix<Rational> newrays = rays.minor(directional,~scalar2set(0));
      Matrix<Rational> newlineality = linspace.minor(All,~scalar2set(0));
      Vector<Integer> newweights;
      
      //Re-map ray indices
      Map<int,int> indexMap; int i = 0;
      for(Entire<Set<int> >::iterator d = entire(directional); !d.at_end(); d++) {
	 indexMap[*d] = i;
	 i++;
      }
      
      //We compute the recession cone of each cell
      Vector<Set<int> > rec_cones;
      for(int mc = 0; mc < cones.rows(); mc++) {
	Set<int> mcDirectional = directional * cones.row(mc);
	//Compute that it has the right dimension
	int mcDim = rank(rays.minor(mcDirectional,All));
	if(mcDirectional.size() > 0 && mcDim == dim) {
	  Set<int> transformCone = attach_operation(mcDirectional,pm::operations::associative_access<Map<int,int>,int>(&indexMap));
	  rec_cones |= transformCone;	  
	  newweights |= weights[mc];
	}
      }
      
      std::ostringstream newdesc;
	newdesc << "Recession fan of \"" << desc << "\"";
            
      //Compute the complexification of the recession cones
      perl::Object cplxify = complexify(newrays,rec_cones,newweights,false);
      //Extract its values and put them into the result
      perl::Object result("WeightedComplex");
	result.take("RAYS") << cplxify.give("RAYS");
	result.take("MAXIMAL_CONES") << cplxify.give("MAXIMAL_CONES");
	result.take("TROPICAL_WEIGHTS") << cplxify.give("TROPICAL_WEIGHTS");
	result.take("LINEALITY_SPACE") << newlineality;
	result.take("DESCRIPTION") << newdesc.str();
      return result;
      
//       Refine fan
//       Matrix<Rational> facetNormals = complex.give("FACET_NORMALS");
//       complex = facetRefinement(complex,facetNormals);
//       
//       //Extract values
//       Matrix<Rational> rays = complex.give("RAYS");
//       Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
//       IncidenceMatrix<> cones = complex.give("MAXIMAL_CONES");
//       Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
//       Set<int> directional = complex.give("DIRECTIONAL_RAYS");
//       int dim = complex.give("CMPLX_DIM"); 
//       
//       
//       //Re-map ray indices
//       Map<int,int> indexMap; int i = 0;
//       for(Entire<Set<int> >::iterator d = entire(directional); !d.at_end(); d++) {
// 	 indexMap[*d] = i;
// 	 i++;
//       }
//       Matrix<Rational> newrays = rays.minor(directional,~scalar2set(0));
//       Matrix<Rational> newlineality = linspace.minor(All,~scalar2set(0));
//       
//       Vector<Set<int> > newcones;
//       Vector<Integer> newweights;
//       
//       //Compute the recession cones of all cones
//       for(int mc = 0; mc < cones.rows(); mc++) {
// 	//Compute the directional rays of the cone
// 	Set<int> mcDirectional = directional * cones.row(mc);
// 	//Compute that it has the right dimension
// 	int mcDim = rank(rays.minor(mcDirectional,All));
// 	if(mcDirectional.size() > 0 && mcDim == dim) {
// 	  //Compute the image of the cone set under the index map
// 	  Set<int> transformCone =
// 	    attach_operation(mcDirectional,pm::operations::associative_access<Map<int,int>,int>(&indexMap));
// 	  //Check if this cone already exists
// 	  bool found = false;
// 	  int tcSize = transformCone.size();
// 	  for(int i = 0; i < newcones.dim(); i++) {
// 	    if((newcones[i]*transformCone).size() == tcSize) {
// 	      found = true;
// 	      newweights[i] += weights[mc];
// 	      break;
// 	    }
// 	  }
// 	  //If the cone doesn't exist yet, add it
// 	  if(!found) {
// 	    newcones = newcones | transformCone;
// 	    newweights |= weights[mc];
// 	  }
// 	}
//       }
//       
//       //Create the result
//       perl::Object result("WeightedComplex");
// 	result.take("RAYS") << newrays;
// 	result.take("MAXIMAL_CONES") << newcones;
// 	result.take("LINEALITY_SPACE") << newlineality;
// 	result.take("TROPICAL_WEIGHTS") << newweights; 
// 	result.take("USES_HOMOGENEOUS_C") << false;
// 	
//       return result;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see perl wrapper
    Integer degree(perl::Object fan) {
      //Extract values
      int dim = fan.give("CMPLX_DIM");
      
      //If it's 0-dimension, just sum up weights
      if(dim == 0) {
	Vector<Integer> weights = fan.give("TROPICAL_WEIGHTS");
	return accumulate(weights,operations::add());
      }
      
      //Otherwise compute the appropriate intersection product
      int ambient_dim = fan.give("CMPLX_AMBIENT_DIM");
      int codim = fan.give("CMPLX_CODIMENSION");
      //Create the function matrix for L^n_codim
      Matrix<Rational> functionmatrix = unit_matrix<Rational>(ambient_dim);
	functionmatrix = zero_vector<Rational>(functionmatrix.cols()) / functionmatrix;
	functionmatrix |= zero_vector<Rational>(functionmatrix.rows());
      perl::Object function("MinMaxFunction");
	function.take("FUNCTION_MATRIX") << functionmatrix;
	function.take("USES_MIN") << false;
      perl::Object div = fan;
      for(int i = 1; i <= ambient_dim - codim; i++) {
	div = divisorByPLF(div,function);
      }
      return degree(div);
    }
    
//     /**
//       @brief Takes a list of tropical fans and checks whether each one has degree 1
//       @param std::vector<perl::Object> fans A list of tropical fans
//       @return int -1, if all are smooth, the index of the first non-smooth local fan otherwise
//     */
//     bool check_smoothness(std::vector<perl::Object> fans) {
//       //Check if all local vertex fans have degree 1
//       for(unsigned int i = 0; i < fans.size(); i++) {
// 	if(degree(fans[i]) != 1) return i;
//       }
//       
//       return -1;
//     }
    
    // ------------------------- PERL WRAPPERS ---------------------------------------------------
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Computes the intersection product of two tropical cycles in a common ambient vector space V."
		      "# @param WeightedComplex X The first tropical variety"
		      "# @param WeightedComplex Y The second tropical variety. Should have the same actual ambient "
		      "# dimension (homogeneous coordinates don't count) as X."
		      "# @return WeightedComplex Z The intersection product X*Y in V = R^n, where n is the ambient "
		      "# dimension of X and Y. The result has homogeneous coordinates, if and only if X or Y has homog. " "# coordinates.",
		      &cycle_intersection,"cycle_intersection(WeightedComplex, WeightedComplex)");

    UserFunction4perl("# @category Tropical geometry"
		      "# Computes the recession fan of a tropical variety"
		      "# @param WeightedComplex complex A tropical variety. If it is a fan, the complex itself is returned"
		      "# @return WeightedComplex A tropical fan, the recession fan of the complex",
		      &recession_fan, "recession_fan(WeightedComplex)");

    UserFunction4perl("# @category Tropical geometry"
		      "# Computes the degree of a tropical variety as the degree of the "
		      "# 0-dimensional complex obtained when intersecting "
		      "# the variety with an appropriate linear space L^n_k"
		      "# @param WeightedComplex complex"
		      "# @return Int",
		      &degree, "degree(WeightedComplex)");
    
//     Function4perl(&check_smoothness,"check_smoothness(;@)");
    
}}