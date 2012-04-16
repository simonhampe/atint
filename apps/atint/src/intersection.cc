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
#include "polymake/atint/complexify.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/refine.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint {

    using polymake::polytope::cdd_interface::solver;
  
    using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
//     using namespace atintlog::dotrace;
    
    //Documentation see perl wrapper
    perl::Object selective_cycle_intersection(perl::Object X, perl::Object Y) {
      //Extract values
      int Xcodim = X.give("CMPLX_CODIMENSION");
      int Ycodim = Y.give("CMPLX_CODIMENSION");
      int Xambi  = X.give("CMPLX_AMBIENT_DIM");
      
      //dbgtrace << "Checking codimension" << endl;
      
      //If the codimensions of the varieties add up to something larger then CMPLX_AMBIENT_DIM, return the 0-cycle 
      if(Xcodim + Ycodim > Xambi) {
	return CallPolymakeFunction("zero_cycle");
      }
      
      //dbgtrace << "Homogenizing where necessary" << endl;
      
      //Make sure,both are homogeneous
      bool x_uses_homog = X.give("USES_HOMOGENEOUS_C");
      bool y_uses_homog = Y.give("USES_HOMOGENEOUS_C");
      
      if(!x_uses_homog) X = X.CallPolymakeMethod("homogenize");
      if(!y_uses_homog) Y = Y.CallPolymakeMethod("homogenize");
      
      //dbgtrace << "Computing diagonal" << endl;
      
      //Create diagonal equalities and diagonal functions
      Matrix<Rational> diagLin = zero_vector<Rational>() | (zero_vector<Rational>() | 
	  ( - unit_matrix<Rational>(Xambi)) | unit_matrix<Rational>(Xambi));
      	
      perl::ListResult psi = ListCallPolymakeFunction("atint::diagonal_functions",Xambi);
      
      //dbgtrace << "Computing lineality space intersection" << endl;
      
      //The lineality space of the intersection product is the intersection of the lineality spaces
      //Compute the intersection of the two spaces
      Matrix<Rational> x_lineality = X.give("LINEALITY_SPACE");
      Matrix<Rational> y_lineality = Y.give("LINEALITY_SPACE");
      Matrix<Rational> r_lineality(0,Xambi + 1);      
      if(x_lineality.rows() != 0 && y_lineality.rows() != 0){
	//We compute the kernel of (x_lineality | -y_lineality)
 	Matrix<Rational> i_lineality = T(x_lineality  / (-y_lineality));
 	Matrix<Rational> dependence =  null_space(i_lineality);
	r_lineality = dependence.minor(All,sequence(0,x_lineality.rows())) * x_lineality;
	//r_lineality = zero_vector<Rational>() | r_lineality;
      }
      int r_lineality_dim = rank(r_lineality);
      
      //Compute the expected dimension of the intersection product (in the possibly homogeneous coordinates)
      // Substract the lineality dimension, since we compute without it
      int expectedDimension = Xambi - (Xcodim + Ycodim) - r_lineality_dim;
     
      //dbgtrace << "Result: " << r_lineality << endl;
      
      //dbgtrace << "Computing product complex" << endl;
            
      //Compute the product complex
      std::vector<perl::Object> XandY;
	XandY.push_back(X); XandY.push_back(Y);
      perl::Object Z = compute_product_complex(XandY);
      
      //Reduce to those cones that intersect the diagonal at least in the expectedDimension -----------------------
      
      solver<Rational> sv;
      IncidenceMatrix<> maximalCones = Z.give("MAXIMAL_CONES");
      Matrix<Rational> rays = Z.give("RAYS");
      Vector<Integer> weights = Z.give("TROPICAL_WEIGHTS");
      Matrix<Rational> lineality = Z.give("LINEALITY_SPACE");
      IncidenceMatrix<> local_restriction = Z.give("LOCAL_RESTRICTION");
      
      //If we only have lineality spaces, we're done
      if(rays.rows() == 0) {
	perl::Object result("WeightedComplex");
	  result.take("RAYS") << Matrix<Rational>(0,r_lineality.cols());
	  result.take("MAXIMAL_CONES") << maximalCones;
	  result.take("LINEALITY_SPACE") << r_lineality;
	  result.take("TROPICAL_WEIGHTS") << Z.give("TROPICAL_WEIGHTS");
	  result.take("USES_HOMOGENEOUS_C") << true;
	return result;
      }
      
      //Now intersect with the functions -------------------------------------------------------------
      for(int i = 0; i <= Xambi; i++) {
	
	//Remove cones that have too low dimension with the diagonal before computing a divisor
	//dbgtrace << "Reducing to cones that intersect diagonal in appropriate dimension" << endl;
	//dbgtrace << "Have " << maximalCones.rows() << " maximal cones" << endl;
	
	Set<int> remainingCones;
	Set<int> usedRays;
	for(int mc = 0; mc < maximalCones.rows(); mc++) {
	  //Compute intersection with diagonal 
	  std::pair<Matrix<Rational>, Matrix<Rational> > mcFacets = 
	    sv.enumerate_facets(zero_vector<Rational>() | rays.minor(maximalCones.row(mc),All), zero_vector<Rational>() | lineality, true,false);
	    
	    Matrix<Rational> inrays = sv.enumerate_vertices(mcFacets.first, mcFacets.second / diagLin, true,true).first.minor(All,~scalar2set(0));
	  if(rank(inrays) - 1 >= expectedDimension) {
	      remainingCones += mc;
	      usedRays += maximalCones.row(mc);
	  }
	}//END for all maximal cones mc
	//dbgtrace << "Remaining cones: " << remainingCones << endl;
	//dbgtrace << "Used rays: " << usedRays << endl;
	
	rays = rays.minor(usedRays,All);
	maximalCones = maximalCones.minor(remainingCones,usedRays);
	weights = weights.slice(remainingCones);
      	
	//Check for all local restrictions, if there are still any cones containing them
	//If not, discard it
	if(local_restriction.rows() > 0) {
	  //dbgtrace << "Local restriction " << local_restriction << endl;
	  IncidenceMatrix<> new_local_restriction = local_restriction.minor(All,usedRays);
	  Set<int> removable_locals;
	  for(int lr = 0; lr < local_restriction.rows(); lr++) {
	    //If we lose any rays, we lose the cone
	    if(new_local_restriction.row(lr).size() < local_restriction.row(lr).size()) {
	      removable_locals += lr;
	    }
	    else {
	      //Otherwise find a maximal cone containing the row lr
	      bool found_maximal_container = false;
	      for(int mc = 0; mc < maximalCones.rows(); mc++) {
		if( (maximalCones.row(mc) * new_local_restriction.row(lr)).size() == new_local_restriction.row(lr).size()) {
		    found_maximal_container = true; break;
		}
	      }
	      
	      if(!found_maximal_container) removable_locals += lr;
	      
	    }
	  }//END iterate local cones
	  
	  //dbgtrace << "To be removed: " << removable_locals << endl;
	  
	  //If no local restriction remains, return the zero cycle
	  //(Since there is nothing around the local cones)
	  if(removable_locals.size() == local_restriction.rows()) {
	      return CallPolymakeFunction("zero_cycle");
	  }
	  else {
	    local_restriction = new_local_restriction.minor(~removable_locals,All);
	  }
	}//END adapt local cones
	
	
	
	//First of all we check if any cones are left - Otherwise we return the zero cycle
	if(maximalCones.rows() == 0) {
	    return CallPolymakeFunction("zero_cycle");
	}
	
	//Now intersect the current complex with the linearity domains of the current function
	if(i < Xambi) {
	  //dbgtrace << "Computing divisor of function " << (i+1) << endl;
	  perl::Object currentComplex("WeightedComplex");
	    currentComplex.take("RAYS") << rays;
	    currentComplex.take("MAXIMAL_CONES") << maximalCones;
	    currentComplex.take("TROPICAL_WEIGHTS") << weights;
	    currentComplex.take("USES_HOMOGENEOUS_C") << true;
	    currentComplex.take("LINEALITY_SPACE") << lineality;
	    currentComplex.take("LOCAL_RESTRICTION") << local_restriction;
	  perl::Object divisor = divisor_minmax(currentComplex,psi[i],1); 
	  
	  Matrix<Rational> divrays = divisor.give("RAYS");
	  IncidenceMatrix<> divmaximalCones = divisor.give("MAXIMAL_CONES");
	  Vector<Integer> divweights = divisor.give("TROPICAL_WEIGHTS");
	  Matrix<Rational> divlin = divisor.give("LINEALITY_SPACE");
	  IncidenceMatrix<> divlocal = divisor.give("LOCAL_RESTRICTION");
	    rays = divrays; maximalCones = divmaximalCones; weights = divweights; lineality = divlin; local_restriction = divlocal;
	}
	
      }//END for all functions psi[i]
      
      
      
      perl::Object result("WeightedComplex");
 	result.take("RAYS") << rays.minor(All, sequence(0,Xambi+1));
 	result.take("MAXIMAL_CONES") << maximalCones;
 	result.take("TROPICAL_WEIGHTS") << weights;
 	result.take("LINEALITY_SPACE") << r_lineality;
 	result.take("USES_HOMOGENEOUS_C") << true;
       
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
      std::string desc = complex.description();
      
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
	result.set_description() << newdesc.str();
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
      perl::Object div = divisor_minmax(fan,function,ambient_dim-codim);
//       fan;
//       for(int i = 1; i <= ambient_dim - codim; i++) {
// 	div = divisorByPLF(div,function);
//       }
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
		      "# dimension of X and Y. The result has homogeneous coordinates in any case",
		      &selective_cycle_intersection,"cycle_intersection(WeightedComplex, WeightedComplex)");

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
    //Function4perl(&selective_cycle_intersection,"selective_cycle_intersection(WeightedComplex,WeightedComplex)");
    
//     Function4perl(&check_smoothness,"check_smoothness(;@)");
    
}}