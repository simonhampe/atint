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
 
 This file provides functionality to compute certain special tropical varieties
 */

#include "polymake/client.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/PowerSet.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Array.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Graph.h"
#include "polymake/linalg.h"

namespace polymake { namespace atint {

    using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
    
    //Documentation: see specialvarieties.h
    perl::Object tropical_lnk(const int n, const int k) {
      
      //Ensure that dimensions match
      if(k > n) {
	  throw std::runtime_error("Error creating L^n_k: Fan dimension is larger than ambient dimension");
      }
      if(k < 0 || n < 0) {
	  throw std::runtime_error("Error creating L^n_k: No negative dimensions allowed");
      }
      if(k == 0) {
	  perl::Object result("WeightedComplex");
	  Set<int> nullset;
	    nullset = nullset + 0;
	  Array<Set<int> > nullarray(1);
	    nullarray[0] = nullset;
	  Matrix<Rational> rm(1,n);
	  result.take("RAYS") << rm; 
	  result.take("MAXIMAL_CONES") << nullarray;
	  result.take("USES_HOMOGENEOUS_C") << false;
	  return result;
      }
      
      //Create rays
      Matrix<Rational> rayMatrix(0,n);
      Vector<Rational> e0(n);
      
      for(int i = 1; i <= n; i++) {
	  Vector<Rational> ei(unit_vector<Rational>(n,i-1));
	  rayMatrix = rayMatrix / -ei;
	  e0 = e0 + ei;
      }
      rayMatrix = rayMatrix / e0;
      
      //Create cones
      Set<int> indices;
	for(int i = 0; i <= n; i++) { indices += i;}
      Array<Set<int> >  kSets = pm::Subsets_of_k<Set<int> > ( indices,k );
      
      //Create weights
      Array<int> weights(kSets.size());
      for(int i = 0; i < weights.size(); i++) {
	  weights[i] = 1;
      }
            
      perl::Object fan("WeightedComplex");
	fan.take("RAYS") << rayMatrix;
	fan.take("MAXIMAL_CONES") << kSets;
	fan.take("TROPICAL_WEIGHTS") << weights;
	fan.take("USES_HOMOGENEOUS_C") << false;
	
      return fan;
      
    }
    
    //Documentation see specialvarieties.h
    perl::Object computeBergmanFan(perl::Object fan_skeleton, perl::Object matroid_poly, bool modOutLineality, int projectionCoordinate) {
      //Extract values
      IncidenceMatrix<> maximalCones = fan_skeleton.give("MAXIMAL_CONES");
      Matrix<Rational> facets = matroid_poly.give("FACETS");
      IncidenceMatrix<> facetsThruVertices = matroid_poly.give("FACETS_THRU_VERTICES");
      Matrix<Rational> linearSpan = matroid_poly.give("LINEAR_SPAN");
      Matrix<Rational> rays = fan_skeleton.give("RAYS");
      
      dbgtrace << "Checking for skeleton faces that correspond to loopfree matroids" << endl;
      
      //Compute a list of those n-rank-dimensional faces whose vertices cover [n]
      Vector<Set<int> > listOfFacets;
      for(int mc = 0; mc < maximalCones.rows(); mc++) {
	Vector<Rational> v(rays.cols());
	Set<int> mcSet = maximalCones.row(mc);
	for(Entire<Set<int> >::iterator vindex = entire(mcSet); !vindex.at_end(); vindex++) {
	    v += rays.row(*vindex);
	}
	//Check if the vector has a zero component. If the vertices cover [n], it shouldn't
	bool hasZero = false;
	for(int i = 1; i < v.dim(); i++) { //Start at 1 because of homog. coordinates
	    if(v[i] == 0) {
	      hasZero = true;
	      break;
	    }
	}
	if(!hasZero) {
	    dbglog << "Cone " << mc << " is valid with interior vector " << v << endl;
	    listOfFacets = listOfFacets | maximalCones.row(mc);
	}
      }
      
      dbglog << "Done. " << listOfFacets.dim() << " facets remaining." << endl;
      
      //Now compute normal cones for these faces
      Vector<Set<int> > bergmanCones;
      Vector<Integer > bergmanWeights = ones_vector<Integer>(listOfFacets.dim());
      Matrix<Rational> bergmanLineality = linearSpan.minor(All, ~scalar2set(0));
      
      //For each face: Intersect the set of normal rays of each vertex of the face
      for(int face = 0; face < listOfFacets.dim(); face++) {
	dbgtrace << "computing rays of facet " << face << endl;
	Set<int> raySet = sequence(0,facets.rows());
	dbgtrace << "Starting with " << raySet << endl;
	Set<int> faceSet = listOfFacets[face];
	dbgtrace << "Face has vertices " << faceSet << endl;
	for(Entire<Set<int> >::iterator vertex = entire(faceSet); !vertex.at_end(); ++vertex) {
	  raySet = raySet * facetsThruVertices.row(*vertex);
	}
	dbgtrace << "Remaining rays : " << raySet << endl; 
	//Make this a cone
	bergmanCones = bergmanCones | raySet;
      }
      
      dbgtrace << "Facets are: " << facets.minor(All,~scalar2set(0)) << endl;
      dbgtrace << "Cones are: " << bergmanCones << endl;
      dbgtrace << "Weights are: " << bergmanWeights << endl;
      dbgtrace << "Lineality is: " << bergmanLineality << endl;
      
      Matrix<Rational> bergmanRays = facets.minor(All,~scalar2set(0));
      if(modOutLineality) {
	//Create the projection matrix
	Matrix<Rational> unitMatrix = unit_matrix<Rational>(bergmanRays.cols()-1);
	Matrix<Rational> projectionMatrix(0,unitMatrix.cols());
	
	//Insert a -1's- vector at the right position
	if(projectionCoordinate > 0) {
	    projectionMatrix /= unitMatrix.minor(sequence(0,projectionCoordinate),All);
	}
	projectionMatrix /= - ones_vector<Rational>(unitMatrix.cols());
	if(projectionCoordinate < unitMatrix.rows()) {						
	      projectionMatrix /= unitMatrix.minor(sequence(projectionCoordinate,unitMatrix.rows() -	projectionCoordinate),All);
	}
	
	dbgtrace << "Projection matrix is " << projectionMatrix << endl;
			   
	bergmanRays = bergmanRays * projectionMatrix;
	
	//Apply projection to the lineality space, but make sure the remaining rows are a basis
	bergmanLineality = bergmanLineality * projectionMatrix;
	Set<int> rbasis = basis_rows(bergmanLineality);
	bergmanLineality = bergmanLineality.minor(rbasis,All);
      }
      
      perl::Object result("WeightedComplex");
	result.take("INPUT_RAYS") << bergmanRays;
	result.take("INPUT_CONES") << bergmanCones;
	result.take("TROPICAL_WEIGHTS") << bergmanWeights;
	if(bergmanLineality.rows() > 0) { result.take("LINEALITY_SPACE") << bergmanLineality;}
	
      return result;
    }
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Creates the linear tropical space L^n_k. This tropical fan is defined in the following way: "
		      "# As rays we take -e_i,i=1,...,n, where e_i is the i-th standard basis vector of R^n and "
		      "# e_0 = e_1 + ... + e_n. As maximal cones we take the cones generated by rays {e_i, i in S}, where"
		      "# S runs over all k-element subsets of {0,..,n}."
		      "# @param Int n The ambient dimension of the fan."
		      "# @param Int k The dimension of the fan (should be smaller equal n, otherwise an error is thrown)."
		      "# @return WeightedComplex A PolyhedralFan object representing L^n_k",
		      &tropical_lnk,"tropical_lnk($,$)");       
		      
    Function4perl(&computeBergmanFan,"computeBergmanFan(WeightedComplex, polytope::Polytope,$,$)");
    
}
}