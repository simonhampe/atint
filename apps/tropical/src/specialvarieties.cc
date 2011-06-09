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
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/PowerSet.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Array.h"
#include "polymake/IncidenceMatrix.h"


namespace polymake { namespace tropical {

    //using namespace atint::donotlog;
    //using namespace atint::dolog;
    using namespace atint::dotrace;
    
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
	  perl::Object result("fan::PolyhedralFan");
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
            
      perl::Object fan("fan::PolyhedralFan");
	fan.take("RAYS") << rayMatrix;
	fan.take("MAXIMAL_CONES") << kSets;
	fan.take("TROPICAL_WEIGHTS") << weights;
	fan.take("USES_HOMOGENEOUS_C") << false;
	
      return fan;
      
    }
    
    //Documentation see specialvarieties.h
    perl::Object computeBergmanFan(perl::Object fan_skeleton, perl::Object matroid_poly, bool modOutLineality) {
      //Extract values
      IncidenceMatrix<> maximalCones = fan_skeleton.give("MAXIMAL_CONES");
      Matrix<Rational> facets = matroid_poly.give("FACETS");
      IncidenceMatrix<> facetsThruVertices = matroid_poly.give("FACETS_THRU_VERTICES");
      Matrix<Rational> linearSpan = matroid_poly.give("LINEAR_SPAN");
      Matrix<Rational> rays = fan_skeleton.give("RAYS");
      
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
	for(int i = 1; i < v.dim(); i++) {
	    if(v[i] == 0) {
	      hasZero = true;
	      break;
	    }
	}
	if(!hasZero) {
	    listOfFacets = listOfFacets | maximalCones.row(mc);
	}
      }
      
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
	dbgtrace << "Face has vertices "
	for(Entire<Set<int> >::iterator vertex = entire(faceSet); !vertex.at_end(); ++vertex) {
	  raySet = raySet * facetsThruVertices.row(*vertex);
	}
	//Make this a cone
	bergmanCones = bergmanCones | raySet;
      }
      
      dbgtrace << "Facets are: " << facets.minor(All,~scalar2set(0)) << endl;
      dbgtrace << "Cones are: " << bergmanCones << endl;
      dbgtrace << "Weights are: " << bergmanWeights << endl;
      dbgtrace << "Lineality is: " << bergmanLineality << endl;
      
      //TODO: Mod out lineality
      
      perl::Object result("fan::PolyhedralFan");
	result.take("INPUT_RAYS") << facets.minor(All,~scalar2set(0));
	result.take("INPUT_CONES") << bergmanCones;
	result.take("TROPICAL_WEIGHTS") << bergmanWeights;
	result.take("LINEALITY_SPACE") << bergmanLineality;
	
      return result;
    }
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Creates the linear tropical space L^n_k. This tropical fan is defined in the following way: "
		      "# As rays we take -e_i,i=1,...,n, where e_i is the i-th standard basis vector of R^n and "
		      "# e_0 = e_1 + ... + e_n. As maximal cones we take the cones generated by rays {e_i, i in S}, where"
		      "# S runs over all k-element subsets of {0,..,n}."
		      "# @param Int n The ambient dimension of the fan."
		      "# @param Int k The dimension of the fan (should be smaller equal n, otherwise an error is thrown)."
		      "# @return fan::PolyhedralFan A PolyhedralFan object representing L^n_k",
		      &tropical_lnk,"tropical_lnk($,$)");       
		      
    Function4perl(&computeBergmanFan,"computeBergmanFan(fan::PolyhedralFan, polytope::Polytope,$)");
    
//     UserFunction4perl("# @category Tropical geometry"
// 		      "# Creates the bergman fan of a given matroid fan."
// 		      "# @param matroid::Matroid m A matroid"
// 		      "# @param Bool modOutLineality Optional argument. If set to TRUE, the lineality space is divided out before returning the "
// 		      "# fan. The next parameter specifies the exact modalities of the division. By default, this parameter is set to FALSE"
// 		      "# @param int projectionCoordinate Optional argument. An integer in {0,..,n-1}, where n is the number of elements of the matroid. If modOutLineality is set to TRUE, the standard basis vector with index projectionCoordinate is mapped to minus the sum of the remaining standard basis vectors to mod out the lineality space. By default, this is 0."
// 		      "# @return fan::PolyhedralFan The bergman fan of the matroid, possibly with the lineality space divided out",
// 		      &bergman_fan,"bergman_fan($;$=0,$=0)");

  //Function4perl(&bergman_fan_via_polytope,"bergman_fan_via_polytope(polytope::Polytope,$;$=0,$=0)");
}
}