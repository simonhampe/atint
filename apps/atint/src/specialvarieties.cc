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
      
      //Reorder rays: Not all rays might be used in cones, so we have to go through all cones
      // and recompute indices
      Vector<int> newIndices(facets.rows()); //contains the new indices of the rays
	for(int i = 0; i < newIndices.dim(); i++) { newIndices[i] = -1;}
      Vector<Set<int> > reorderedCones; //contains the cones in term of the reordered rays
      int nextindex = 0;
      Matrix<Rational> bergmanRays(0,facets.cols());
      for(int c = 0; c < bergmanCones.dim(); c++) {
	Set<int> crays = bergmanCones[c];
	Set<int> newcone;
	//Copy all rays of the cone 
	for(Entire<Set<int> >::iterator r = entire(crays); !r.at_end(); r++) {
	    if(newIndices[*r] == -1) {
	      newIndices[*r] = nextindex;
	      newcone += nextindex;
	      bergmanRays /= facets.row(nextindex);
	      nextindex++;	      
	    }
	    else {
	      newcone += newIndices[*r];
	    }
	}
	reorderedCones |= newcone;
      }
      bergmanCones = reorderedCones;
      if(bergmanRays.rows() > 0) {
	bergmanRays = bergmanRays.minor(All,~scalar2set(0));
      }
      else {
	bergmanRays = Matrix<Rational>(0,facets.cols()-1);
      }
      
      dbgtrace << "Rays are: " << bergmanRays << endl;
      dbgtrace << "Cones are: " << bergmanCones << endl;
      dbgtrace << "Weights are: " << bergmanWeights << endl;
      dbgtrace << "Lineality is: " << bergmanLineality << endl;
      
      if(bergmanRays.rows() == 0 && bergmanLineality.rows() == 0) {
	return perl::Object("WeightedComplex");
      }
      
      if(modOutLineality) {
	int cols = (bergmanRays.rows() > 0? bergmanRays.cols() : bergmanLineality.cols());
	//Create the projection matrix
	Matrix<Rational> unitMatrix = unit_matrix<Rational>(cols-1);
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
			   
	if(bergmanRays.rows() > 0) bergmanRays = bergmanRays * projectionMatrix;
	
	//Apply projection to the lineality space, but make sure the remaining rows are a basis
	if(bergmanLineality.rows() > 0) {
	  bergmanLineality = bergmanLineality * projectionMatrix;
	  Set<int> rbasis = basis_rows(bergmanLineality);
	  bergmanLineality = bergmanLineality.minor(rbasis,All);
	}
      }
      
      perl::Object result("WeightedComplex");
	result.take("RAYS") << bergmanRays;
	result.take("MAXIMAL_CONES") << bergmanCones;
	result.take("USES_HOMOGENEOUS_C") << false;
	result.take("TROPICAL_WEIGHTS") << bergmanWeights;
	if(bergmanLineality.rows() > 0) { result.take("LINEALITY_SPACE") << bergmanLineality;}
	
      return result;
    }
    
    //Documentation see perl wrapper
    perl::Object halfspace_complex(Rational a, Vector<Rational> g ) {
      //Prepare rays and cones
      Matrix<Rational> rays(0,g.dim());
	rays /= Vector<Rational>(g);
	Matrix<Rational> lineality = null_space(rays);
	rays /= (-g);
      Vector<Set<int> > cones;
	Set<int> first; first += 0;
	Set<int> second; second += 1;
      
	
      //If necessary, compute a vertex
      if(a != 0) {
	//Compute the norm squared of g
	Rational sum = accumulate(attach_operation(g,operations::square()),operations::add());
	Vector<Rational> vertex = (a / sum) * g;
	  vertex = 1 | vertex;
	rays = zero_vector<Rational>() | rays;
	rays /= vertex;
	first += 2;
	second += 2;
	lineality = zero_vector<Rational>() | lineality;
      }
      cones |= first;
      cones |= second;
      Vector<Integer> weights = ones_vector<Integer>(2);
      
      //Create result
      perl::Object fan("WeightedComplex");
	fan.take("RAYS") << rays;
	fan.take("MAXIMAL_CONES") << cones;
	fan.take("LINEALITY_SPACE") << lineality;
	fan.take("USES_HOMOGENEOUS_C") << (a != 0);
	fan.take("TROPICAL_WEIGHTS") << weights;
	
      return fan;      
    }
    
    /**
     @brief Computes all vectors of dimension n with entries +1 and -1. They are sorted such that each vector v has the row index determined by the sum: sum_{i: v_i = 1} 2^i (where i runs from 0 to n-1)
     @param int n The column dimension of the matrix
     @return Matrix<Rational> A 2^n by n matrix containing all +-1-vectors of dimension n
     */
    Matrix<Rational> binaryMatrix(int n) {
      Matrix<Rational> result(0,n);
	result /= (- ones_vector<Rational>(n));
      //Now increase the last row of result by "one" in each iteration and append the new row to result
      Integer iterations = pow(2,n)-1;
      int i = 1;
      while(i <= iterations) {
	//Find the first -1-entry
	int index = 0;
	while(result(i-1,index) == 1) { index++;}
	//Now toggle this and all preceding entries
	Vector<Rational> newrow(result.row(i-1));
	  newrow[index] = 1;
	  for(int j = 0; j < index; j++) newrow[j] = -1;
	result /= newrow;
	i++;
      }
      return result;
    }
    
    /**
     @brief Assumes v is a vector with entries +1 and -1 only. Returns sum_{i: v_i = 1} 2^i (where i runs from 0 to n-1
     */
    inline int binaryIndex(Vector<Rational> v) {
      int result = 0;
      for(int i = 0; i < v.dim(); i++) {
	if(v[i] == 1) result += pow(2,i);
      }
      return result;
    }
    
    //Documentation see perl wrapper
    perl::Object tropical_cube(int n, int k) {
      //Create the cube vertices
      Matrix<Rational> rays = binaryMatrix(n);
      Vector<Set<int> > cones;
      
      //First we treat the special case where k == 0 (or n < k)
      if(k == 0 || n< k) {
	perl::Object result("WeightedComplex");
	  result.take("RAYS") << rays;
	    Vector<Set<int> > singlefaces;
	     for(int i = 0; i < rays.rows(); i++) {
	      Set<int> iset; iset += i;
	      singlefaces |= iset;
	     }
	  result.take("MAXIMAL_CONES") << singlefaces;
	  result.take("TROPICAL_WEIGHTS") << ones_vector<Integer>(rays.rows());
	  result.take("USES_HOMOGENEOUS_C") << true;
	return result;
      }

      //Now create the k-skeleton of the n-cube: For each n-k-set S of 0,..,n-1 and for each vertex
      // v of the n-k-dimensional cube: Insert the entries of v in S and then insert all possible 
      //vertices of the k-dimensional cube in S^c to obtain a k-dimensional face of the cube
      Array<Set<int> > nmkSets = pm::Subsets_of_k<Set<int> > ( sequence(0,n),n-k );
      Matrix<Rational> nmkVertices = binaryMatrix(n-k);
      Matrix<Rational> kVertices = binaryMatrix(k);
      
      for(int s = 0; s < nmkSets.size(); s++) {
	for(int v = 0; v < nmkVertices.rows(); v++) {
	   Set<int> S = nmkSets[s];
	   Set<int> newface;
	   Vector<Rational> vertex(n);
	   vertex.slice(S) = nmkVertices.row(v);
	   for(int w = 0; w < kVertices.rows(); w++) {
	      vertex.slice(~S) = kVertices.row(w);
	      newface += binaryIndex(vertex);
	   }
	   cones |= newface;
	}
      }//End create k-skeleton
      
      int vertexnumber = rays.rows();
      
      //Now we also create the k-1-skeleton of the cube to compute the ray faces
      Array<Set<int> > nmlSets = pm::Subsets_of_k<Set<int> > (sequence(0,n),n-k+1);
      Matrix<Rational> nmlVertices = binaryMatrix(n-k+1);
      Matrix<Rational> lVertices = binaryMatrix(k-1);
      Vector<Set<int> > raycones;
      
      for(int s = 0; s < nmlSets.size(); s++) {
	for(int v = 0; v < nmlVertices.rows(); v++) {
	    Set<int> S = nmlSets[s];
	    Set<int> newface;
	    Vector<Rational> vertex(n);
	    vertex.slice(S) = nmlVertices.row(v);
	    for(int w = 0; w < lVertices.rows(); w++) {
	      vertex.slice(~S) = lVertices.row(w);
	      newface += binaryIndex(vertex);
	    }
	    raycones |= newface;
	}
      }//End create k-1-skeleton
      
      
      //We add a copy of each vertex and consider it a ray.
      //Now, for each face S of the k-1-skeleton, we add a cone that contains for each i in S:
      //The vertex i and its corresponding ray
      if(k > 0) rays = rays / rays;
      int iter = raycones.size();
      for(int c = 0; c < iter; c++) {
	Set<int> newface; 
	Set<int> cubeface = raycones[c];
	for(Entire<Set<int> >::iterator v = entire(cubeface); !v.at_end(); v++) {
	    newface += *v;
	    int rayindex = *v + vertexnumber;
	    newface += rayindex;
	}
	cones |= newface;
      }
      
      //Now create the result
      Vector<Rational> homog_coord = ones_vector<Rational>(vertexnumber) |
				    zero_vector<Rational>(vertexnumber);
      rays = homog_coord | rays; //Homog. Coords!
      Vector<Integer> weights = ones_vector<Integer>(cones.size());
      perl::Object result("WeightedComplex"); 
	result.take("RAYS") << rays;
	result.take("MAXIMAL_CONES") << cones;
	result.take("TROPICAL_WEIGHTS") << weights;
	result.take("USES_HOMOGENEOUS_C") << true;
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
	
    UserFunction4perl("# @category Tropical geometry"
		      "# Creates the halfspace complex defined by an rational vector g and a rational b, i.e. the "
		      "# complex consisting of the two maximal cones g >= a and g <= a"
		      "# @param Rational constant The constant translation a"
		      "# @param Vector<Rational> equation The defining equation g"
		      "# @return WeightedComplex The resultin halfspace complex",
		      &halfspace_complex,"halfspace_complex($,Vector<Rational>)");

    UserFunction4perl("# @category Tropical geometry"
		      "# Creates the tropical cube T^n_k, i.e. the tropical variety obtained by"
		      "# glueing together the 2^n L^n_k obtained by applying all possible sign "
		      "# changes, such that the vertices form the k-skeleton of "
		      "# the n-dimensional cube"
		      "# @param Int n The ambient dimension"
		      "# @param Int k The dimension of the cube"
		      "# @return WeightedComplex cube",
		      &tropical_cube,"tropical_cube($,$)");
    
    Function4perl(&computeBergmanFan,"computeBergmanFan(WeightedComplex, polytope::Polytope,$,$)");
    Function4perl(&binaryMatrix,"binaryMatrix($)");
}
}