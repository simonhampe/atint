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

This file provides (or subsumes) all the functionality necessary to compute
properties of the PolyhedralFan structure extended to a tropical variety in atint.
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/normalvector.h"
#include "polymake/tropical/LoggingPrinter.h"

namespace polymake { namespace fan{ 
  
  //using namespace atint::donotlog;
  //using namespace atint::dolog;
  using namespace atint::dotrace;
  
  /**
    @brief Takes a polyhedral fan and computes its codimension one cones and an incidence matrix indicating which codim one cones lie in which maximal cone. The corresponding properties in the fan are set automatically.
    @param fan::PolyhedralFan fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeCodimensionOne(perl::Object fan) {
    //First we construct the set of all facets 
    Array<IncidenceMatrix<> > maximal_cone_incidence = fan.give("MAXIMAL_CONES_INCIDENCES");
    
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> rays = fan.give("RAYS");
    
    //This will contain the set of indices defining the codim one faces
    Vector<Set<int> > facetArray;
    
    //This will define the codim-1-maximal-cone incidence matrix
    Vector<Set<int> > fIncones;
    
    for(int maxcone = 0; maxcone < maximal_cone_incidence.size(); maxcone++) {
      //This is the incidence matrix for the maximal cone indexed by maxcone
      IncidenceMatrix<> fcts = maximal_cone_incidence[maxcone];
      for(int facet = 0; facet < fcts.rows(); facet++) {
	Set<int> facetToCheck = fcts.row(facet);
	 //If we use homog. coords: Check if this facet intersects x0 = 1, otherwise go to the next one 
	 //More precisely: Check if at least one of its rays has x0-coord != 0
	 if(uses_homog) {
	  Vector<Rational> firstColumn = rays.minor(facetToCheck,All).col(0);
	  if(firstColumn == zero_vector<Rational>(firstColumn.dim())) {
	    continue;
	  }
	 }
	 //Otherwise check if we already have that facet and remember its index
	 int fcIndex = -1;
	 for(int existing = 0; existing < facetArray.dim(); existing++) {
	  if(facetArray[existing] == facetToCheck) {
	    fcIndex = existing;
	    break;
	  }
	 }
	 //Add the facet if necessary and add its maximal-cone indices
	 if(fcIndex == -1) {
	    facetArray = facetArray | facetToCheck;
	    Set<int> singlecone;
	      singlecone = singlecone + maxcone;
	    fIncones = fIncones | singlecone;
	 }
	 else {
	  fIncones[fcIndex] = fIncones[fcIndex] + maxcone;
	 }
      }
    }
    
    fan.take("CODIM_1_FACES") << IncidenceMatrix<>(facetArray);
    fan.take("CODIM_1_IN_MAXIMAL_CONES") << IncidenceMatrix<>(fIncones);
    
  }
  
  /**
   @brief Takes a polyhedral fan and computes a map of lattice normals. The corresponding property in the fan is set automatically.
   @param fan::PolyhedralFan A polyhedral fan, extended by atint to a tropical variety",
   */	  
  void computeLatticeNormals(perl::Object fan) {
    
    //Extract basic properties of fan
    int ambient_dim = fan.give("AMBIENT_DIM");
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    IncidenceMatrix<> codimInc = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    
    //This will contain the result
    Map<int, Map<int, Vector<Integer> > > latticeNormals;
    
    //This equation is added to all cone linear span matrices (and stands for intersecting
    // with (x0 = 1), if we use homogeneous coordinates, to ensure that the lattice normal
    // is of the form (0,...)
    Vector<Rational> intereq(unit_vector<Rational>(ambient_dim,0));
    
    //Compute all the linear spans of the cones before, so we don't do it several times
    Vector<Matrix<Rational> > codimone;
    Vector<Matrix<Rational> > maximal;
    for(int facet = 0; facet < codimOneCones.rows(); facet++) {
      codimone = codimone | null_space(rays.minor(codimOneCones.row(facet),All) / linspace);
    }
    for(int mcone = 0; mcone < maximalCones.rows(); mcone++) {
      maximal = maximal | null_space(rays.minor(maximalCones.row(mcone),All) / linspace);
    }
    
    //Go through all codim one faces
    for(int facet = 0; facet < codimone.dim(); facet++) {
      //Create map for this facet
      latticeNormals[facet] = Map<int, Vector<Integer> >();
      
      //Construct the dual of the linear span of the cone:
      //= kernel of generators
      Matrix<Rational> facetmatrix = codimone[facet];
      if(uses_homog) {
	facetmatrix = facetmatrix / intereq;
      }
      Set<int> adjacentCones = codimInc.row(facet);
      //Go through all adjacent cones
      for(Entire<Set<int> >::iterator e=entire(adjacentCones); !e.at_end(); ++e) {
	Matrix<Rational> maxmatrix = maximal[*e];
	if(uses_homog) {
	  maxmatrix = maxmatrix / intereq;
	}
	//Extract the additional ray by taking any index in the maximal cone
	// and not in the codim 1 cone
	int additionalRayIndex = *(maximalCones.row(*e) - codimOneCones.row(facet)).begin();
	Vector<Rational> additionalRay(rays.row(additionalRayIndex));
	
	//Compute normalvector
	(latticeNormals[facet])[*e] = tropical::latticeNormal(facetmatrix,maxmatrix,additionalRay);
      }
    }
    
    fan.take("LATTICE_NORMALS") << latticeNormals;
  }

  void computeLatticeNormalSum(perl::Object fan) {
    //Extract all necessary properties
     Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
    int ambient_dim = fan.give("AMBIENT_DIM");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Array<Integer> weights = fan.give("TROPICAL_WEIGHTS");
    IncidenceMatrix<> codimInc = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    
    //This will contain the result
    Matrix<Integer> summatrix(0,ambient_dim);
    
    //Iterate over all codim one faces
    for(int facet = 0; facet < codimOneCones.rows(); facet++) {
      //This will contain the weighted sum of the lattice normals
      Vector<Integer> result = zero_vector<Integer>(ambient_dim);
      Set<int> adjacentCones = codimInc.row(facet);
      //Go through all adjacent cones
      for(Entire<Set<int> >::iterator e=entire(adjacentCones); !e.at_end(); ++e) {
	result = result +(latticeNormals[facet])[*e] * weights[*e];
      }
      summatrix = summatrix / result;
    }
     
    fan.take("LATTICE_NORMAL_SUM") << summatrix;
     
  }

  void computeIfBalanced(perl::Object fan) {
    //Extract all necessary properties
    Matrix<Rational> summatrix = fan.give("LATTICE_NORMAL_SUM");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    
    //Now iterate over all codim one cones
    for(int facet = 0; facet < codimOneCones.rows(); facet++) {
      Matrix<Rational> vtau = rays.minor(codimOneCones.row(facet),All) / linspace;
      int r = rank(vtau);
      if(rank(vtau/ (summatrix.row(facet))) > r) {
	fan.take("IS_BALANCED") << false;
	return;
      }
    }
    
    fan.take("IS_BALANCED") << true; 
  }

  
// ------------------------- PERL WRAPPERS ---------------------------------------------------

Function4perl(&computeCodimensionOne,"computeCodimensionOne(fan::PolyhedralFan)");

Function4perl(&computeLatticeNormals, "computeLatticeNormals(fan::PolyhedralFan)");

Function4perl(&computeLatticeNormalSum, "computeLatticeNormalSum(fan::PolyhedralFan)");

Function4perl(&computeIfBalanced, "computeIfBalanced(fan::PolyhedralFan)");
  
}}