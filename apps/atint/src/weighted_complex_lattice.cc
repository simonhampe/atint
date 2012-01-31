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
 * 
 */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/normalvector.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a polyhedral fan and computes a map of lattice normals. The corresponding property in the fan is set automatically.
   @param WeightedComplex A polyhedral fan, extended by atint to a tropical variety
   */	  
  void computeLatticeNormals(perl::Object fan) {
    
    //Extract basic properties of fan
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    IncidenceMatrix<> codimInc = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    
    //This will contain the result
    Map<int, Map<int, Vector<Integer> > > latticeNormals;
    if(codimOneCones.rows() == 0) {
      fan.take("LATTICE_NORMALS") << latticeNormals;
      return;
    }
    
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
	(latticeNormals[facet])[*e] = latticeNormal(facetmatrix,maxmatrix,additionalRay);
      }
    }
    
    fan.take("LATTICE_NORMALS") << latticeNormals;
  }

  ///////////////////////////////////////////////////////////////////////////////////////

  /**
  @brief Takes a polyhedral fan and computes the weighted sum of the lattice normals at each codimension one face. Sets the corresponding property LATTICE_NORMAL_SUM automatically.
  @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeLatticeNormalSum(perl::Object fan) {
    //Extract all necessary properties
     Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
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

  ///////////////////////////////////////////////////////////////////////////////////////

  /**
  @brief Takes a polyhedral fan and computes whether the fan is balanced as a weighted polyhedral complex. Sets the corresponding property IS_BALANCED automatically.
  @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeIfBalanced(perl::Object fan) {
    //Extract all necessary properties
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    Matrix<Rational> summatrix = fan.give("LATTICE_NORMAL_SUM");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.exists("LINEALITY_SPACE")? 
				  fan.give("LINEALITY_SPACE") :
				  Matrix<Rational>(0,ambient_dim);
    
    //Now iterate over all codim one cones
    for(int facet = 0; facet < codimOneCones.rows(); facet++) {
      Matrix<Rational> vtau = rays.minor(codimOneCones.row(facet),All) / linspace;
      int r = rank(vtau);
      if(rank(vtau/ (summatrix.row(facet))) > r) {
	fan.take("IS_BALANCED") << false;
	dbgtrace << "Not balanced at codimension one face no. " << facet << endl;
	return;
      }
    }
    
    fan.take("IS_BALANCED") << true; 
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
  @brief Takes a polyhedral fan and computes the function vectors for the lattice normals (and their sums). Sets the corresponding properties LATTICE_NORMAL_FCT_VECTOR, LATTICE_NORMAL_SUM_FCT_VECTOR, BALANCED_FACES automatically.
  @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeFunctionVectors(perl::Object fan) {
    //Extract properties from the fan
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    int lineality_dim = fan.give("LINEALITY_DIM");
    Matrix<Rational> linealitySpace = fan.give("LINEALITY_SPACE");
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    
  
    Matrix<Rational> rays = (!uses_homog)? fan.give("RAYS") : fan.give("CMPLX_RAYS");
    Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
    Matrix<Rational> normalsums = fan.give("LATTICE_NORMAL_SUM");
    IncidenceMatrix<> codimOneCones = (!uses_homog)? fan.give("CODIM_1_FACES") : fan.give("CMPLX_CODIM_1_FACES");
    IncidenceMatrix<> maximalCones = (!uses_homog)? fan.give("MAXIMAL_CONES") : fan.give("CMPLX_MAXIMAL_CONES");
    IncidenceMatrix<> coneIncidences = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    
    
    //Result variables
    Map<int, Map<int, Vector<Rational> > > summap;
    Matrix<Rational> summatrix;
    Vector<bool> balancedFaces(codimOneCones.rows());
    
    //Iterate over all codim 1 faces
    for(int fct = 0; fct < codimOneCones.rows(); fct++) {
      dbgtrace << "Facet: " << fct << endl;
      summap[fct] = Map<int, Vector<Rational> >();
      
      Set<int> adjacentCones = coneIncidences.row(fct);
      for(Entire<Set<int> >::iterator mc = entire(adjacentCones); !mc.at_end(); ++mc) {
	dbgtrace << "Maxcone " << *mc << endl;
	Vector<Rational> normalvector((latticeNormals[fct])[*mc]);
	//Compute the representation of the normal vector
	(summap[fct])[*mc]= functionRepresentationVector(maximalCones.row(*mc),normalvector,
				  ambient_dim, uses_homog, rays, linealitySpace, lineality_dim);
      }
      
      //Now compute the representation of the sum of the normals
      try {
	summatrix = summatrix / functionRepresentationVector(codimOneCones.row(fct), normalsums.row(fct),
				  ambient_dim, uses_homog, rays, linealitySpace, lineality_dim);
	balancedFaces[fct] = true;
      }
      catch(std::runtime_error &e) { //This goes wrong, if X is not balanced at a given codim 1 face
	summatrix /= zero_vector<Rational>(rays.rows() + lineality_dim);
	balancedFaces[fct] = false;
      }
    }
    
    //Set fan properties
    fan.take("LATTICE_NORMAL_FCT_VECTOR") << summap;
    fan.take("LATTICE_NORMAL_SUM_FCT_VECTOR") << summatrix; 
    fan.take("BALANCED_FACES") << balancedFaces;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
 /* 
  Vector<Integer> projection_lattice_normal(perl::Object tau, perl::Object sigma, Matrix<Integer> latticeBasis) {
    int dim = sigma.give("CONE_DIM");
    Set<int> I = basis_cols(latticeBasis);
    Matrix<Rational> S = inv(T(latticeBasis.minor(All,I)));
    Matrix<Rational> P(dim,latticeBasis.cols());
    P.minor(All,I) = S;
    
    pm::cout << "Trafo " << P << endl;
    
    Matrix<Rational> taurays = tau.give("RAYS");
      taurays = T(P * T(taurays));
    Matrix<Rational> sigmarays = sigma.give("RAYS");
      sigmarays = T(P * T(sigmarays));
      
    pm::cout << "Tau: " << taurays << endl;
    pm::cout << "sigma: " << sigmarays << endl;
      
    perl::Object trtau("polytope::Cone");
      trtau.take("RAYS") << taurays;
    perl::Object trsigma("polytope::Cone");
      trsigma.take("RAYS") << sigmarays;
      
    Vector<Integer> ln = latticeNormalByCone(trtau,trsigma);
    pm::cout << "ln: " << ln << endl;
    
    return T(latticeBasis)*ln;
    
    
  }*/
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a WeightedComplex and computed the properties LATTICE_GENERATORS and LATTICE_BASES
   */
  void computeLatticeBases(perl::Object complex) {
//     //Extract properties
//     Matrix<Rational> rays = complex.give("RAYS");
//     Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
//     IncidenceMatrix<> cones = complex.give("MAXIMAL_CONES");
//     bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
//     bool is_unimodular = complex.give("IS_UNIMODULAR");
//     Set<int> directional = complex.give("DIRECTIONAL_RAYS");
//     
//     Matrix<Integer> generators;
//     Vector<Set<int> > bases;
//     
//     //Iterate all cones
//     for(int mc = 0; mc < cones.rows(); mc++) {
//       //Compute a lattice basis for the cone
//       
//     }
    
  }
  
  
  
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&computeLatticeNormals, "computeLatticeNormals(WeightedComplex)");

  Function4perl(&computeLatticeNormalSum, "computeLatticeNormalSum(WeightedComplex)");

  Function4perl(&computeIfBalanced, "computeIfBalanced(WeightedComplex)");

  Function4perl(&computeFunctionVectors, "computeFunctionVectors(WeightedComplex)");
  
//   Function4perl(&projection_lattice_normal,"projection_lattice_normal(polytope::Cone,polytope::Cone,Matrix<Integer>)");
  
}}