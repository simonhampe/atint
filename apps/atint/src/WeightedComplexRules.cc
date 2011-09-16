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
#include "polymake/AccurateFloat.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/normalvector.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
  
  using polymake::polytope::cdd_interface::solver;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Takes a polyhedral fan and computes its codimension one cones and an incidence matrix indicating which codim one cones lie in which maximal cone. The corresponding properties in the fan are set automatically.
    @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeCodimensionOne(perl::Object fan) {
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> rays = fan.give("RAYS");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    
    dbgtrace << "Computing all facets..." << endl;
    
    //First we construct the set of all facets 
    //Array<IncidenceMatrix<> > maximal_cone_incidence = fan.give("MAXIMAL_CONES_INCIDENCES");
    //Compute the rays-in-facets for each cone directly
    Vector<IncidenceMatrix<> > maximal_cone_incidence;
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      Set<int> mset = maximalCones.row(mc);
      //Extract inequalities
      dbgtrace << "Computing facets for cone set " << mset << endl;
      Matrix<Rational> facets = solver<Rational>().enumerate_facets(
	    zero_vector<Rational>() | rays.minor(mset,All),
	    zero_vector<Rational>() | linspace).first.minor(All, ~scalar2set(0));;
      dbgtrace << "Done. Checking rays..." << endl;
      //For each inequality, check which rays lie in it
      Vector<Set<int> > facetIncidences;
      for(int row = 0; row < facets.rows(); row++) {
	Set<int> facetRays;
	for(Entire<Set<int> >::iterator m = entire(mset); !m.at_end(); m++) {
	  if(facets.row(row) * rays.row(*m) == 0) {
	    facetRays += *m;
	  }
	}
	facetIncidences |= facetRays;
      }
      dbgtrace << "Done." << endl;
      maximal_cone_incidence |= IncidenceMatrix<>(facetIncidences);
    }
    
    dbgtrace << "Check for doubles and useless facets..." << endl;
    
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
  
//  /**
//*    @brief Takes a polyhedral fan and computes its codimension one cones and an incidence matrix indicating which codim one cones lie in which maximal cone. The corresponding properties in the fan are set automatically. This function differs from the function computeCodimensionOne in that it does not use the property MAXIMAL_CONES_INCICDENCES, but uses the FACETS_THRU_VERTICES property of polytope::Cone (which seems to make everything a lot faster)
//     @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
//     @param std::vector<perl::Object> cones An array of polytope::Cone objects (whose RAYS_IN_FACETS should have been precomputed) that represent the maximal cones of fan (s.t. the i-th cone is the i-th maximal cone)*/
//   */
//   void computeCodimensionOneViaCones(perl::Object fan, std::vector<perl::Object> cones) {
//     
//     //Extract needed properties
//     Matrix<Rational> rays = fan.give("RAYS");
//     IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
//     bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
//     
//     //This will contain the array of codim-1-cones, which we use to construct the incidence matrix
//     Vector<Set<int> > facetArray;
//     //This will define the codim-1-maximal-cone incidence matrix
//     Vector<Set<int> > facetsInCones;
//     
//     dbgtrace << "Going through all maximal cones" << endl;
//     
//     //Iterate through all cone objects
//     for(int i = 0; i < maximalCones.rows(); i++) {
//       dbgtrace << "Considering maximal cone no. " << i << endl;
//       //Put its rays in a vector
//       Vector<int> sortedRayIndices(maximalCones.row(i));
//       //Retrieve the facets of the cone
//       IncidenceMatrix<> raysInFacets = cones[i].give("RAYS_IN_FACETS");
//       //Now go through all facets of the cone and map them back to facets of the fan
//       for(int irow = 0; irow < raysInFacets.rows(); irow++) {
// 	dbgtrace << "Considering its facet no. " << irow << endl;
// 	//Convert ray indices back to fan indices
// 	Set<int> potentialFacet;
// 	Set<int> coneRays = raysInFacets.row(irow);
// 	for(Entire<Set<int> >::iterator iray = entire(coneRays); !iray.at_end(); iray++) {
// 	    potentialFacet = potentialFacet + sortedRayIndices[*iray];
// 	}
// 	dbgtrace << "Converted indices. Now checking for existence" << endl;
// 	//If we use homog. coords: Check if this facet intersects x0 = 1, otherwise go to the next one 
// 	//More precisely: Check if at least one of its rays has x0-coord != 0
// 	if(uses_homog) {
// 	  Vector<Rational> firstColumn = rays.minor(potentialFacet,All).col(0);
// 	  if(firstColumn == zero_vector<Rational>(firstColumn.dim())) {
// 	    continue;
// 	  }
// 	}
// 	//Otherwise check if we already have that facet and remember its index
// 	int fcIndex = -1;
// 	for(int existing = 0; existing < facetArray.dim(); existing++) {
// 	if(facetArray[existing] == potentialFacet) {
// 	  fcIndex = existing;
// 	  break;
// 	}
// 	}
// 	//Add the facet if necessary and add its maximal-cone indices
// 	if(fcIndex == -1) {
// 	  facetArray = facetArray | potentialFacet;
// 	  Set<int> singlecone;
// 	    singlecone = singlecone + i;
// 	  facetsInCones = facetsInCones | singlecone;
// 	}
// 	else {
// 	  facetsInCones[fcIndex] = facetsInCones[fcIndex] + i;
// 	}
// 	
// 	
//       }
//       
//     }
    
//     //Finally: Insert values
//     fan.take("CODIM_1_FACES") << IncidenceMatrix<>(facetArray);
//     fan.take("CODIM_1_IN_MAXIMAL_CONES") << IncidenceMatrix<>(facetsInCones);    
//   }
  
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
  @brief  This method takes a set of row indices for [[RAYS]] (or [[CMPLX_RAYS]] in the homogeneous case) and a vector that is supposed to be in the span of these row vectors and the lineality space (or their affine space, see description of [[LATTICE_NORMAL_FCT_VECTOR]]). It then computes the corresponding represenation in these vectors
  @param s a set of row indices of [[RAYS]] (or [[CMPLX_RAYS]])
  @param v a vector supposed to lie in the span (or affine span) of [[RAYS]] (or [[CMPLX_RAYS]]) + [[LINEALITY_SPACE]]
  @param ambient_dim The ambient dimension of the fan
  @param uses_homog Whether the fan uses homogeneous coordinates
  @param rays The matrix of [[RAYS]] or [[CMPLX_RAYS]]
  @param linealitySpace A matrix of generators of the lineality space
  @param lineality_dim The dimension of the lineality space
  @return A vector of length [[N_RAYS]] + [[LINEALITY_DIM]] (or [[CMPLX_RAYS]]->rows() + [[LINEALITY_DIM]] in the homogeneous case) with linear coefficients of a representation in the generators chosen via s. The last elements always refer to the lineality space.
  */
  Vector<Rational> functionRepresentationVector(const Set<int> &rayIndices, const Vector<Rational> &v,
						int ambient_dim, bool uses_homog, 
						const Matrix<Rational> &rays,
						const Matrix<Rational> &linealitySpace,
						int lineality_dim) {
    dbgtrace << "Starting representation computation" << endl;
    //Put ray indices in fixed order
    Array<int> fixedIndices(rayIndices);
    //Matrix of generators
    Matrix<Rational> m(0,ambient_dim);
    //First affine ray (not used in non-homog. coords)
    Vector<Rational> baseray;
    int baseRayIndex = -1; //Index of baseray in fixedIndices
    
    //Compute matrix of generators
    if(!uses_homog) {
      m = m / rays.minor(rayIndices,All);
    }
    else {
      for(int r = 0; r < fixedIndices.size(); r++) {
	Vector<Rational> rayVector = rays.row(fixedIndices[r]);
	//Add all directional rays as is
	if(rayVector[0] == 0) {
	  m = m / rayVector;
	}
	//Use relative differences in the homog. case
	else {
	  if(baseRayIndex == -1) {
	    baseray = rayVector;
	    baseRayIndex = r;
	  }
	  else {
	    m = m / (rayVector - baseray);
	  }
	}
      }
    }
    if(lineality_dim > 0) {
      m = m / linealitySpace;
    }
    
    //If there were no row indices, i.e. m = 0-space, just enter a zero row
    if(m.rows() == 0) {
      m = m / zero_vector<Rational>(ambient_dim);
    }
    
    dbgtrace << "Generator matrix is " << m << endl;
    dbgtrace << "Vector is " << v << endl;
    
    //Now compute the representation
    Vector<Rational> repv = linearRepresentation(v,m);
    
    dbgtrace << "Representation vector: " << repv << endl;
    
    if(repv.dim() == 0) {
      throw std::runtime_error("Error: vector not in linear span of generators");
    }
    
    //Insert coefficients at correct places
    Vector<Rational> result(lineality_dim + rays.rows());
    for(int r = 0; r < fixedIndices.size(); r++) {
      if(r != baseRayIndex) {
	//If a ray came after the baseray, its matrix row index is one lower then its array index.
	int matrixindex = (baseRayIndex == -1)? r : (r > baseRayIndex? r-1 : r);
	result[fixedIndices[r]] = repv[matrixindex];
	dbgtrace << "Inserting " << repv[matrixindex] << " at " << fixedIndices[r] << endl;
	//if this is an affine ray, substract its coefficient at the baseray
	if(rays(fixedIndices[r],0) != 0 && uses_homog) {
	  result[fixedIndices[baseRayIndex]] -= repv[matrixindex];
	}
      }
    }
    
    dbgtrace << "Result vector is " << result << endl;
    
    //Insert linspace coefficients at the end
    int repvSize = repv.dim();
    for(int lingen = 0; lingen < lineality_dim; lingen++) {
     result[rays.rows() + lingen] = repv[repvSize - lineality_dim + lingen];
    }
    dbgtrace << "Done." << endl;
    return result;    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
  @brief Takes a polyhedral fan and computes the function vectors for the lattice normals (and their sums). Sets the corresponding properties LATTICE_NORMAL_FCT_VECTOR, LATTICE_NORMAL_SUM_FCT_VECTOR automatically.
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
      summatrix = summatrix / functionRepresentationVector(codimOneCones.row(fct), normalsums.row(fct),
				  ambient_dim, uses_homog, rays, linealitySpace, lineality_dim);
    }
    
    //Set fan properties
    fan.take("LATTICE_NORMAL_FCT_VECTOR") << summap;
    fan.take("LATTICE_NORMAL_SUM_FCT_VECTOR") << summatrix; 
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Takes a polyhedral fan and computes the ray data of the corresponding polyhedral complex. Sets the corresponding properties CMPLX_RAYS, CMPLX_MAXIMAL_CONES, CMPLX_CODIM_1_FACES automatically.
    @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeComplexData(perl::Object fan) {
    //Extract properties of fan
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    
    if(!uses_homog) {
      fan.take("CMPLX_RAYS") << fan.give("RAYS");
      fan.take("CMPLX_MAXIMAL_CONES") << fan.give("MAXIMAL_CONES");
      fan.take("CMPLX_CODIM_1_FACES") << fan.give("CODIM_1_FACES");
      return;
    }
    
    Matrix<Rational> rays = fan.give("RAYS");
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    IncidenceMatrix<> facet_incidences = fan.give("CODIM_1_IN_MAXIMAL_CONES");
      facet_incidences = T(facet_incidences);
      
    //Result variables
    Matrix<Rational> cmplxrays(0,ambient_dim);
    Vector<Set<int> > maxcones(maximalCones.rows());
    Vector<Set<int> > codimone(codimOneCones.rows());
    
    dbgtrace << "Dividing rays..." << endl;
    
    //Divide the set of rays into those with x0 != 0 and those with x0 = 0
    Set<int> affineRays;
    Set<int> directionalRays;
    Map<int,int> newAffineIndices; //This maps the old ray indices to the new ones in cmplxrays
    for(int r = 0; r < rays.rows(); r++) {
      if(rays.row(r)[0] == 0) {
	directionalRays = directionalRays + r;
      }
      else {
	affineRays = affineRays + r;
	cmplxrays = cmplxrays / rays.row(r);
	newAffineIndices[r] = cmplxrays.rows()-1;
      }
    }
    
    dbgtrace << "Affine rays: " << affineRays << ", directional rays: " << directionalRays << endl;
    
    //Insert the indices of the new affine rays for each cone
    for(int co = 0; co < codimOneCones.rows(); co++) {
      Set<int> corays = codimOneCones.row(co) * affineRays;
      codimone[co] = Set<int>();
      for(Entire<Set<int> >::iterator e = entire( corays); !e.at_end(); ++e) {
	codimone[co] = codimone[co] + newAffineIndices[*e];
      }
    }
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      Set<int> mcrays = maximalCones.row(mc) * affineRays;
      maxcones[mc] = Set<int>();
      for(Entire<Set<int> >::iterator e = entire( mcrays); !e.at_end(); ++e) {
	maxcones[mc] = maxcones[mc] + newAffineIndices[*e];
      }
    }
    
    dbgtrace << "Added affine rays to cones" << endl;
    
    //Now we go through the directional rays and compute the connected component for each one
    for(Entire<Set<int> >::iterator r = entire(directionalRays); !r.at_end(); ++r) {
      
      dbgtrace << "Computing components of ray " << *r << endl;
      
      //List of connected components of this ray, each element is a component
      //containing the indices of the maximal cones
      Vector<Set<int> > connectedComponents;
      //The inverse of the component matrix, i.e. maps cone indices to row indices of connectedComponents
      Map<int,int> inverseMap;
      
      //Compute the set of maximal cones containing r
      Set<int> rcones;
      for(int mc = 0; mc < maximalCones.rows(); mc++) {
	if(maximalCones.row(mc).contains(*r)) {
	  rcones = rcones + mc;
	}
      }
      
      dbgtrace << "Computed set of cones containing r:" << rcones << endl;
      
      //For each such maximal cone, compute its component (if it hasnt been computed yet).
      for(Entire<Set<int> >::iterator mc = entire(rcones); !mc.at_end(); ++mc) {
	if(!inverseMap.exists(*mc)) {
	  dbgtrace << "Creating new component" << endl;
	  //Create new component
	  Set<int> newset; newset = newset + *mc;
	  connectedComponents = connectedComponents | newset;
	  inverseMap[*mc] = connectedComponents.dim()-1;
	  
	  //Do a breadth-first search for all other cones in the component
	  std::list<int> queue;
	    queue.push_back(*mc);
	    //Semantics: Elements in that queue have been added but their neighbours might not
	    dbgtrace << "Calculating component" << endl;
	  while(queue.size() != 0) {
	    int node = queue.front(); //Take the first element and find its neighbours
	      queue.pop_front();
	    for(Entire<Set<int> >::iterator othercone = entire(rcones); !othercone.at_end(); othercone++) {
	      //We only want 'homeless' cones
	      if(!inverseMap.exists(*othercone)) {
		//This checks whether both cones share a ray with x0=1
		if((maximalCones.row(node) * maximalCones.row(*othercone) * affineRays).size() > 0) {
		    //Add this cone to the component
		    connectedComponents[connectedComponents.dim()-1] += *othercone;
		    inverseMap[*othercone] = connectedComponents.dim()-1;
		    queue.push_back(*othercone);
		}
	      }
	    }
	  }
	  
	}
      } //END computation of connected components
      
      dbgtrace << "Connected components:\n" << connectedComponents << endl;
      
      //Now add r once for each connected component to the appropriate cones
      for(int cc = 0; cc < connectedComponents.dim(); cc++) {
	cmplxrays = cmplxrays / rays.row(*r);
	int rowindex = cmplxrays.rows()-1;
	Set<int> ccset = connectedComponents[cc];
	dbgtrace << "Inserting for component " << cc+1 << endl;
	for(Entire<Set<int> >::iterator mc = entire(ccset); !mc.at_end(); ++mc) {
	  maxcones[*mc] = maxcones[*mc] + rowindex;
	  //For each facet of mc that contains r, add rowindex
	  Set<int> fcset;
	  //If there are maximal cones not intersecting x0 = 1, they have no facets
	  //in facet_incidences, hence the following check
	  if(*mc < facet_incidences.rows()) {
	    fcset = facet_incidences.row(*mc);
	  }
	  for(Entire<Set<int> >::iterator fct = entire(fcset); !fct.at_end(); ++fct) {
	    if(codimOneCones.row(*fct).contains(*r)) {
	      codimone[*fct] = codimone[*fct] + rowindex;
	    }
	  }
	}
      }
      
    }//END iterate over all rays
    
    dbgtrace << "Done computing rays, inserting values..." << endl;
    
    //Insert values
    fan.take("CMPLX_RAYS") << cmplxrays;
    fan.take("CMPLX_MAXIMAL_CONES") << maxcones;
    fan.take("CMPLX_CODIM_1_FACES") << codimone;
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the data necessary to visualize the polyhedral complex. Will be called by CMPLX_VISUAL.
   @param fan The weighted polyhedral complex (in homogeneous coordinates!)
   @param scale The scaling factor for the directional rays that are added to the affine rays
   @param showWeights Whether the weights of the cells should be displayed in the center
   @return A perl array containg:
   A list of rational polytopes representing the complex.
   For each cell of the complex, there is a polytope obtained by adding all directional rays (with a scaling factor) to all
   affine rays.
   2)  A polytope::PointConfiguration that will contain the center of each cell as vertex, labelled with the corresponding weight. This is only computed if showWeights is true.
  */
  perl::ListReturn computeVisualPolyhedra(const perl::Object &fan, const Rational &scale, bool showWeights) { 
    //Extract values
    bool weightsExist = fan.exists("TROPICAL_WEIGHTS");
    Array<Integer> weights;
      if(weightsExist) {
	weights = fan.give("TROPICAL_WEIGHTS");
      }
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linealitySpace = fan.give("LINEALITY_SPACE");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    
    dbgtrace << "Extracted values" << endl;
    
    //First separate affine and directional rays
    Set<int> affineRays;
    Set<int> directionalRays;
    for(int r = 0; r < rays.rows(); r++) {
      if(rays.row(r)[0] == 0) {
	directionalRays = directionalRays + r;
      }
      else {
	affineRays = affineRays + r;
      }
    }
    
    dbgtrace << "Separated rays" << endl;
    
    //Create a polytope for each cone
    perl::ListReturn result;
    
    //This will contain the cell centers with the weight labels
    perl::Object weightCenters("polytope::PointConfiguration");
    Matrix<Rational> centermatrix(0,ambient_dim);
    Vector<std::string> centerlabels;
      
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      
      dbgtrace << "Computing geometry for cone " << mc << endl;
      
      //First, create a point matrix by adding all directional rays to all affine rays
      Matrix<Rational> v(0,ambient_dim);
      
      Set<int> maxAffine = maximalCones.row(mc) * affineRays;
      Set<int> maxDirectional = maximalCones.row(mc) * directionalRays;
      
      for(Entire<Set<int> >::iterator aRay = entire(maxAffine); !aRay.at_end(); ++aRay) {
	v = v / rays.row(*aRay);
	for(Entire<Set<int> >::iterator dRay = entire(maxDirectional); !dRay.at_end(); ++dRay) {
	    v = v / (rays.row(*aRay) + scale * rays.row(*dRay));
	}
	for(int linrow = 0; linrow < linealitySpace.rows(); linrow++) {
	    v = v / (rays.row(*aRay) + scale * linealitySpace.row(linrow));
	    v = v / (rays.row(*aRay) - scale * linealitySpace.row(linrow));
	}
      }
      
      //Then create a rational polytope for labelling
      dbgtrace << "Points " << v << endl;
      perl::Object ratPolytope("polytope::Polytope<Rational>");
	ratPolytope.take("POINTS") << v;
      
      result << ratPolytope;
            
      //If necessary, compute centers and label with weights
      if(showWeights && weightsExist) {
	Vector<Rational> center = ratPolytope.give("VERTEX_BARYCENTER");
	centermatrix = centermatrix / center;
	std::ostringstream wlabel;
	wlabel << "# " << mc << ": " << weights[mc];
	centerlabels = centerlabels | wlabel.str();
      }      
      
      
      dbgtrace << "Done." << endl;
      
    }
    
    if(showWeights && weightsExist) {
      dbgtrace << "Computed weight labels, inserting them" << endl;
      weightCenters.take("POINTS") << centermatrix;
      weightCenters.take("LABELS") << centerlabels;
      dbgtrace << "Done" << endl;
    }
    
    result << weightCenters;
     
    dbgtrace << "Done." << endl;
      
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Computes the polyhedral data necessary for visualization with a bounding box
    @param perl::Object fan The polyhedral complex to be visualized, in homogeneous coordinates
    @param bool isRelative true, iff the bounding box is given relative to the complex
    @param bool showWeights If true, the barycenters of the polytopes are computed for weight labelling
    @param Rational bbDistance The relative distance of the border of the bounding box to the affine part of the complex (should be a positive number).
    @param bool onlyBoundingBox If true, only the relative bounding box with respect to the given bbDistance is computed and returned
    @param Matrix<Rational> bBox The absolute bounding box needed if isRelative is false. Given by two row vectors indicating the extreme points of the box
    @return A perl::ListReturn containing 
    1) the list of polytopes to be rendered
    2) A polytope::PointConfiguration that will contain the center of each cell as vertex, labelled with the corresponding weight. This is only computed if showWeights is true, but is contained in the ListReturn in any case.
    If however, onlyBoundingBox is true, the ListReturn will only contain a Matrix<Rational> specifying the relative  bounding box.
  */
  perl::ListReturn computeBoundedVisual(perl::Object fan, bool isRelative, bool showWeights,Rational bbDistance, bool onlyBoundingBox, Matrix<Rational> bBox) {
    
    //Extract values
    int ambient_dim = fan.give("CMPLX_AMBIENT_DIM");
    Matrix<Rational> rays = fan.give("RAYS");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    Matrix<Rational> facetNormals = fan.give("FACET_NORMALS");
    Matrix<Rational> facetNormalsInCones = fan.give("MAXIMAL_CONES_FACETS");
    Matrix<Rational> linearSpan = fan.give("LINEAR_SPAN_NORMALS");
    Matrix<Rational> linearSpanInCones = fan.give("MAXIMAL_CONES_LINEAR_SPAN_NORMALS");
    int fan_dim = fan.give("CMPLX_DIM");
    
     //First separate affine and directional rays
    Set<int> affineRays;
    Set<int> directionalRays;
    for(int r = 0; r < rays.rows(); r++) {
      if(rays.row(r)[0] == 0) {
	directionalRays = directionalRays + r;
      }
      else {
	affineRays = affineRays + r;
      }
    }  
    
    dbgtrace << "Computing bounding box..." << endl;
    
    //Compute facets of the bounding box
    Matrix<Rational> bbFacets(0,ambient_dim+1);
    
    //For each coordinate, contain minimum and maximum
    Vector<Rational> minCoord(ambient_dim);
    Vector<Rational> maxCoord(ambient_dim);
    //If bounding mode is relative, determine the maximal/minimal coordinates of the affine rays
    if(isRelative) {
      for(Entire<Set<int> >::iterator aff = entire(affineRays); !aff.at_end(); aff++) {
	for(int i = 0; i < ambient_dim; i++) {
	  Rational val = rays(*aff,i+1);
	  if(val > maxCoord[i]) maxCoord[i] = val;
	  if(val < minCoord[i]) minCoord[i] = val;
	}
      }
      //Now add the bbDistance to all values
      for(int i = 0; i < ambient_dim; i++) {
	maxCoord[i] += bbDistance;
	minCoord[i] -= bbDistance;
      }
      if(onlyBoundingBox) {
	Matrix<Rational> bb(0,ambient_dim);
	bb /= minCoord; 
	bb /= maxCoord;
	perl::ListReturn smallResult;
	  smallResult << bb;
	return smallResult;
      }
    }
    //otherwise take min and max from the given bounding box
    else {
      for(int i = 0; i < ambient_dim; i++) {
	maxCoord[i] = bBox(0,i) > bBox(1,i)? bBox(0,i) : bBox(1,i);
	minCoord[i] = bBox(0,i) < bBox(1,i)? bBox(0,i) : bBox(1,i);
      }
    }
    //Now make these coordinates into facets
    for(int i = 0; i < ambient_dim; i++) {
      Vector<Rational> facetVector = unit_vector<Rational>(ambient_dim,i);
      bbFacets /= (maxCoord[i] | -facetVector);
      bbFacets /= (-minCoord[i] | facetVector);
    }
    
    dbgtrace << "Done." << endl;
    
    perl::ListReturn result;
    
    //This will contain the cell centers with the weight labels
    perl::Object weightCenters("polytope::PointConfiguration");
    Matrix<Rational> centermatrix(0,ambient_dim);
    Vector<std::string> centerlabels;
    Array<Integer> weights;
    if(showWeights) {
      weights = fan.give("TROPICAL_WEIGHTS");
    }
    
    //Now compute all polyhedra to be rendered
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      dbgtrace << "Computing polytope of cone " << mc << endl;
      //Compute the facets ans equalities of the current cone and add the bbox facets
      Matrix<Rational> facets(0,ambient_dim+1);
      Matrix<Rational> linspan = linearSpan.minor(linearSpanInCones.row(mc),All);
	linspan = linspan;
      for(int fn = 0; fn < facetNormalsInCones.cols(); fn++) {
	if(facetNormalsInCones(mc,fn) == 1) {
	    facets /= facetNormals.row(fn);
	}
	if(facetNormalsInCones(mc,fn) == -1) {
	    facets /= (-facetNormals.row(fn));
	}
      }
      facets /= bbFacets;
      //facets = facets;
      
      dbgtrace << "Facets are " << facets << "Equalities are " << linspan << endl;
      
      //Compute the polytope vertices from that
      Matrix<Rational> polyRays = solver<Rational>().enumerate_vertices(zero_vector<Rational>()| facets, zero_vector<Rational>() | linspan).first;
      polyRays = polyRays.minor(All,~scalar2set(0));
      //We have to make sure that the polytope has
      //at least dim +1 vertices after cutting, otherwise its a point set or graph to the
      //visualization and all the Facet options don't work
      if(polyRays.rows() >= fan_dim+1) {
	perl::Object polytope("polytope::Polytope<Rational>");
	  polytope.take("VERTICES") << polyRays; //The polytope shouldn't have a lineality space
	result << polytope;
	
	//If weight labels should be displayed, compute the vertex barycenter of the polytope and
	// label it
	if(showWeights) {
	  Vector<Rational> barycenter = average(rows(polyRays));
	  centermatrix = centermatrix / barycenter;
	  std::ostringstream wlabel;
	  wlabel << "# " << mc << ": " << weights[mc];
	  centerlabels = centerlabels | wlabel.str();
	}
      }
//       }
//       catch(...) { //An error should only occur if the polytope is empty. Then just omit it
// 	dbgtrace << "Cone " << mc << " not in bounding box. Omitting." << endl;
// 	pm::cout << "Test" << endl;
//       }
      
    }
    
    if(showWeights) {
      weightCenters.take("POINTS") << centermatrix;
      weightCenters.take("LABELS") << centerlabels;
    }
    result << weightCenters;
    
    return result;
  }
  
// ------------------------- PERL WRAPPERS ---------------------------------------------------

Function4perl(&computeCodimensionOne,"computeCodimensionOne(WeightedComplex)");

//Function4perl(&computeCodimensionOneViaCones,"computeCodimensionOneViaCones(WeightedComplex;@)");

Function4perl(&computeLatticeNormals, "computeLatticeNormals(WeightedComplex)");

Function4perl(&computeLatticeNormalSum, "computeLatticeNormalSum(WeightedComplex)");

Function4perl(&computeIfBalanced, "computeIfBalanced(WeightedComplex)");

Function4perl(&computeComplexData, "computeComplexData(WeightedComplex)");

Function4perl(&computeFunctionVectors, "computeFunctionVectors(WeightedComplex)");

Function4perl(&computeVisualPolyhedra, "computeVisualPolyhedra(WeightedComplex, Rational, $)");

Function4perl(&computeBoundedVisual, "computeBoundedVisual(WeightedComplex, $, $, Rational,$, Matrix<Rational>)");

}}
