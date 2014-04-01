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
#include "polymake/atint/WeightedComplexRules.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  using polymake::polytope::cdd_interface::solver;
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  bool is_coneset_compatible(const Set<int> &cone, const IncidenceMatrix<> &local_restriction) {
    for(int i = 0; i < local_restriction.rows(); i++) { 
	Set<int> inter = cone * local_restriction[i];
	if(inter.size() == local_restriction.row(i).size()) {
	  return true;
	}
    }
    return false;
  }
  
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
    IncidenceMatrix<> local_restriction = fan.give("LOCAL_RESTRICTION");
    //Special case: fan is only origin
    if(maximalCones.rows() == 1) {
      if(maximalCones.row(0).size() == 0) {
	fan.take("CODIM_1_FACES") << IncidenceMatrix<>();
	fan.take("CODIM_1_IN_MAXIMAL_CONES") << IncidenceMatrix<>();
	return;
      }
    }
    
    CodimensionOneResult r = calculateCodimOneData(rays, maximalCones, uses_homog, linspace, local_restriction);
    
    fan.take("CODIM_1_FACES") << r.codimOneCones;
    fan.take("CODIM_1_IN_MAXIMAL_CONES") << r.codimOneInMaximal;    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  CodimensionOneResult calculateCodimOneData(const Matrix<Rational> &rays, const IncidenceMatrix<> &maximalCones, bool uses_homog, const Matrix<Rational> &linspace, const IncidenceMatrix<> &local_restriction) {
    //dbgtrace << "Computing all facets..." << endl;
    
    //First we construct the set of all facets 
    //Array<IncidenceMatrix<> > maximal_cone_incidence = fan.give("MAXIMAL_CONES_INCIDENCES");
    //Compute the rays-in-facets for each cone directly
    Vector<IncidenceMatrix<> > maximal_cone_incidence;
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      Set<int> mset = maximalCones.row(mc);
      //Extract inequalities
      //dbgtrace << "Computing facets for cone set " << mset << endl;
      Matrix<Rational> facets = solver<Rational>().enumerate_facets(
	    zero_vector<Rational>() | rays.minor(mset,All),
	    zero_vector<Rational>() | linspace).first.minor(All, ~scalar2set(0));;
      //dbgtrace << "Done. Checking rays..." << endl;
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
      //dbgtrace << "Done." << endl;
      maximal_cone_incidence |= IncidenceMatrix<>(facetIncidences);
    }
    
    //dbgtrace << "Check for doubles and useless facets..." << endl;
    
    //This will contain the set of indices defining the codim one faces
    Vector<Set<int> > facetArray;
    
    //This will define the codim-1-maximal-cone incidence matrix
    Vector<Set<int> > fIncones;
    
    for(int maxcone = 0; maxcone < maximal_cone_incidence.size(); maxcone++) {
      //This is the incidence matrix for the maximal cone indexed by maxcone
      IncidenceMatrix<> fcts = maximal_cone_incidence[maxcone];
      for(int facet = 0; facet < fcts.rows(); facet++) {
	Set<int> facetToCheck = fcts.row(facet);
	//If there is a local restriction, check if the facet is compatible
	if(local_restriction.rows() > 0) {
	  if(!is_coneset_compatible(facetToCheck, local_restriction)) continue;
	}
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
	  facetArray |= facetToCheck;
	  Set<int> singlecone;
	    singlecone = singlecone + maxcone;
	  fIncones |= singlecone;
	}
	else {
	fIncones[fcIndex] = fIncones[fcIndex] + maxcone;
	}
      }
    }
    
    CodimensionOneResult r;
    	r.codimOneCones = IncidenceMatrix<>(facetArray);
	r.codimOneInMaximal = IncidenceMatrix<>(fIncones);
    return r;    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Takes a polyhedral fan and computes the ray data of the corresponding polyhedral complex. Sets the corresponding properties CMPLX_RAYS, CMPLX_MAXIMAL_CONES, CMPLX_CODIM_1_FACES automatically.
    @param WeightedComplex fan A polyhedral fan, extended by atint to a tropical variety
  */
  void computeComplexData(perl::Object fan) {
    //Extract properties of fan
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> rays = fan.give("RAYS");
    
    if(!uses_homog) {
      fan.take("CMPLX_RAYS") << rays;
      fan.take("CMPLX_MAXIMAL_CONES") << fan.give("MAXIMAL_CONES");
      fan.take("CMPLX_CODIM_1_FACES") << fan.give("CODIM_1_FACES");
      fan.take("CMPLX_CONVERSION_VECTOR") << Vector<int>(sequence(0,rays.rows()));
      return;
    }
        
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    IncidenceMatrix<> facet_incidences = fan.give("CODIM_1_IN_MAXIMAL_CONES");
      facet_incidences = T(facet_incidences);
    IncidenceMatrix<> local_restriction = fan.give("LOCAL_RESTRICTION");
      
    //Result variables
    Matrix<Rational> cmplxrays(0,ambient_dim);
    Vector<Set<int> > maxcones(maximalCones.rows());
    Vector<Set<int> > codimone(codimOneCones.rows());
    Vector<int> conversion;
    
    //dbgtrace << "Dividing rays..." << endl;
    
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
	conversion |= r;
	newAffineIndices[r] = cmplxrays.rows()-1;
      }
    }
    
    //dbgtrace << "Affine rays: " << affineRays << ", directional rays: " << directionalRays << endl;
    
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
    
    //If there is a local restriction, we keep only compatible affine rays
    //for equivalence computation
    if(local_restriction.rows() > 0) {
//       affineRays *= accumulate(rows(local_restriction),operations::add());
      Set<int> compatible;
      for(Entire<Set<int> >::iterator af = entire(affineRays); !af.at_end(); af++) {
	Set<int> single; single += *af;
	if(is_coneset_compatible(single, local_restriction)) {
	    compatible += *af;
	}
      }
      affineRays = compatible;
    }
    
    //dbgtrace << "Added affine rays to cones" << endl;
    
    //Now we go through the directional rays and compute the connected component for each one
    for(Entire<Set<int> >::iterator r = entire(directionalRays); !r.at_end(); ++r) {
      
      //dbgtrace << "Computing components of ray " << *r << endl;
      
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
      
      //dbgtrace << "Computed set of cones containing r:" << rcones << endl;
      
      //For each such maximal cone, compute its component (if it hasnt been computed yet).
      for(Entire<Set<int> >::iterator mc = entire(rcones); !mc.at_end(); ++mc) {
	if(!inverseMap.exists(*mc)) {
	  //dbgtrace << "Creating new component" << endl;
	  //Create new component
	  Set<int> newset; newset = newset + *mc;
	  connectedComponents |= newset;
	  inverseMap[*mc] = connectedComponents.dim()-1;
	  
	  //Do a breadth-first search for all other cones in the component
	  std::list<int> queue;
	    queue.push_back(*mc);
	    //Semantics: Elements in that queue have been added but their neighbours might not
	    //dbgtrace << "Calculating component" << endl;
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
      
      //dbgtrace << "Connected components:\n" << connectedComponents << endl;
      
      //Now add r once for each connected component to the appropriate cones
      for(int cc = 0; cc < connectedComponents.dim(); cc++) {
	cmplxrays = cmplxrays / rays.row(*r);
	conversion |= (*r);
	int rowindex = cmplxrays.rows()-1;
	Set<int> ccset = connectedComponents[cc];
	//dbgtrace << "Inserting for component " << cc+1 << endl;
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
    
    //dbgtrace << "Done computing rays, inserting values..." << endl;
    
    //Insert values
    fan.take("CMPLX_RAYS") << cmplxrays;
    fan.take("CMPLX_MAXIMAL_CONES") << maxcones;
    fan.take("CMPLX_CODIM_1_FACES") << codimone;
    fan.take("CMPLX_CONVERSION_VECTOR") << conversion;
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  IncidenceMatrix<> computeAllCones(const Matrix<Rational> &rays, const IncidenceMatrix<> &maximalCones, bool uses_homog){
    Vector<Set<int> > result;
    IncidenceMatrix<> mcones = maximalCones;
    while(mcones.rows() > 0) {
	for(int r = 0; r < mcones.rows(); r++) {
	    result |= mcones.row(r);
	}
	//Check if only the empty cone is left. In that case, terminate
	if(mcones.rows() == 1 && mcones.row(0).size() == 0) break;
	CodimensionOneResult cr = calculateCodimOneData(rays, mcones, uses_homog, Matrix<Rational>(0,rays.cols()),IncidenceMatrix<>());
	mcones = cr.codimOneCones;  
    }
    return IncidenceMatrix<>(result);
  }
  
// ------------------------- PERL WRAPPERS ---------------------------------------------------

Function4perl(&computeCodimensionOne,"computeCodimensionOne(WeightedComplex)");

Function4perl(&computeComplexData, "computeComplexData(WeightedComplex)");

Function4perl(&computeAllCones, "computeAllCones(Matrix<Rational>, IncidenceMatrix, $)");


}}
