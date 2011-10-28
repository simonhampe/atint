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
#include "polymake/Array.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/linalg.h"
#include "polymake/atint/normalvector.h"
#include "polymake/atint/WeightedComplexRules.h"
#include "polymake/atint/refine.h"

namespace polymake { namespace atint {

    using polymake::polytope::cdd_interface::solver;

    using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see header -------------------------------------------------------------
    inline Rational functionValue(Matrix<Rational> functionMatrix, Vector<Rational> point, bool uses_min, bool uses_homog) {
      //Remove the first coordinate and add a 1 at the end of point for the constant coefficient
      if(uses_homog) point = point.slice(~scalar2set(0));
      point |= 1;
      Vector<Rational> listOfValues = functionMatrix * point;
      return uses_min? 	accumulate(Set<Rational>(listOfValues), operations::min())
			: accumulate(Set<Rational>(listOfValues), operations::max());
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //TODO: Replace
    //Documentation see header -------------------------------------------------------------
    perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan) {
      
      dbglog << "Starting refinement..." << endl;
      
      dbglog << "Extracting values" << endl;
      
      //Extract values
      dbglog << "Computing rays " << endl;
      Matrix<Rational> rays = fan.give("RAYS");
      dbglog << "Done." << endl;
      Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
      int ambient_dim = rays.cols() < linspace.cols() ? linspace.cols() : rays.cols();
      IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
      bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
      int dimension = fan.give("CMPLX_DIM");
	//If fan is zero-dimensional, no refinement is necessary anyway, we just return fan	
	if(dimension == 0) return fan;
	//We actually want to compute with the homogeneous dimension
	if(uses_homog) dimension++;
      
      Array<Integer> weights;
      bool weightsExist = false;
      if(fan.exists("TROPICAL_WEIGHTS")) {
	weights = fan.give("TROPICAL_WEIGHTS");
	weightsExist = true;
      }
      
      
      Matrix<Rational> compRays = completeFan.give("RAYS");
      IncidenceMatrix<> compMaximalCones = completeFan.give("MAXIMAL_CONES");
      Matrix<Rational> compLinealitySpace = completeFan.give("LINEALITY_SPACE");

      dbglog << "Done. Computing intersection" << endl;
      
      //Result variables
      Matrix<Rational> newRays(0,ambient_dim);
      Vector<Set<int> > newCones;
      Set<Set<int> > newConeSet;//This is kept to check for doubles. It has the same content as maximalCones
      Matrix<Rational> newLineality(0,ambient_dim);
	bool computedLineality = false;
      Vector<Integer> newWeights;
            
      //If the complete fan is only a lineality space, it is the whole space, so we return fan
      if(compMaximalCones.rows() == 0) {
	return fan;
      }
      
      //If the fan, however, is only a lineality space, we have to take care that it is considered
      // as a maximal cone of fan (note that at this point, the lineality dimension must then be > 0,
      // since the fan dimension is > 0
      bool fanHasOnlyLineality = maximalCones.rows() == 0;
      int fanUpperBound = fanHasOnlyLineality? 1 : maximalCones.rows();
      
      dbgtrace << "Intersecting cones" << endl;
      
      //Now intersect all pairs of maximal cones.
      for(int fanIndex = 0; fanIndex < fanUpperBound; fanIndex++) {
	//Compute inequalities and equations for the fan cone
	std::pair<Matrix<Rational>, Matrix<Rational> > fanEqs =
	    solver<Rational>().enumerate_facets(fanHasOnlyLineality? Matrix<Rational>(0,ambient_dim) :
				zero_vector<Rational>() | rays.minor(maximalCones.row(fanIndex),All),
				zero_vector<Rational>() | linspace,true);
	for(int cIndex = 0; cIndex < compMaximalCones.rows(); cIndex++) {
	    std::pair<Matrix<Rational>, Matrix<Rational> > compEqs =
		solver<Rational>().enumerate_facets(
			      zero_vector<Rational>() | compRays.minor(compMaximalCones.row(cIndex),All),
			      zero_vector<Rational>() | compLinealitySpace,true);
	    
	    Matrix<Rational> interinequalities = fanEqs.first / compEqs.first;
	    Matrix<Rational> interequalities = fanEqs.second / compEqs.second;
	    std::pair<Matrix<Rational>, Matrix<Rational> > s = 	
		  solver<Rational>().enumerate_vertices(interinequalities, fanEqs.second / compEqs.second,true);
	    s.first = s.first.minor(All, range(1,s.first.cols()-1));
	    s.second = s.second.minor(All,range(1,s.second.cols()-1));
	    
	    dbgtrace << "Computed rays of intersection" << endl;
		  
	    int coneDimension = rank(s.first) + rank(s.second);
	    
	    //Only consider intersections that have the correct dimension
	    if(coneDimension == dimension) {
		dbgtrace << "Rays: " << s.first << endl;
		dbgtrace << "Lineality: " << s.second << endl;
		//The first time we compute this, the lineality space of the intersection is copied as the lineality space of the resulting fan
		if(!computedLineality) {
		  newLineality = s.second;
		  //Make sure the lineality space has the right dimension
		  if(newLineality.cols() <= 0) {
		    newLineality = Matrix<Rational>(0,ambient_dim);
		  }
		  computedLineality = true;
		  dbgtrace << "Setting lineality space to " << newLineality << endl;
		  
		}
		
		//Now go through the rays of the intersection to find those already added
		Matrix<Rational> interRays = s.first;//inter.give("RAYS");
		
		dbgtrace << "Checking rays" << endl;
		
		Set<int> remainingRows;
		Set<int> newMaximalCone;
		for(int ray = 0; ray < interRays.rows(); ray++) {
		  int raysIndex = -1;
		  //Compare it to all existing rays
		  for(int testRay = 0; testRay < newRays.rows(); testRay++) {
		    if(interRays.row(ray) == newRays.row(testRay)) {
		      raysIndex = testRay;
		      break;
		    }
		  }
		  //If we found the ray, add the index to the maximal cone, otherwise
		  //keep track of the row to create new rays later
		  if(raysIndex >= 0) {
		    newMaximalCone += raysIndex;
		  }
		  else {
		    remainingRows += ray;
		  }
		}//END go through all intersection rays
		
		//Now add remaining rays and the corresponding cone (if it doesn't exist yet)
		bool addCone = false;
		if(remainingRows.size() > 0) {
		  addCone = true;
		}
		else {
		  addCone = !(newConeSet.contains(newMaximalCone));
		}
		
		if(addCone) {
		    dbgtrace << "Adding cone" << endl;		  
		    newRays = newRays / interRays.minor(remainingRows,All);
		    for(int i = newRays.rows() - remainingRows.size(); i < newRays.rows(); i++) {
		      newMaximalCone += i;
		    }
		    newCones = newCones | newMaximalCone;
		    newConeSet += newMaximalCone;
		    if(weightsExist) {
			newWeights = newWeights | weights[fanIndex];
		    }
		}
	    }//End check intersection cone
	    
	   
	}//END iterate over complete fan cones
      }//END iterate over fan cones
      
      //Return result
      perl::Object result("WeightedComplex");
	result.take("RAYS") << newRays;
	result.take("MAXIMAL_CONES") << newCones;
	result.take("USES_HOMOGENEOUS_C") << uses_homog;
	result.take("LINEALITY_SPACE") << newLineality;
	  int cmplx_dim = uses_homog? dimension-1 : dimension;
	result.take("CMPLX_DIM") << cmplx_dim;
	if(weightsExist) result.take("TROPICAL_WEIGHTS") << newWeights;
      return result;
      
      
    }

    ///////////////////////////////////////////////////////////////////////////////////////

    //TODO: Replace
    //Documentation see header -------------------------------------------------------------
    perl::Object divisorByValueVector(perl::Object fan, Vector<Rational> values) {
      //Extract properties
      bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
      Matrix<Rational> rays = uses_homog? fan.give("CMPLX_RAYS") : fan.give("RAYS");
      int lineality_dim = fan.give("LINEALITY_DIM");
      Matrix<Rational> lineality_space = fan.give("LINEALITY_SPACE");
      Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
      IncidenceMatrix<> codimOneCones = uses_homog? fan.give("CMPLX_CODIM_1_FACES") : fan.give("CODIM_1_FACES");
      IncidenceMatrix<> maximalCones = uses_homog? fan.give("CMPLX_MAXIMAL_CONES") : fan.give("MAXIMAL_CONES");
      IncidenceMatrix<> coneIncidences = fan.give("CODIM_1_IN_MAXIMAL_CONES");
      Array<Integer> tropicalWeights = fan.give("TROPICAL_WEIGHTS");
      Map<int, Map<int, Vector<Rational> > > lnFunctionVector = fan.give("LATTICE_NORMAL_FCT_VECTOR");
      Matrix<Rational> lsumFunctionVector = fan.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
      
      //If the fan only consists of the lineality space, the divisor is empty
      int noOfRays = rays.rows();
      if(noOfRays == 0) {
	return CallPolymakeFunction("zero_cycle");
      }
      
      dbgtrace << "Making sure values have right length" << endl;
      
      //Append zero values, if necessary
      while(values.dim() < noOfRays + lineality_dim) {
	values = values | 0;
      }
      
      //Slice away superfluous values
      if(values.dim() > noOfRays + lineality_dim) {
	values = values.slice(0,noOfRays + lineality_dim);
      }
      
      dbgtrace << "Done." << endl;
      
      //Result variables
      Vector<Set<int> > newcones;
      Vector<Integer> newweights;
      
      dbgtrace << "Computing facet weights" << endl;
      
      //Go through each facet and compute its weight. Only add it, if the weight is nonzero
      Matrix<Rational> newrays(0,rays.cols()); //Will contain the rays that remain after removing weight-0-facets
      Vector<bool> hasBeenAdded(rays.rows()); //True, if a facet containing that ray has been added with > 0 weight
      Map<int,int> oldToNew; //Maps old ray indices to ray indices in newrays
      
      for(int co = 0; co < codimOneCones.rows(); co++) {
	dbgtrace << "Codim 1 face " << co << endl;
	Integer coweight(0);
	Set<int> adjacentCones = coneIncidences.row(co);
	for(Entire<Set<int> >::iterator mc = entire(adjacentCones); !mc.at_end(); ++mc) {
	  dbgtrace << "Maximal cone " << *mc << endl;
	  coweight = coweight + tropicalWeights[*mc] * (lnFunctionVector[co])[*mc] * values;
	}
	//Now substract the value of the lattice normal sum
	dbgtrace << "Substracting sum" << endl;
	coweight = coweight - lsumFunctionVector.row(co) * values;
	
	if(coweight != 0) {
	    newweights = newweights | coweight;
	    //Go through the rays, check which one has been added and take its index, otherwise create a new row
	    Set<int> oldConeRays = codimOneCones.row(co);
	    Set<int> newConeRays;
	    for(Entire<Set<int> >::iterator r = entire(oldConeRays); !r.at_end(); r++) {
	      if(hasBeenAdded[*r]) {
		newConeRays += oldToNew[*r];
	      }
	      else {
		hasBeenAdded[*r] = true;
		newrays /= rays.row(*r);
		oldToNew[*r] = newrays.rows()-1;
		newConeRays += (newrays.rows()-1);
	      }
	    }
	    newcones = newcones | newConeRays;
	}
      }
      
      //If the lineality space is the codim 1 cone and its weight is 0, the divisor is still empty
      if(newweights.dim() == 0) {
	return CallPolymakeFunction("zero_cycle");
      }
      
      dbgtrace << "Done\n" << endl;
      
      //Return result
      perl::Object result("WeightedComplex");
	result.take("RAYS") << newrays;
	result.take("USES_HOMOGENEOUS_C") << uses_homog;
	result.take("MAXIMAL_CONES") << newcones;
	result.take("TROPICAL_WEIGHTS") << newweights;
	result.take("LINEALITY_SPACE") << lineality_space; 
      return result;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////

    //Documentation see header -------------------------------------------------------------
    perl::Object divisorByValueMatrix(perl::Object complex, Matrix<Rational> values) {
      //This value carries all the intermediate results.
      perl::Object result = complex;      
      
      //Now we extract the values that we will later recompute by hand or that don't change at all
      
      bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
      Matrix<Rational> rays = complex.give("RAYS");
      Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
      Matrix<Rational> lineality_space = complex.give("LINEALITY_SPACE");
      int lineality_dim = complex.give("LINEALITY_DIM");
      
      Matrix<Rational> lineality_values = values.minor(All,~(sequence(0,values.cols() - lineality_dim)));
      
      //Prepare the additional variables that will be used in all but the first iteration to recompute the
      //values vector
      
      //Contains at position i the row index of ray i in the ray matrix of the iteration before
      Vector<int> newRaysToOldRays; 
      
      //Contains the CMPLX_MAXIMAL_CONES of the iteration before 
      IncidenceMatrix<> cmplx_oldcones;
      
      //Tells which new maximal cone is contained in which old maximal cone (this is essentially the
      //CODIM_1_IN_MAXIMAL_CONES without the rows for cones of weight 0)
      IncidenceMatrix<> newConesInOld;
      
      //For each cmplx_ray in the LAST iteration, this tells which should be the appropriate
      //column index in values for function value computation
      Vector<int> cmplx_origins = sequence(0,values.cols() - lineality_dim);
      
      //Contains the conversion vector for the last iteration (this one we recompute during
      //value recomputation)
      Vector<int> old_conversion;
      
      //Now we iterate through the matrix rows 
      for(int r = 0; r < values.rows(); r++) {
	//First we recompute values that we can't/won't compute by hand

	IncidenceMatrix<> codimOneCones = result.give("CODIM_1_FACES");
	  if(codimOneCones.rows() == 0) return CallPolymakeFunction("zero_cycle");
	IncidenceMatrix<> coneIncidences = result.give("CODIM_1_IN_MAXIMAL_CONES");
      
	Map<int, Map<int, Vector<Rational> > > lnFunctionVector = result.give("LATTICE_NORMAL_FCT_VECTOR");
	Matrix<Rational> lsumFunctionVector = result.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
	Map<int, Map<int, Vector<Integer> > > latticeNormals = result.give("LATTICE_NORMALS");
	
	//Now we compute the correct value vector:
	Vector<Rational> currentValues;
	if(r == 0 || !uses_homog) {
	  currentValues = values.row(r);
	}
	else {
	  Matrix<Rational> cmplx_rays = result.give("CMPLX_RAYS");
	  Vector<int> conversion_vector = result.give("CMPLX_CONVERSION_VECTOR");
	  //Compute the maximal cones containing each cmplx_ray
	  IncidenceMatrix<> cmplx_cones_t = result.give("CMPLX_MAXIMAL_CONES");
	    cmplx_cones_t = T(cmplx_cones_t);
	    
	  Vector<int> newcmplx_origins;  
	  for(int cr = 0; cr < cmplx_rays.rows(); cr++) {
	    //Find the corresponding cmplx_ray in the last iteration
	    int mc = *(cmplx_cones_t.row(cr).begin()); //A cone containing the ray
	    int oc = *(newConesInOld.row(mc).begin()); //An old cone containing mc
	    //Now find the cmplx_ray of the old cone, such that 
	    Set<int> ocrays = cmplx_oldcones.row(oc);
	    for(Entire<Set<int> >::iterator ocr = entire(ocrays); !ocr.at_end(); ocr++) {
	      //If the old ray (in non-complex counting in the old iteration) is the same as 
	      //the new ray (in non-complex counding in the new iteration, we can
	      //copy its function column index
	      if(old_conversion[*ocr] == newRaysToOldRays[conversion_vector[cr]]) {
		currentValues |= values(r,cmplx_origins[*ocr]);
		newcmplx_origins |= cmplx_origins[*ocr];
		break;
	      }
	    }
	  }
	  cmplx_origins = newcmplx_origins;
	  //Finally append lineality values
	  currentValues |= lineality_values.row(r);
	}
	
	//Then we compute the divisor
	Vector<Integer> newweights; //Contains the new weights
	Set<int> usedCones; //Contains the codim 1 cones with weight != 0
	Set<int> usedRays; //Contains the rays in used cones
	//Go through each facet and compute its weight. 
	for(int co = 0; co < codimOneCones.rows(); co++) {
// 	  dbgtrace << "Codim 1 face " << co << endl;
	  Integer coweight(0);
	  Set<int> adjacentCones = coneIncidences.row(co);
	  for(Entire<Set<int> >::iterator mc = entire(adjacentCones); !mc.at_end(); ++mc) {
// 	    dbgtrace << "Maximal cone " << *mc << endl;
	    coweight = coweight + weights[*mc] * (lnFunctionVector[co])[*mc] * currentValues;
	  }
	  //Now substract the value of the lattice normal sum
// 	  dbgtrace << "Substracting sum" << endl;
	  coweight = coweight - lsumFunctionVector.row(co) * currentValues;
	  if(coweight != 0) {
	    newweights = newweights | coweight;	  
	    usedCones += co;
	    usedRays += codimOneCones.row(co);
	  }
	}//END iterate co-1-cones
	
	//Compute the new-to-old maps used for recomputing the value vector in the next iteration
	if(r != values.rows()-1) {
	  
	  newConesInOld = coneIncidences.minor(usedCones,All);	  
	  //Need this to avoid compiler errors about ambiguos overloads of .give
	  IncidenceMatrix<> oc = result.give("CMPLX_MAXIMAL_CONES"); 
	    cmplx_oldcones = oc;
	  Vector<int> ocv = result.give("CMPLX_CONVERSION_VECTOR");
	    old_conversion = ocv;
	  newRaysToOldRays = Vector<int>();
	  for(Entire<Set<int> >::iterator orays = entire(usedRays); !orays.at_end(); orays++) {
	    newRaysToOldRays |= (*orays);
	  }
	}
	
	//Now recompute the rays and maximal cones for re-initialization of the result
	rays = rays.minor(usedRays,All);
	weights = newweights;
	result = perl::Object("WeightedComplex");
	  result.take("USES_HOMOGENEOUS_C") << uses_homog;
	  result.take("RAYS") << rays;
	  result.take("MAXIMAL_CONES") << codimOneCones.minor(usedCones,usedRays);
	  result.take("TROPICAL_WEIGHTS") << weights;
	  result.take("LINEALITY_SPACE") << lineality_space;
	
	
      } //END iterate function rows
      
      return result;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //TODO: Replace   
    //Documentation see header -------------------------------------------------------------
    perl::Object divisorByPLF(perl::Object fan, perl::Object function) {
      dbglog << "Preparing computations" << endl;
      
      //Homogenize the fan and refine it
      bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
      if(!uses_homog) fan = fan.CallPolymakeMethod("homogenize");
      perl::Object linearityDomains = function.give("DOMAIN");
      dbglog << "Refining fan" << endl;
      fan = intersect_complete_fan(fan, linearityDomains);
      
      dbglog << "Extracting values for divisor computation" << endl;
      
      //Extract values
      Matrix<Rational> rays = fan.give("CMPLX_RAYS");	
      Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
      Array<Set<int> > maximal = fan.give("CMPLX_MAXIMAL_CONES");
      Matrix<Rational> fmatrix = function.give("FUNCTION_MATRIX");
      bool uses_min = function.give("USES_MIN");
      
      dbglog << "Extracted values" << endl;
      
      //Now compute function values
      Vector<Rational> values;
      
      //Separate rays
      Set<int> affine, directional;
      for(int r = 0; r < rays.rows(); r++) {
	if(rays(r,0) == 0) {
	  directional += r;
	}
	else {
	  affine += r;
	}
      }
      
      dbglog << "Computing values" << endl;
	
      for(int index = 0; index < rays.rows(); index++) {
	//If it is an affine ray, simply compute the function value at that point
	if(rays(index,0) == 1) {
	    values |= functionValue(fmatrix, rays.row(index),uses_min,true);
	}
	//Otherwise find an affine ray that shares a maximal cone with this ray (if there is none, set the value to 0, it doesn't matter anyway)
	else {
	    for(int mc = 0; mc < maximal.size(); mc++) {
	      if(maximal[mc].contains(index)) {
		Set<int> sharedRays = maximal[mc] * affine;
		if(sharedRays.size() > 0) {
		    //If we have found such an affine ray, compute the value of affine + direction and substract the value of the affine ray and stop searching
		    int sray = *(sharedRays.begin());
		    values |= (functionValue(fmatrix, rays.row(sray) + rays.row(index),uses_min,true) - 
			      functionValue(fmatrix, rays.row(sray),uses_min,true));
		    break;
		}
	      }
	    }
	}
      }
      
      //Finally we add the function values on the lineality space
      for(int index = 0; index < linspace.rows(); index++) {
	values |= functionValue(fmatrix, linspace.row(index), uses_min,true);
      }
      
      return divisorByValueVector(fan, values);
      
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see header
    perl::Object divisor_minmax(perl::Object complex, perl::Object function, int k) {
      dbglog << "Preparing computations" << endl;
      
      //Homogenize the fan and refine it
      bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
      if(!uses_homog) complex = complex.CallPolymakeMethod("homogenize");
      perl::Object linearityDomains = function.give("DOMAIN");
      
      dbglog << "Refining fan" << endl;
      RefinementResult r = refinement(complex,linearityDomains,false,false,true,true);
      
      //Extract values
      Matrix<Rational> rays = r.complex.give("CMPLX_RAYS");	
      Matrix<Rational> linspace = r.complex.give("LINEALITY_SPACE");
      Array<Set<int> > maximal = r.complex.give("CMPLX_MAXIMAL_CONES");
      Vector<int> assocRep = r.associatedRep;
      
      Matrix<Rational> fmatrix = function.give("FUNCTION_MATRIX");
      bool uses_min = function.give("USES_MIN");
      
      dbglog << "Extracted values" << endl;
      
      //Now compute function values
      Vector<Rational> values;    
      for(int r = 0; r < rays.rows(); r++) {
	//If it is an affine ray, simply compute the function value at that point
	if(rays(r,0) == 1) {
	    values |= functionValue(fmatrix, rays.row(r),uses_min,true);
	}
	//Otherwise take the function difference between (x+this ray) and x for an associated vertex x
	else {
	    values |=  (functionValue(fmatrix, rays.row(assocRep[r]) + rays.row(r),uses_min,true) - 
			      functionValue(fmatrix, rays.row(assocRep[r]),uses_min,true));
	}
      }
      //Finally we add the function values on the lineality space
      for(int index = 0; index < linspace.rows(); index++) {
	values |= functionValue(fmatrix, linspace.row(index), uses_min,true);
      }
      
      //Glue together to a value matrix
      Matrix<Rational> vmatrix(0,values.dim());
      for(int l = 1; l <= k; l++) {
	vmatrix /= values;
      }
      
      return divisorByValueMatrix(r.complex,vmatrix);
      
    }
    
    //Documentation see perl wrapper
    perl::Object divisor_nr(perl::Object complex, perl::Object function, int k) {
      //Extract function
      Vector<Rational> rvalues = function.give("RAY_VALUES");
      Vector<Rational> lvalues = function.give("LIN_VALUES");
      Vector<Rational> values = rvalues | lvalues;
      int power = function.give("POWER");

      //Create function matrix
      Matrix<Rational> fmatrix(0,values.dim());
      for(int l = 1; l <= (k==-1? power : k); l++) {
	fmatrix /= values;
      }
      
      return divisorByValueMatrix(complex,fmatrix);
    }
    
// ------------------------- PERL WRAPPERS ---------------------------------------------------
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete "
		      "# (i.e. its support is the whole ambient space) to compute the intersection correctly."
		      "# @param WeightedComplex fan An arbitrary weighted polyhedral fan"
		      "# @param fan::PolyhedralFan completeFan A complete polyhedral fan"
		      "# @return WeightedComplex The intersection of both fans (whose support is equal to the support of fan). The "
		      "# resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, "
		      "# the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned." ,
		      &intersect_complete_fan,"intersect_complete_fan(WeightedComplex, fan::PolyhedralFan)");
    
    Function4perl(&divisorByValueVector,"divisorByValueVector(WeightedComplex, Vector<Rational>)");  
    
    Function4perl(&divisorByValueMatrix, "divisorByValueMatrix(WeightedComplex, Matrix<Rational>)");
    
    Function4perl(&divisor_minmax,"divisor_minmax(WeightedComplex,MinMaxFunction;$=1)");
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Computes the divisor of a MinMaxFunction on a given tropical variety. The result will be "
		      "# in homogeneous coordinates, whether the tropical variety uses them or not. The function "
		      "# should be given on the affine coordinates of the variety, NOT the homogeneous ones."
		      "# @param WeightedComplex fan A tropical variety, on which the divisor is computed"
		      "# @param MinMaxFunction A function whose DOMAIN should be equal to the affine coordinate "
		      "# space of the variety, i.e. AMBIENT_DIM-1, if the variety uses homogeneous coordinates, " 
		      "# AMBIENT_DIM otherwise."
		      "# @return The corresponding divisor as a tropical variety in homogeneous coordinates.",
		      &divisorByPLF, "divisorByPLF(WeightedComplex,MinMaxFunction)");
    
   UserFunction4perl("# @category Tropical geometry"
		     "# Works exactly as divisor(WeightedComplex, RationalFunction;Int). Should be called ONLY,"
		     "# when the function f is defined on a DOMAIN equal to X (in the sense that all properties"
		     "# like RAYS, MAXIMAL_CONES, etc. agree. Being equal as varieties is not sufficient). In this"
		     "# case this function will in general be faster.",
		     &divisor_nr,"divisor_nr(WeightedComplex,RationalFunction;$=-1)");
   
}
}