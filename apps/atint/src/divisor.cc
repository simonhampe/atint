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

    using namespace atintlog::donotlog;
// // 	    using namespace atintlog::dolog;
//     using namespace atintlog::dotrace;
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see header -------------------------------------------------------------
    Rational functionValue(Matrix<Rational> functionMatrix, Vector<Rational> point, bool uses_min, bool uses_homog) {
      //Remove the first coordinate and add a 1 at the end of point for the constant coefficient
      if(uses_homog) point = point.slice(~scalar2set(0));
      point |= 1;
      Vector<Rational> listOfValues = functionMatrix * point;
      return uses_min? 	accumulate(Set<Rational>(listOfValues), operations::min())
			: accumulate(Set<Rational>(listOfValues), operations::max());
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see header -------------------------------------------------------------
    perl::Object intersect_container(perl::Object fan, perl::Object completeFan, bool forceLatticeComputation) {
	 RefinementResult r = refinement(fan,completeFan,false,false,false,true,forceLatticeComputation);
	 return r.complex;
    }

    ///////////////////////////////////////////////////////////////////////////////////////

    //Documentation see header -------------------------------------------------------------
    perl::Object divisorByValueMatrix(perl::Object complex, Matrix<Rational> values, bool uses_min) {
      //This value carries all the intermediate results.
      perl::Object result = complex;      
      
      //Now we extract the values that we will later recompute by hand or that don't change at all
      
      bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
      Matrix<Rational> rays = complex.give("RAYS");
      Matrix<Rational> crays = complex.give("CMPLX_RAYS");
      Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
      Matrix<Rational> lineality_space = complex.give("LINEALITY_SPACE");
      int lineality_dim = complex.give("LINEALITY_DIM");
      IncidenceMatrix<> local_restriction = complex.give("LOCAL_RESTRICTION");
      
      Matrix<Integer> lattice_generators = complex.give("LATTICE_GENERATORS");
      IncidenceMatrix<> lattice_bases = complex.give("LATTICE_BASES");
      
  //dbgtrace << "Rays: " << crays << endl;
      //dbgtrace << "Values: " << values << endl;
      
      //Do a compatibility check on the value matrix to avoid segfaults in the case of faulty input
      if(values.cols() != crays.rows() + lineality_space.rows()) {
	  throw std::runtime_error("Value matrix is not compatible with variety. Aborting computation");
      }
      
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
      Vector<int> cmplx_origins (sequence(0,values.cols() - lineality_dim));
      
      //Contains the conversion vector for the last iteration (this one we recompute during
      //value recomputation)
      Vector<int> old_conversion;
      
      //Only uses in the fan case: For all iterations but the first it contains the set of rays of
      //the last iteration that remained after computing the divisor
      Set<int> remainingFanRays;
      
      //When computing the codim-one-weights, this contains the correct function value vector for the current iteration
      //When computing the new function vector for the current iteration, this means it contains the function
      //values of the old iteration
      Vector<Rational> currentValues;
      
      //Now we iterate through the matrix rows 
      for(int r = 0; r < values.rows(); r++) {
	//dbgtrace << "Computing on row " << r << endl;
	//First we recompute values that we can't/won't compute by hand

	IncidenceMatrix<> codimOneCones = result.give("CODIM_1_FACES");
	  if(codimOneCones.rows() == 0) return CallPolymakeFunction("zero_cycle");
	IncidenceMatrix<> coneIncidences = result.give("CODIM_1_IN_MAXIMAL_CONES");
      
	Map<int, Map<int, Vector<Integer> > > latticeNormals = result.give("LATTICE_NORMALS");
	Map<int, Map<int, Vector<Rational> > > lnFunctionVector = result.give("LATTICE_NORMAL_FCT_VECTOR");
	Matrix<Rational> lsumFunctionVector = result.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
	Vector<bool> balancedFaces = result.give("BALANCED_FACES");
	
	//dbgtrace << "Balanced faces: " << balancedFaces << endl;
	
	//Recompute the lattice bases
	Vector<Set<int> > new_lattice_bases;
	for(int co = 0; co < codimOneCones.rows(); co++) {
	  new_lattice_bases |= lattice_bases.row(*(coneIncidences.row(co).begin()));
	  //dbgtrace << "Co " << co << ": " << new_lattice_bases[new_lattice_bases.dim()-1] << " from " << *(coneIncidences.row(co).begin()) << 	  endl;
	}
	lattice_bases = new_lattice_bases;
	//dbgtrace << "lb: " << lattice_bases << endl;
	
	//Now we compute the correct value vector:
	
	if(r == 0 || !uses_homog) {
	  if(r == 0) {
	    currentValues = values.row(r);
	  }
	  //We treat the fan case specially, since there aren't any new rays, just rays disappearing
	  else {
	    //Add values for remaining rays
	    Set<int> value_set;
	    for(Entire<Set<int> >::iterator ry = entire(remainingFanRays); !ry.at_end(); ry++) {
	      value_set += cmplx_origins[*ry];
	    }
	    cmplx_origins = cmplx_origins.slice(remainingFanRays);
	    currentValues = values.row(r).slice(value_set);
	    currentValues |= lineality_values.row(r);
	  }
	}
	else {
	  currentValues = Vector<Rational>();
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
	    //its corresponding ray is equal to the corresponding ray of the new ray
	    Set<int> ocrays = cmplx_oldcones.row(oc);
	    for(Entire<Set<int> >::iterator ocr = entire(ocrays); !ocr.at_end(); ocr++) {
	      //If the old ray (in non-complex counting in the old iteration) is the same as 
	      //the new ray (in non-complex counting) in the new iteration, we can
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
	//dbgtrace << "Value vector is: " << currentValues << endl;
	
	//Then we compute the divisor
	Vector<Integer> newweights; //Contains the new weights
	Set<int> usedCones; //Contains the codim 1 cones with weight != 0
	Set<int> usedRays; //Contains the rays in used cones
	//Go through each facet and compute its weight. 
	for(int co = 0; co < codimOneCones.rows(); co++) {
	  if(balancedFaces[co]) { //Only compute values at balanced codim-1-cones
  	  //dbgtrace << "Codim 1 face " << co << endl;
	    Rational coweight(0); //Have to take rational since intermediate values may be rational
	    Set<int> adjacentCones = coneIncidences.row(co);
	    for(Entire<Set<int> >::iterator mc = entire(adjacentCones); !mc.at_end(); ++mc) {
  	    //dbgtrace << "Maximal cone " << *mc << endl;
	      coweight = coweight + weights[*mc] * (lnFunctionVector[co])[*mc] * currentValues;
 	      //dbgtrace <<(lnFunctionVector[co])[*mc] * currentValues << endl; 
	    }
	    //Now substract the value of the lattice normal sum
   	  //dbgtrace << "Substracting sum" << endl;
 	    //dbgtrace << lsumFunctionVector.row(co) * currentValues << endl;
	    coweight = coweight - lsumFunctionVector.row(co) * currentValues;
	    if(coweight != 0) {
	      //Invert weight sign for min people.
	      if(uses_min) coweight = - coweight;
	      newweights = newweights | Integer(coweight);	  
	      usedCones += co;
	      usedRays += codimOneCones.row(co);
	    }
	  }
	}//END iterate co-1-cones
	
	//dbgtrace << "Remaining " << usedCones.size() << " of " << codimOneCones.rows() << " codim one cones" << endl;
	//dbgtrace << "Removing " << sequence(0,codimOneCones.rows()) - usedCones << endl;
	
	//dbgtrace << "Computed codim one weights" << endl;
	//dbgtrace << "Weights are " << newweights << endl;
	
	//Compute the new-to-old maps used for recomputing the value vector in the next iteration
	if(r != values.rows()-1) {
	  remainingFanRays = usedRays;
	  newConesInOld = coneIncidences.minor(usedCones,All);	  
	  result.give("CMPLX_MAXIMAL_CONES") >> cmplx_oldcones; 
	  result.give("CMPLX_CONVERSION_VECTOR") >> old_conversion;
	  newRaysToOldRays = Vector<int>();
	  for(Entire<Set<int> >::iterator orays = entire(usedRays); !orays.at_end(); orays++) {
	    newRaysToOldRays |= (*orays);
	  }
	}
	//dbgtrace << "newConesInOld: " << newConesInOld << endl;
	//dbgtrace << "newRaysToOldRays:" << newRaysToOldRays << endl;
	
	//Now recompute the rays and maximal cones for re-initialization of the result
	rays = rays.minor(usedRays,All);
	weights = newweights;
	IncidenceMatrix<> newMaximal = codimOneCones.minor(usedCones,usedRays);
	//Recompute local restriction cones
	if(local_restriction.rows() > 0) {
	  //We need to adapt rays indices and remove old maximal local cones
	  // and codimension one cones that have weight 0
	  //Also we remove all local cones that lose rays
	  //dbgtrace << "Local restriction before: " << local_restriction << endl;
	  IncidenceMatrix<> maxCones = result.give("MAXIMAL_CONES");
	  Set<int> removableCones;
	  Set<int> weightzerocones = sequence(0,codimOneCones.rows()) - usedCones;
	  Set<int> codimToReplace; //Indices of used codim one cones that are local
	  for(int lc = 0; lc < local_restriction.rows(); lc++) {
	    //If the local cone loses any rays, remove it
	    if((local_restriction.row(lc) * usedRays).size() < local_restriction.row(lc).size()) {
		removableCones += lc;
		continue;
	    }
	    bool found_cone = false;
	    for(int mc = 0; mc < maxCones.rows(); mc++) {
	      if((local_restriction.row(lc) * maxCones.row(mc)).size() == maxCones.row(mc).size()) {
		removableCones += lc;
		found_cone = true; break;
	      }
	    }
	    for(Entire<Set<int> >::iterator cz = entire(weightzerocones); !cz.at_end() && !found_cone; cz++) {
	      if((local_restriction.row(lc) * codimOneCones.row(*cz)).size() == codimOneCones.row(*cz).size()) {
		removableCones += lc;
		break;
	      }
	    }
	  }
	  
	  //Remove cones
	  local_restriction = local_restriction.minor(~removableCones, usedRays);
	  
// 	  //Replace all local codim one cones
// 	  Vector<Set<int> > interior_cones;
// 	  for(Entire<Set<int> >::iterator rc = entire(codimToReplace); !rc.at_end(); rc++) {
// 	    //Compute all intersections with remaining local cones
// 	    for(int lr = 0; lr < local_restriction.rows(); lr++) {
// 	      Set<int> isection = codimOneCones.row(*rc) * local_restriction.row(lr);
// 	      if(isection.size() > 0) {
// 		//Check for doubles and sets that contain this set (or vice versa)
// 		bool iscontained = false;
// 		Set<int> containedInSet;
// 		
// 	      }
// 	    }
// 	    
// 	  }
 	  
	  //dbgtrace << "Adapted local cones: " << local_restriction << endl;
	}//END adapt local restriction	
	
	result = perl::Object("WeightedComplex");
	  result.take("USES_HOMOGENEOUS_C") << uses_homog;
	  result.take("RAYS") << rays;
	  result.take("MAXIMAL_CONES") << newMaximal;
	  result.take("TROPICAL_WEIGHTS") << weights;
	  result.take("LINEALITY_SPACE") << lineality_space;
	  result.take("LOCAL_RESTRICTION") << local_restriction;
	  result.take("LATTICE_GENERATORS") << lattice_generators;
	    lattice_bases = lattice_bases.minor(usedCones,All);
	  result.take("LATTICE_BASES") << lattice_bases;//(lattice_bases.minor(usedCones,All));
	
      } //END iterate function rows

      //dbgtrace << "Done. Returning divisor" << endl;
      
      return result;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////

    //Kept for backward compatibility
    //Documentation see header -------------------------------------------------------------
    perl::Object divisorByValueVector(perl::Object fan, Vector<Rational> values,bool uses_min) {
      pm::cout << "divisorByValue is considered deprecated. Use \"divisor\" and \"function_value\" instead. For details please consult the documentation." << endl;
      Matrix<Rational> vmatrix(0,values.dim());
	vmatrix /= values;
      return divisorByValueMatrix(fan,vmatrix,uses_min);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see header
    perl::Object divisor_minmax(perl::Object complex, perl::Object function, int k) {
      //dbgtrace << "Preparing computations" << endl;
      
      //Homogenize the fan and refine it
      bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
      if(!uses_homog) complex = complex.CallPolymakeMethod("homogenize");
      perl::Object linearityDomains = function.give("NORMAL_FAN");
      
      //dbgtrace << "Refining fan" << endl;
      RefinementResult r = refinement(complex,linearityDomains,false,false,true,true,true);
      
      //Extract values
      Matrix<Rational> rays = r.complex.give("CMPLX_RAYS");	
      Matrix<Rational> linspace = r.complex.give("LINEALITY_SPACE");
//       Array<Set<int> > maximal = r.complex.give("CMPLX_MAXIMAL_CONES");
      Vector<int> assocRep = r.associatedRep;
      
      Matrix<Rational> fmatrix = function.give("FUNCTION_MATRIX");
      bool uses_min = function.give("USES_MIN");
      int power = function.give("POWER");
      
      power = (k < 1? power : k);
      
      //dbgtrace << "Extracted values" << endl;
      
      //Now compute function values
      Vector<Rational> values;    
      bool basepoint_found = false; //Save the first vertex
      Vector<Rational> basepoint;
      for(int r = 0; r < rays.rows(); r++) {
	//If it is an affine ray, simply compute the function value at that point
	if(rays(r,0) == 1) {
	    values |= functionValue(fmatrix, rays.row(r),uses_min,true);
	    if(!basepoint_found) {
	      basepoint_found = true;
	      basepoint = rays.row(r);
	    }
	}
	//Otherwise take the function difference between (x+this ray) and x for an associated vertex x
	else {
	    values |=  (functionValue(fmatrix, rays.row(assocRep[r]) + rays.row(r),uses_min,true) - 
			      functionValue(fmatrix, rays.row(assocRep[r]),uses_min,true));
	}
      }
      
      //Finally we add the function values on the lineality space
      for(int index = 0; index < linspace.rows(); index++) {
	values |= 
	    (functionValue(fmatrix, basepoint + linspace.row(index), uses_min,true) - 
	    functionValue(fmatrix,basepoint, uses_min,true));
// 	values |= functionValue(fmatrix, linspace.row(index), uses_min,true);
      }
      
      //Glue together to a value matrix
      Matrix<Rational> vmatrix(0,values.dim());
      for(int l = 1; l <= power; l++) {
	vmatrix /= values;
      }
      
      
      return divisorByValueMatrix(r.complex,vmatrix,uses_min);
      
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Kept for backward compatibility
    //Documentation see header -------------------------------------------------------------
    perl::Object divisorByPLF(perl::Object fan, perl::Object function) {
      pm::cout << "divisorByPLF is considered deprecated. Use \"divisor\" instead. For details please consult the documentation." << endl;
      return divisor_minmax(fan,function,1);
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    //Documentation see perl wrapper
    perl::Object divisor_nr(perl::Object complex, perl::Object function, int k) {
      //Extract function
      Vector<Rational> rvalues = function.give("RAY_VALUES");
      Vector<Rational> lvalues = function.give("LIN_VALUES");
      bool uses_min = function.give("USES_MIN");
      Vector<Rational> values = rvalues | lvalues;
      int power = function.give("POWER");

      //Create function matrix
      Matrix<Rational> fmatrix(0,values.dim());
      for(int l = 1; l <= (k==-1? power : k); l++) {
	fmatrix /= values;
      }
      
      return divisorByValueMatrix(complex,fmatrix,uses_min);
    }
 
    ///////////////////////////////////////////////////////////////////////////////////////
    
    perl::Object divisor_rational(perl::Object complex, perl::Object function, int k) {
      //Homogenize the fan if necessary and then refine it
      bool cmplx_uses_homog = complex.give("USES_HOMOGENEOUS_C");
      perl::Object domain = function.give("DOMAIN");
      bool fct_uses_homog = domain.give("USES_HOMOGENEOUS_C");
      Vector<Rational> rvalues = function.give("RAY_VALUES");
      Vector<Rational> lvalues = function.give("LIN_VALUES");
      bool uses_min = function.give("USES_MIN");
      int power = function.give("POWER");
      Vector<Rational> values = rvalues | lvalues;
      
      //Before we do anything, we check that the function is sane to avoid /0 divisions and segfaults
      Matrix<Rational> fnrays = domain.give("CMPLX_RAYS");
      Matrix<Rational> fnlin = domain.give("LINEALITY_SPACE");
      if(rvalues.dim() != fnrays.rows() || lvalues.dim() != fnlin.rows()) {
	  throw std::runtime_error("Function is not valid: Wrong number of values. Aborting computation. ");
      }
      
      if(fct_uses_homog && !cmplx_uses_homog) {
	complex = complex.CallPolymakeMethod("homogenize");
      }
      //dbgtrace << "Refining... " << endl;
      RefinementResult r = refinement(complex,domain, false,true,false,true,true);
      //dbgtrace << "Done. " << endl;
      
      
      //Compute the ray values on the new complex
      Matrix<Rational> rayRep = r.rayRepFromY;
      Matrix<Rational> linRep = r.linRepFromY;
      
      Vector<Rational> newvalues;
      for(int ray = 0; ray < rayRep.rows(); ray++) {
	newvalues |= (rayRep.row(ray) * values);
      }
      for(int lin = 0; lin < linRep.rows(); lin++) {
	newvalues |= (linRep.row(lin) * lvalues);
      }
      
      //Glue together to a matrix
      Matrix<Rational> fmatrix(0,newvalues.dim());
      for(int l = 1; l <= (k == -1? power : k); l++) {
	fmatrix /= newvalues;
      }
      
      return divisorByValueMatrix(r.complex, fmatrix,uses_min);
      
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////
    
    /**
      @brief Computes the RAY_VALUES and LIN_VALUES of a MinMaxFunction on its DOMAIN. Sets the corresponding properties automatically in the MinMaxFunction object
      @param perl::Object function
    */
    void computeMinMaxFunctionData(perl::Object function) {
      //Extract the functions data
      perl::Object domain = function.give("DOMAIN");
      Matrix<Rational> rays = domain.give("CMPLX_RAYS");
      Matrix<Rational> linspace = domain.give("LINEALITY_SPACE");
      Matrix<Rational> fmatrix = function.give("FUNCTION_MATRIX");
      bool uses_min = function.give("USES_MIN");
      
      //Compute the associated rays
      RefinementResult r = refinement(domain,domain,false,false,true,false);
      
      Vector<int> assocRep = r.associatedRep;
      //Now compute function values
      Vector<Rational> values;    
      bool basepoint_found = false; //Save the first vertex
      Vector<Rational> basepoint;
      for(int r = 0; r < rays.rows(); r++) {
	//If it is an affine ray, simply compute the function value at that point
	if(rays(r,0) == 1) {
	    values |= functionValue(fmatrix, rays.row(r),uses_min,true);
	    if(!basepoint_found) {
	      basepoint_found = true;
	      basepoint = rays.row(r);
	    }
	}
	//Otherwise take the function difference between (x+this ray) and x for an associated vertex x
	else {
	    values |=  (functionValue(fmatrix, rays.row(assocRep[r]) + rays.row(r),uses_min,true) - 
			      functionValue(fmatrix, rays.row(assocRep[r]),uses_min,true));
	}
      }
      function.take("RAY_VALUES") << values;
      Vector<Rational> linValues;
      
      //Finally we add the function values on the lineality space
      for(int index = 0; index < linspace.rows(); index++) {
	linValues |= (functionValue(fmatrix, basepoint + linspace.row(index), uses_min,true) - 
			  functionValue(fmatrix,basepoint, uses_min,true));
      }
      function.take("LIN_VALUES") << linValues;
    }
	
    ///////////////////////////////////////////////////////////////////////////////////////
	
    //Documentation see header
    perl::Object add_rational_functions(perl::Object f, perl::Object g) {
      //First, if any of the two is homogeneous and the other is not, we homogenize
      perl::Object fDomain = f.give("DOMAIN");
      perl::Object gDomain = g.give("DOMAIN");
      bool fhomog = fDomain.give("USES_HOMOGENEOUS_C");
      bool ghomog = gDomain.give("USES_HOMOGENEOUS_C");
      if(fhomog || ghomog) {
	if(!fhomog) {
	  f = f.CallPolymakeMethod("homogenize");
	  fDomain = f.give("DOMAIN");
	}
	if(!ghomog) {
	  g = g.CallPolymakeMethod("homogenize");
	  gDomain = g.give("DOMAIN");
	}
      }
      
      bool uses_min = f.give("USES_MIN");
      
      //Then compute the common refinement of the domains
      RefinementResult r = refinement(fDomain,gDomain,true,true,false,true);
	perl::Object nDomain = r.complex;
	Matrix<Rational> x_rayrep = r.rayRepFromX;
	Matrix<Rational> y_rayrep = r.rayRepFromY;
	Matrix<Rational> x_linrep = r.linRepFromX;
	Matrix<Rational> y_linrep = r.linRepFromY;
	
	Vector<Rational> f_rayval = f.give("RAY_VALUES");
	Vector<Rational> g_rayval = g.give("RAY_VALUES");
	Vector<Rational> f_linval = f.give("LIN_VALUES");
	Vector<Rational> g_linval = g.give("LIN_VALUES");
	
	Vector<Rational> fval = f_rayval | f_linval;
	Vector<Rational> gval = g_rayval | g_linval;
	
      Matrix<Rational> rays = nDomain.give("CMPLX_RAYS");
      Matrix<Rational> linspace = nDomain.give("LINEALITY_SPACE");
      
      //Now compute ray values
      Vector<Rational> rValues;
      for(int r = 0; r < rays.rows(); r++) {
	rValues |= (x_rayrep.row(r) * fval) + (y_rayrep.row(r) * gval);
      }
      //Now compute lin values
      Vector<Rational> lValues;
      for(int l = 0; l < linspace.rows(); l++) {
	lValues |= (x_linrep.row(l) * f_linval) + (y_linrep.row(l) * g_linval);
      }
      
      //Return result
      perl::Object func("RationalFunction");
	func.take("DOMAIN") << nDomain;
	func.take("RAY_VALUES") << rValues;
	func.take("LIN_VALUES") << lValues;
	func.take("USES_MIN") << uses_min;
	
      return func;
      
    }
    
// ------------------------- PERL WRAPPERS ---------------------------------------------------
    
    UserFunction4perl("# @category Basic polyhedral operations"
		      "# DEPRECATED. Use intersect_container instead." ,
		      &intersect_container,"intersect_complete_fan(WeightedComplex, fan::PolyhedralFan;$=0)");
    
    UserFunction4perl("# @category Basic polyhedral operations"
		      "# Takes two fans and computes the intersection of both. The function"
		      "# relies on the fact that the latter fan contains the first fan to "
		      "# compute the refinement correctly"
		      "# The function copies [[TROPICAL_WEIGHTS]] and [[LATTICE_BASES]]"
		      "# if they exist"
		      "# @param WeightedComplex fan An arbitrary weighted polyhedral fan"
		      "# @param fan::PolyhedralFan container A polyhedral fan containing the "
		      "# first one (as a set)"
		      "# @param Bool forceLatticeComputation Whether the properties"
		      "# [[LATTICE_BASES]] and [[LATTICE_GENERATORS]] of fan should be computed"
		      "# before refining. False by default."
		      "# @return WeightedComplex The intersection of both fans (whose support is equal to the support of fan). The "
		      "# resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, "
		      "# the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned." ,
		      &intersect_container,"intersect_container(WeightedComplex, fan::PolyhedralFan;$=0)");
    
    Function4perl(&divisorByValueVector,"divisorByValueVector(WeightedComplex, Vector<Rational>,$)");  
    
    Function4perl(&divisorByValueMatrix, "divisorByValueMatrix(WeightedComplex, Matrix<Rational>, $)");
    
    Function4perl(&divisor_minmax,"divisor_minmax(WeightedComplex,MinMaxFunction;$=1)");
    
    Function4perl(&divisor_rational, "divisor_rational(WeightedComplex,RationalFunction; $=1)");
    
    Function4perl(&computeMinMaxFunctionData,"computeMinMaxFunctionData(MinMaxFunction)");
    
    Function4perl(&add_rational_functions,"add_rational_functions(RationalFunction,RationalFunction)");
    
    UserFunction4perl("# @category Divisors"
		      "# NOTE: Deprecated. Use divisor(..) instead"
		      "# Computes the divisor of a MinMaxFunction on a given tropical variety. The result will be "
		      "# in homogeneous coordinates, whether the tropical variety uses them or not. The function "
		      "# should be given on the affine coordinates of the variety, NOT the homogeneous ones."
		      "# @param WeightedComplex fan A tropical variety, on which the divisor is computed"
		      "# @param MinMaxFunction A function whose DOMAIN should be equal to the affine coordinate "
		      "# space of the variety, i.e. AMBIENT_DIM-1, if the variety uses homogeneous coordinates, " 
		      "# AMBIENT_DIM otherwise."
		      "# @return The corresponding divisor as a tropical variety in homogeneous coordinates.",
		      &divisorByPLF, "divisorByPLF(WeightedComplex,MinMaxFunction)");
    
   UserFunction4perl("# @category Divisors"
		     "# Works exactly as divisor(WeightedComplex, RationalFunction;Int). Should be called ONLY,"
		     "# when the function f is defined on a DOMAIN equal to X (in the sense that all properties"
		     "# like RAYS, MAXIMAL_CONES, etc. agree. Being equal as varieties is not sufficient). In this"
		     "# case this function will in general be faster.",
		     &divisor_nr,"divisor_nr(WeightedComplex,RationalFunction;$=-1)");
   
   
}
}
