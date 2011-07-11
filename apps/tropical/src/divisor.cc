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
#include "polymake/tropical/divisor.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"


namespace polymake { namespace tropical {

    using polymake::polytope::cdd_interface::solver;
  
    using namespace atint::donotlog;
    //using namespace atint::dolog;
    //using namespace atint::dotrace;
    
    perl::ListReturn intersection_test(perl::Object c1, perl::Object c2) {
      Matrix<Rational> facets = c1.give("FACETS | INEQUALITIES");
      Matrix<Rational> f2 = c2.give("FACETS | INEQUALITIES");
	facets /= f2;
      Matrix<Rational> linear = c1.give("LINEAR_SPAN | EQUATIONS");
      Matrix<Rational> l2 =  c2.give("LINEAR_SPAN | EQUATIONS");
	linear /= l2;
      std::pair<Matrix<Rational>, Matrix<Rational> > s = solver<Rational>().enumerate_vertices(facets,linear,true);
      
      perl::ListReturn result;
	result << s.first;
	result << s.second;      
      return result;
    }
    
    //Documentation see header
    perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan) {
      
      //First we compute ambient dimensions and check that they are equal
      int ambient_dim =  fan.CallPolymakeMethod("ambient_dim_fix"); //FIXME: Replace
      int comp_ambient_dim = completeFan.CallPolymakeMethod("ambient_dim_fix"); //FIXME: Replace
            
      if(ambient_dim != comp_ambient_dim) {
	throw std::runtime_error("Cannot intersect fans in different ambient dimension.");
      }
      
      //If fan is zero-dimensional it is just a point, so it has no maximal cones.
      //But since no refinement is necessary anyway, we just return fan
      int dimension = fan.CallPolymakeMethod("dim_fix");
      if(dimension == 0) return fan;
      
      //Extract values
      Matrix<Rational> rays = fan.give("RAYS");
      Matrix<Rational> lineality_space = fan.give("LINEALITY_SPACE");
      Array<Integer> weights;
      bool weightsExist = false;
      if(fan.exists("TROPICAL_WEIGHTS")) {
	weights = fan.give("TROPICAL_WEIGHTS");
	weightsExist = true;
      }
      IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
      bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
      
      Matrix<Rational> compRays = completeFan.give("RAYS");
      IncidenceMatrix<> compMaximalCones = completeFan.give("MAXIMAL_CONES");
      Matrix<Rational> compLinealitySpace = completeFan.give("LINEALITY_SPACE");
      
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
      
      //Now intersect all pairs of maximal cones.
      for(int fanIndex = 0; fanIndex < fanUpperBound; fanIndex++) {
	for(int cIndex = 0; cIndex < compMaximalCones.rows(); cIndex++) {
	    //Create the cones and compute the intersection
	    perl::Object fanCone("polytope::Cone<Rational>");
	    perl::Object compCone("polytope::Cone<Rational>");
	      fanCone.take("RAYS") << (fanHasOnlyLineality?
					Matrix<Rational>(1,ambient_dim) :
					rays.minor(maximalCones.row(fanIndex),All));
	      fanCone.take("LINEALITY_SPACE") << fan.give("LINEALITY_SPACE");
	      compCone.take("RAYS") << compRays.minor(compMaximalCones.row(cIndex),All);
	      compCone.take("LINEALITY_SPACE") << compLinealitySpace;
	    perl::Object inter = CallPolymakeFunction("intersection",fanCone,compCone);
	    int coneDimension = inter.give("CONE_DIM");
	    
	    //Only consider intersections that have the correct dimension
	    if(coneDimension == dimension) {
	      //The first time we compute this, the lineality space of the intersection is copied as the lineality space of
		//the resulting fan
		if(!computedLineality) {
		  inter.give("LINEALITY_SPACE") >> newLineality;
		  computedLineality = true;
		  dbgtrace << "Setting lineality space to " << newLineality << endl;
		  
		}
		
		//Now go through the rays of the intersection to find those already added
		Matrix<Rational> interRays = inter.give("RAYS");
		
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
		if(remainingRows.size() < 0) {
		  addCone = true;
		}
		else {
		  addCone = !(newConeSet.contains(newMaximalCone));
		}
		
		if(addCone) {
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
      perl::Object result("fan::PolyhedralFan");
	result.take("RAYS") << newRays;
	result.take("MAXIMAL_CONES") << newCones;
	result.take("USES_HOMOGENEOUS_C") << uses_homog;
	result.take("LINEALITY_SPACE") << newLineality;
	if(weightsExist) result.take("TROPICAL_WEIGHTS") << newWeights;
      return result;
      
      
    }
    
    perl::Object divisorByValueVector(perl::Object fan, Vector<Rational> values) {
      //Extract properties
      bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
      Matrix<Rational> rays = uses_homog? fan.give("CMPLX_RAYS") : fan.give("RAYS");
      int lineality_dim = fan.give("LINEALITY_DIM");
      Matrix<Rational> lineality_space = fan.give("LINEALITY_SPACE");
      Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
      IncidenceMatrix<> codimOneCones = uses_homog? fan.give("CMPLX_CODIM_1_FACES") : fan.give("CODIM_1_FACES");
      IncidenceMatrix<> coneIncidences = fan.give("CODIM_1_IN_MAXIMAL_CONES");
      Array<Integer> tropicalWeights = fan.give("TROPICAL_WEIGHTS");
      Map<int, Map<int, Vector<Rational> > > lnFunctionVector = fan.give("LATTICE_NORMAL_FCT_VECTOR");
      Matrix<Rational> lsumFunctionVector = fan.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
      
      //If the fan only consists of the lineality space, the divisor is empty
      int noOfRays = rays.rows();
      if(noOfRays == 0) {
	return perl::Object("fan::PolyhedralFan");
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
	    newcones = newcones | codimOneCones.row(co);
	    newweights = newweights | coweight;
	}
      }
      
      dbgtrace << "Done\n" << endl;
      
      //Return result
      perl::Object result("fan::PolyhedralFan");
	if(uses_homog) {
	    result.take("INPUT_HOM_RAYS") << rays;
	}
	else {
	    result.take("INPUT_RAYS") << rays;
	}
	result.take("INPUT_CONES") << newcones;
	result.take("TROPICAL_WEIGHTS") << newweights;
	result.take("LINEALITY_SPACE") << lineality_space; 
	//FIXME: Why do I need this? It seems, INPUT_RAYS does not get computed properly otherwise
	result.give("INPUT_RAYS");
      return result;
    }
    
    //perl::Object divisorByPLF(perl::Object fan, perl::Object function)
    
// ------------------------- PERL WRAPPERS ---------------------------------------------------
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete "
		      "# (i.e. its support is the whole ambient space) to compute the intersection correctly."
		      "# @param fan::PolyhedralFan fan An arbitrary polyhedral fan"
		      "# @param fan::PolyhedralFan completeFan A complete polyhedral fan"
		      "# @return fan::PolyhedralFan The intersection of both fans (whose support is equal to the support of fan). The "
		      "# resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, "
		      "# the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned." ,
		      &intersect_complete_fan,"intersect_complete_fan(fan::PolyhedralFan, fan::PolyhedralFan)");
    
    Function4perl(&divisorByValueVector,"divisorByValueVector(fan::PolyhedralFan, Vector<Rational>)");  
    
    Function4perl(&intersection_test,"intersection_test(polytope::Cone, polytope::Cone)");
}
}