#include "polymake/client.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Array.h"

namespace polymake { namespace tropical {

    //using namespace atint::donotlog;
    //using namespace atint::dolog;
    using namespace atint::dotrace;
    
    //Documentation see header
    perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan) {
      
      Integer fanAmbientDim = fan.CallPolymakeMethod("ambient_dim_fix"); //FIXME: Replace both when possible
      Integer compAmbientDim = completeFan.CallPolymakeMethod("ambient_dim_fix");
      
      //First check that ambient dimensions coincide
      if(fanAmbientDim != compAmbientDim) {
	dbgtrace << "Fan ambient dim is " << fanAmbientDim << endl;
	dbgtrace << "Complete fan ambient dim is " << compAmbientDim << endl;
	
	throw std::runtime_error("Cannot intersect fans in different ambient dimension.");
      }
      
      
      //If fan is zero-dimensional it is just a point, so it has no maximal cones.
      //But since no refinement is necessary anyway, we just return fan
      //if(fan.give("DIM") == 0) return fan;
      
      //Prepare containers for new values
      Matrix<Rational> rays(0,fanAmbientDim);
      Matrix<Rational> linealitySpace;
	bool computedLineality = false;
      Vector<Set<int> > maximalCones;
      Set<Set<int> > maximalConeSet; //This is kept to check for doubles. It has the same content as maximalCones
      Vector<int> tropicalWeights(0);
      
      //Read out old values
      Matrix<Rational> fanOldrays = fan.give("RAYS");
      Matrix<Rational> compOldrays = completeFan.give("RAYS");
      Array<Integer> oldWeights;
      bool weightsExist = false;
	if(fan.exists("TROPICAL_WEIGHTS")) {
	  oldWeights = fan.give("TROPICAL_WEIGHTS");
	  weightsExist = true;
	}
      IncidenceMatrix<> fanMaximals = fan.give("MAXIMAL_CONES");
      IncidenceMatrix<> compMaximals = completeFan.give("MAXIMAL_CONES");
      Integer fanDimension = fan.CallPolymakeMethod("dim_fix"); //FIXME: Replace this when possible
            
      //If the complete fan is only a lineality space, it is the whole space, so we return fan
      if(compMaximals.rows() == 0) {
	return fan;
      }
      
      //If the fan, however, is only a lineality space, we have to take care that it is considered
      // as a maximal cone of fan (note that at this point, the lineality dimension must then be > 0,
      // since the fan dimension is > 0
      bool fanHasOnlyLineality = fanMaximals.rows() == 0;
      int fanUpperBound = fanHasOnlyLineality? 1 : fanMaximals.rows();
     
      dbgtrace << "Starting to intersect maximal cone pairs" << endl;
      
      //Now intersect all pairs of maximal cones
      for(int tindex = 0; tindex < fanUpperBound; tindex ++) {
	for(int cindex = 0; cindex < compMaximals.rows(); cindex ++) {
	  
	  dbgtrace << "Intersecting cones " << tindex << " and " << cindex << endl;
	  
	  //Create the actual cones and compute the intersection
	    perl::Object tcone("polytope::Cone<Rational>");
	    perl::Object ccone("polytope::Cone<Rational>");
	      tcone.take("RAYS") << (fanHasOnlyLineality? 
				Matrix<Rational>(1,fanAmbientDim) : 
				fanOldrays.minor(fanMaximals.row(tindex),All));
				//Have to give it a ray to tell polymake about the ambient dimension
	      ccone.take("RAYS") << compOldrays.minor(compMaximals.row(cindex),All);
	      tcone.take("LINEALITY_SPACE") << fan.give("LINEALITY_SPACE");
	      ccone.take("LINEALITY_SPACE") << completeFan.give("LINEALITY_SPACE");

	    perl::Object inter = CallPolymakeFunction("intersection",tcone,ccone);
	    Integer coneDimension = inter.give("CONE_DIM");
	    
	    dbgtrace << "Cone dimension is " << coneDimension << endl;
	    
	    //Check if they intersect in the correct dimension
	    if(coneDimension == fanDimension) {
		
		dbgtrace << "Dimensions match. Checking rays" << endl;
	      
		//The first time we compute this, the lineality space of the intersection is copied as the lineality space of
		//the resulting fan
		if(!computedLineality) {
		  inter.give("LINEALITY_SPACE") >> linealitySpace;
		  computedLineality = true;
		  
		  dbgtrace << "Setting lineality space to " << linealitySpace << endl;
		  
		}
			
		//Now read off the rays and find the ones that have already been added to rays
		Matrix<Rational> interRays = inter.give("RAYS");
		
		dbgtrace << "Ray matrix is " << interRays << endl;
		
		Set<int> remainingRows;
		Set<int> newMaximalCone;
		//Go through all rays of the new cone
		for(int row = 0; row < interRays.rows(); row++) {
		  int raysIndex = -1;
		  //Compare it to all existing rays
		  for(int testRay = 0; testRay < rays.rows(); testRay ++) {
		    if(rays.row(testRay) == interRays.row(row)) {
			raysIndex = testRay;
			break;
		    }
		  }
		  //If we found one, add the index to the maximal cone, otherwise keep track of the row to create new rays later
		  if(raysIndex >= 0) {
		    newMaximalCone += raysIndex;
		  }
		  else {
		    remainingRows += row;
		  }
		}
		
		dbgtrace << "Remaining rows are " << remainingRows << ", Already existing are " << newMaximalCone << endl;
		
 		//Now add remaining rays and the corresponding cone (if it doesn't exist yet)
		bool addCone = false;
		if(remainingRows.size() > 0) addCone = true;
		else {
		  addCone = !(maximalConeSet.contains(newMaximalCone));
		}
		
		dbgtrace << "Now we " << (addCone? "" : "don't ") << "add the cone." << endl;
		
		if(addCone) {
		  
		  dbgtrace << "Cone does not exist. Adding it" << endl;
		  
		  rays = rays / interRays.minor(remainingRows,All);
		  for(int i = rays.rows() - remainingRows.size(); i < rays.rows(); i++) {
		    newMaximalCone = newMaximalCone + i;
		  }
		  maximalCones = maximalCones | newMaximalCone;
		  maximalConeSet = maximalConeSet + newMaximalCone;
		  if(weightsExist) {
		    tropicalWeights = tropicalWeights | oldWeights[tindex];
		  }	      
		}
	    }
	    
	}	
      }//End iteration over all maximal cone pairs
            
      //Return result
      perl::Object result("fan::PolyhedralFan");
	result.take("RAYS") << rays;
	result.take("MAXIMAL_CONES") << maximalCones;
	result.take("USES_HOMOGENEOUS_C") << fan.give("USES_HOMOGENEOUS_C");
	result.take("LINEALITY_SPACE") << linealitySpace;
	if(weightsExist) result.take("TROPICAL_WEIGHTS") << tropicalWeights;
      return result;
      
    }
    
    UserFunction4perl("# @category Tropical geometry"
		      "# Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete "
		      "# (i.e. its support is the whole ambient space) to compute the intersection correctly."
		      "# @param fan::PolyhedralFan fan An arbitrary polyhedral fan"
		      "# @param fan::PolyhedralFan completeFan A complete polyhedral fan"
		      "# @return fan::PolyhedralFan The intersection of both fans (whose support is equal to the support of fan). The "
		      "# resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, "
		      "# the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned." ,
		      &intersect_complete_fan,"intersect_complete_fan(fan::PolyhedralFan, fan::PolyhedralFan)");
    
}
}