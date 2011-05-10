#include "polymake/client.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Array.h"

namespace polymake { namespace tropical {

    using namespace atint::donotlog;
    //using namespace atint::dolog;
    //using namespace atint::dotrace;
    
    //Documentation see header
    perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan) {
      
      //First check that ambient dimensions coincide
      if(fan.give("AMBIENT_DIM") != completeFan.give("AMBIENT_DIM")) {
	throw std::runtime_error("Cannot intersect fans in different ambient dimension.");
      }
      
      //If fan is zero-dimensional it is just a point, so it has no maximal cones.
      //But since no refinement is necessary anyway, we just return fan
      if(fan.give("DIM") == 0) return fan;
      
      //Prepare containers for new values
      Matrix<Rational> rays(0,fan.give("AMBIENT_DIM"));
      Matrix<Rational> linealitySpace;
	bool computedLineality = false;
      IncidenceMatrix<> maximalCones;
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
      
      
      //Now intersect all pairs of maximal cones
      for(int tindex = 0; tindex < fanMaximals.rows(); tindex ++) {
	for(int cindex = 0; cindex < fanMaximals.rows(); cindex ++) {
	    
	  //Create the actual cones and compute the intersection
	    perl::Object tcone("polytope::Cone<Rational>");
	    perl::Object ccone("polytope::Cone<Rational>");
	      tcone.take("RAYS") << fanOldrays.minor(fanMaximals.row(tindex),All);
	      ccone.take("RAYS") << compOldrays.minor(compMaximals.row(cindex),All);
	      tcone.take("LINEALITY_SPACE") << fan.give("LINEALITY_SPACE");
	      ccone.take("LINEALITY_SPACE") << completeFan.give("LINEALITY_SPACE");
	    perl::Object inter = CallPolymakeFunction("polytope::intersection",tcone,ccone);
	    //The first time we compute this, the lineality space of the intersection is copied as the lineality space of
	    //the resulting fan
	    if(!computedLineality) {
	      inter.give("LINEALITY_SPACE") >> linealitySpace;
	      computedLineality = true;
	    }
	    
	    //Now read off the rays and find the ones that have already been added to rays
	    Matrix<Rational> interRays = inter.give("RAYS");
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
	    
	    //Now add remaining rays and the corresponding cone (if it doesn't exist yet)
	    bool addCone = false;
	    if(remainingRows.size() > 0) addCone = true;
	    else {
	      addCone = !(maximalConeSet.contains(newMaximalCone));
	    }
	    if(addCone) {
	      
	    }
	    
	}	
      }
      
      
    }
    
}}