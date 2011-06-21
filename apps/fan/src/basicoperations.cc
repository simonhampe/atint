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

This file provides c++ realizations of basic operations on polyhedral complexes
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"

namespace polymake { namespace fan{ 
  
  //using namespace atint::donotlog;
  //using namespace atint::dolog;
  using namespace atint::dotrace;
  
  perl::Object product_complex(Array<perl::Object> complexes) {
    //This will contain the sets describing the maximal cones
    Vector<Set<int> > maximalCones;
      //Initially it will contain the empty set
      maximalCones = maximalCones | Set<int>();
    
    //Will contain the rays
    Matrix<Rational> rayMatrix;
    //Will contain the lineality space
    Matrix<Rational> linMatrix;
    //Will contain the weights (if any)
    Vector<Integer> weights;
    
    bool product_uses_homog = false;
    bool product_has_weights = false;
      
    dbgtrace << "Starting to create product" << endl;
    
    //Iterate over all complexes, each step creates the product of the first i complexes
    for(int i = 0; i < complexes.size(); i++) {
      dbgtrace << "Considering factor " << i << endl;
      //Extract properties
      bool uses_homog = complexes[i].give("USES_HOMOGENEOUS_C");
      bool uses_weights = false;
      Matrix<Rational> prerays = complexes[i].give("RAYS");
      Matrix<Rational> prelin = complexes[i].give("LINEALITY_SPACE");
      IncidenceMatrix<> premax = complexes[i].give("MAXIMAL_CONES");
      Array<Integer> preweights;
      if(complexes[i].exists("TROPICAL_WEIGHTS")) {
	uses_weights = true;
	product_has_weights = true;
	preweights = complexes[i].give("TROPICAL_WEIGHTS");
      }
      
      dbgtrace << "Computing rays" << endl;
      
      //Compute new rays and lineality space
      int additionalColumns = uses_homog ? prerays.cols() -1 : prerays.cols();
      int currentNumberOfRays = rayMatrix.rows();
      //Step 1: Append additional zeros to the ray matrix and lin matrix
      rayMatrix = rayMatrix | Matrix<Rational>(rayMatrix.rows(), additionalColumns);
      linMatrix = linMatrix | Matrix<Rational>(linMatrix.rows(), additionalColumns);
      //Step 2: Add a homogeneous coordinate in front if necessary
      if(uses_homog && !product_uses_homog) {
	rayMatrix = zero_vector<Rational>(rayMatrix.rows())| rayMatrix;
	linMatrix = zero_vector<Rational>(linMatrix.rows())| linMatrix;
	product_uses_homog = true;
      }
      //Step 3: Insert a zero's matrix in prerays /prelin between the first and the remaining columns
      // if it uses homogeneous coords, otherwise just prepend
      if(uses_homog) {
	prerays = prerays.col(0) | Matrix<Rational>(prerays.rows(), rayMatrix.cols() - additionalColumns -1) | prerays;
	prelin = prelin.col(0) | Matrix<Rational>(prelin.rows(), linMatrix.cols() - additionalColumns -1) | prelin;
      }
      else {
	prerays = Matrix<Rational>(prerays.rows(), rayMatrix.cols() - additionalColumns) | prerays;
	prelin = Matrix<Rational>(prelin.rows(), linMatrix.cols() - additionalColumns) | prelin;
      }
      //Step 4: Glue together
      rayMatrix = rayMatrix / prerays;
      linMatrix = linMatrix / prelin;
      
      dbgtrace << "New rays are : " << rayMatrix << endl;
      dbgtrace << "New lineality space is : " << linMatrix << endl;
      
      //Compute new cone sets 
      //Step 1: Replace each set by the set shifted by the currentNumberOfRays
      Vector<Set<int> > shiftedSets;
      for(int maxcone = 0; maxcone < premax.rows(); maxcone++) {
	Set<int> newset;
	Set<int> oldset = premax.row(maxcone);
	for(Entire<Set<int> >::iterator element = entire(oldset); !element.at_end(); ++element) {
	    newset += (*element + currentNumberOfRays);
	}
	shiftedSets = shiftedSets | newset;
      }
      //Step 2: Iterate over all pairs of old sets and new sets and build the union. Also compute weights
      //here
      Vector<Set<int> > coneReplacer;
      Vector<Integer> weightReplacer;
      for(int old = 0; old < maximalCones.dim(); old++) {
	for(int n = 0; n < shiftedSets.dim(); n++) {
	    coneReplacer = coneReplacer | (maximalCones[old] + shiftedSets[n]);
	    //If one of them has weights, add them
	    if(product_has_weights || uses_weights) {
	      weightReplacer = weightReplacer | ( (product_has_weights? weights[old] : Integer(0)) +
						   (uses_weights? preweights[n] : Integer(0)));
	    }
	}
      }
      maximalCones = coneReplacer;
      if(product_has_weights || uses_weights) {
	product_has_weights = true;
	weights = weightReplacer;
      }
      
      dbgtrace << "New cones are : " << maximalCones << endl;
      if(product_has_weights) {
	dbgtrace << "New weights are: " << weights << endl;
      }
    }
       
    
    perl::Object result("fan::PolyhedralFan");
      result.take("RAYS") << rayMatrix;
      result.take("MAXIMAL_CONES") << maximalCones;
      result.take("LINEALITY_SPACE") << linMatrix;
      result.take("USES_HOMOGENEOUS_C") << product_uses_homog;
      if(product_has_weights) result.take("TROPICAL_WEIGHTS") << weights;
      
    return result;
    
  }
  
  Function4perl(&product_complex,"compute_product_complex(;@)");
  
}}