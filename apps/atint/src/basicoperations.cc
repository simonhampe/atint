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
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint{ 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  //Documentation see header
  inline std::pair<Set<int>, Set<int> > separateRays(Matrix<Rational> m, Set<int> &affine, Set<int> &directional, bool uses_homog) {
    affine.clear();
    directional.clear();
    if(!uses_homog) {
      directional = sequence(0,m.rows());
      //We use the dummy index -1 to indicate that a homogenized fan has just the origin as affine ray
      affine += -1;
      return std::pair<Set<int>,Set<int> >(affine,directional);
    }
    
    for(int row = 0; row < m.rows(); row++) {
      if(m(row,0) == 0) {
	directional += row;
      }
      else {
	affine += row;
      }
    }
    return std::pair<Set<int>, Set<int> >(affine,directional);
  }
  
  /**
    @brief Takes a list of WeightedComplex objects (that may be weighted and that may use homogeneous coordinates) and computes the cartesian product of these. If any complex uses homogeneous coordinates, so will the result. If any complex has weights, all non-weighted complexes will be treated as having constant weight 1.
    @return The cartesian product of the complexes
  */
  perl::Object compute_product_complex(std::vector<perl::Object> complexes) {
      dbgtrace << "Generating container variables for result" << endl;
    
      perl::Object firstComplex = complexes[0];
      //This will contain the sets describing the maximal cones
      Vector<Set<int> > maximalCones = firstComplex.give("MAXIMAL_CONES");
      //Will contain the rays
      Matrix<Rational> rayMatrix = firstComplex.give("RAYS");
      //Will contain the lineality space
      Matrix<Rational> linMatrix = firstComplex.give("LINEALITY_SPACE");
      //Will contain the weights (if any)
      Vector<Integer> weights;
      bool product_has_weights = false;
      if(firstComplex.lookup("TROPICAL_WEIGHTS") >> weights) {
	product_has_weights = true;
      }
      bool product_uses_homog = firstComplex.give("USES_HOMOGENEOUS_C");
      int product_dim = rayMatrix.cols() > linMatrix.cols() ? rayMatrix.cols() : linMatrix.cols();
      //Sort rays by affine and directional
      Set<int> product_affine;
      Set<int> product_directional;
      separateRays(rayMatrix, product_affine, product_directional, product_uses_homog);
      
      dbgtrace << "Iterating over all " << complexes.size() -1 << " complexes: " << endl;
      
      for(unsigned int i = 1; i < complexes.size(); i++) {
	dbgtrace << "Considering complex nr. " << i+1 << endl;
	//Extract properties
	
	bool uses_homog = complexes[i].give("USES_HOMOGENEOUS_C");
	bool uses_weights = false;
	Matrix<Rational> prerays = complexes[i].give("RAYS");
	Matrix<Rational> prelin = complexes[i].give("LINEALITY_SPACE");
	IncidenceMatrix<> premax = complexes[i].give("MAXIMAL_CONES");
	
	Array<Integer> preweights;
	if(complexes[i].exists("TROPICAL_WEIGHTS")) {
	  uses_weights = true;
	  preweights = complexes[i].give("TROPICAL_WEIGHTS");
	}
	//Sort rays
	Set<int> complex_affine;
	Set<int> complex_directional;
	separateRays(prerays,complex_affine, complex_directional,uses_homog);
	//If this fan uses homog. coordinates, strip away the first column of rays and linear space
	if(uses_homog) {
	  if(prerays.rows() > 0) prerays = prerays.minor(All,~scalar2set(0));
	  if(prelin.rows() > 0) prelin = prelin.minor(All,~scalar2set(0));
	}
	int dim = prerays.cols() > prelin.cols() ? prerays.cols() : prelin.cols();
		
	dbgtrace << "Creating ray matrix" << endl;
	//dbgtrace << "Affine rays of product are " << rayMatrix.minor(product_affine,All) << endl;
	//dbgtrace << "Affine rays of complex are " << prerays.minor(complex_affine,All) << endl;
	
	//Create new ray matrix
	Matrix<Rational> newRays(0,product_dim + dim);
	//First create affine rays
	Map<std::pair<int,int>,int> affineIndices; //For index conversion
	if(product_uses_homog || uses_homog) {
	  for(Entire<Set<int> >::iterator prays = entire(product_affine); !prays.at_end(); prays++)  {
	      Vector<Rational> pRay;
	      if(*prays >= 0) pRay = rayMatrix.row(*prays);
	      else pRay = zero_vector<Rational>(product_dim);
	      for(Entire<Set<int> >::iterator crays = entire(complex_affine); !crays.at_end(); crays++) {
		Vector<Rational> cRay;
		if(*crays >= 0) cRay = prerays.row(*crays);
		else cRay = zero_vector<Rational>(dim);
		newRays = newRays / (pRay | cRay);
		affineIndices[std::make_pair(*prays,*crays)] = newRays.rows()-1;
	      }
	  }
	}
	Set<int> newAffine = sequence(0, newRays.rows());
	
	dbgtrace << "New affine rays read " << newRays << endl;
	
	//If the product is non-homog. and the complex is homog., we have to add a column of ones in front
	//We then also add a zero column in front of the old product rays for the directional rays
	if(!product_uses_homog && uses_homog && newRays.rows() > 0) {
	  newRays = ones_vector<Rational>(newRays.rows()) | newRays;
	  product_dim++;
	  if(rayMatrix.rows() > 0) rayMatrix = zero_vector<Rational>(rayMatrix.rows()) | rayMatrix;
	  if(linMatrix.rows() > 0) linMatrix = zero_vector<Rational>(linMatrix.rows()) | linMatrix;
	}
	//Now add the directional rays of both cones
	Map<int,int> pdirIndices;
	Map<int,int> cdirIndices; //For index conversion
	Vector<Rational> product_zero = zero_vector<Rational>(product_dim);
	Vector<Rational> complex_zero = zero_vector<Rational>(dim);
	for(Entire<Set<int> >::iterator prays = entire(product_directional); !prays.at_end(); prays++)  {
	  newRays = newRays / (rayMatrix.row(*prays) | complex_zero);
	  pdirIndices[*prays] = newRays.rows()-1;
	}
	for(Entire<Set<int> >::iterator crays = entire(complex_directional); !crays.at_end(); crays++) {
	  newRays = newRays / (product_zero | prerays.row(*crays));
	  cdirIndices[*crays] = newRays.rows()-1;
	}
	Set<int> newDirectional = sequence(newRays.rows()-1,newRays.rows() - product_affine.size());
	
	dbgtrace << "Creating lineality matrix" << endl;
	
	//Create new lineality matrix
	if(prelin.rows() > 0) {
	  prelin = Matrix<Rational>(prelin.rows(),product_dim) | prelin;
	}
	if(linMatrix.rows() > 0) {
	  linMatrix = linMatrix | Matrix<Rational>(linMatrix.rows(), dim);
	}
	
	dbgtrace << "Prelin = " << prelin << "\nlinMatrix = " << linMatrix << endl;
	
	//Copy values
	rayMatrix = newRays;
	linMatrix = linMatrix.rows() == 0? prelin : (prelin.rows() == 0? linMatrix : linMatrix / prelin);
	product_dim = rayMatrix.cols() > linMatrix.cols() ? rayMatrix.cols() : linMatrix.cols();
		
	dbgtrace << "Creating cones" << endl;
	
	//Now create the new cones and weights:
	Vector<Set<int> > newMaxCones;
	Vector<Integer> newWeights;
	//Make sure, we have at least one "cone" in each fan, even if it is empty
	if(maximalCones.dim() == 0) { maximalCones = maximalCones | Set<int>();}
	if(premax.rows() == 0) { premax = premax / Set<int>();}
	for(int pmax = 0; pmax < maximalCones.dim(); pmax++) {
	    Set<int> product_cone = maximalCones[pmax];
	    if(!product_uses_homog && uses_homog) {
	      product_cone = product_cone + (-1);
	    }
	    for(int cmax = 0; cmax < premax.rows(); cmax++) {
	      Set<int> complex_cone = premax.row(cmax);
	      if(!uses_homog && product_uses_homog) {
		complex_cone = complex_cone + (-1);
	      }
	      Set<int> newcone;
	      Set<int> pAffine = product_cone * product_affine;
	      Set<int> pDirectional = product_cone * product_directional;
	      Set<int> cAffine = complex_cone * complex_affine;
	      Set<int> cDirectional = complex_cone * complex_directional;
	      //First add the affine rays: For each pair of affine rays add the corresponding index from
	      //affineIndices
	      for(Entire<Set<int> >::iterator pa = entire(pAffine); !pa.at_end(); pa++) {
		for(Entire<Set<int> >::iterator ca = entire(cAffine); !ca.at_end(); ca++) {
		    newcone = newcone + affineIndices[std::make_pair(*pa,*ca)];
		}		
	      }
	      //Now add the directional indices
	      for(Entire<Set<int> >::iterator pd = entire(pDirectional); !pd.at_end(); pd++) {
		newcone = newcone + pdirIndices[*pd];
	      }
	      for(Entire<Set<int> >::iterator cd = entire(cDirectional); !cd.at_end(); cd++) {
		newcone = newcone + cdirIndices[*cd];
	      }
	      newMaxCones = newMaxCones | newcone;
	      //Compute weight
	      if(product_has_weights || uses_weights) {
		newWeights = newWeights | (product_has_weights? weights[pmax] : Integer(1)) * (uses_weights? preweights[cmax] : Integer(1));
	      }
	    }
	}
	
	//Copy values
	product_affine = newAffine;
	product_directional = newDirectional;
	maximalCones = newMaxCones;
	weights = newWeights;
	product_uses_homog = product_uses_homog || uses_homog;
	product_has_weights = product_has_weights ||  uses_weights;
      }
    
      //Fill fan with result
      perl::Object result("WeightedComplex");
	result.take("RAYS") << rayMatrix;
	result.take("LINEALITY_SPACE") << linMatrix;
	result.take("MAXIMAL_CONES") << maximalCones;
	if(product_has_weights) result.take("TROPICAL_WEIGHTS") << weights;
	result.take("USES_HOMOGENEOUS_C") << product_uses_homog;
	
      return result;
         
  }
  
  /**
   @brief Takes the ray matrix of a complex and returns two sets, describing which rows are affine and which are directional rays
   @param Matrix<Rational> m The ray matrix of the complex
   @param bool uses_homog Whether the complex uses homogeneous coordinates (if not, the result is trivial)
   @return perl::ListReturn A list containing two sets. In this order: Affine rays (vertices), directional rays   
   */
  perl::ListReturn separateRayMatrix(Matrix<Rational> m, bool uses_homog) {
    Set<int> affine, directional;
    separateRays(m,affine, directional,uses_homog);
    if(!uses_homog) {
      affine = Set<int>();
    }
    perl::ListReturn result;
      result << affine;
      result << directional;
    return result;
  }
  
  Function4perl(&separateRayMatrix,"separateRayMatrix(Matrix<Rational>,$)");
  
  Function4perl(&compute_product_complex,"compute_product_complex(;@)");

}}