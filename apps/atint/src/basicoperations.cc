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
#include "polymake/atint/specialvarieties.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/WeightedComplexRules.h"

namespace polymake { namespace atint{ 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  inline std::pair<Set<int>, Set<int> > separateRays(Matrix<Rational> m, Set<int> &affine, Set<int> &directional, bool uses_homog) {
    affine.clear();
    directional.clear();
    if(!uses_homog) {
      if(m.rows() > 0) directional = sequence(0,m.rows());
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
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
      Vector<Set<int> > local_restriction = firstComplex.give("LOCAL_RESTRICTION");
      //int product_dim = rayMatrix.cols() > linMatrix.cols() ? rayMatrix.cols() : linMatrix.cols();
      int product_dim = rayMatrix.rows() > 0? rayMatrix.cols() : linMatrix.cols();
      //Sort rays by affine and directional
      Set<int> product_affine;
      Set<int> product_directional;
      separateRays(rayMatrix, product_affine, product_directional, product_uses_homog);
	dbgtrace << "Product affine " << product_affine << endl;
	dbgtrace << "Product directional " << product_directional << endl;
      
      
      dbgtrace << "Iterating over all " << complexes.size() -1 << " complexes: " << endl;
      
      for(unsigned int i = 1; i < complexes.size(); i++) {
	dbgtrace << "Considering complex nr. " << i+1 << endl;
	//Extract properties
	
	bool uses_homog = complexes[i].give("USES_HOMOGENEOUS_C");
	bool uses_weights = false;
	Matrix<Rational> prerays = complexes[i].give("RAYS");
	Matrix<Rational> prelin = complexes[i].give("LINEALITY_SPACE");
	Vector<Set<int> > premax = complexes[i].give("MAXIMAL_CONES");
	Vector<Set<int> > pre_local_restriction = complexes[i].give("LOCAL_RESTRICTION");
	
	Array<Integer> preweights;
	if(complexes[i].exists("TROPICAL_WEIGHTS")) {
	  uses_weights = true;
	  preweights = complexes[i].give("TROPICAL_WEIGHTS");
	}
	//Sort rays
	Set<int> complex_affine;
	Set<int> complex_directional;
	separateRays(prerays,complex_affine, complex_directional,uses_homog);
	  dbgtrace << "Affine: " << complex_affine << endl;
	  dbgtrace << "Directional: " << complex_directional << endl;
	//If this fan uses homog. coordinates, strip away the first column of rays and linear space
	if(uses_homog) {
	  if(prerays.rows() > 0) prerays = prerays.minor(All,~scalar2set(0));
	  if(prelin.rows() > 0) prelin = prelin.minor(All,~scalar2set(0));
	  dbgtrace << "prerays: " << prerays << endl;
	  dbgtrace << "rcols: " << prerays.cols() << endl;
	  dbgtrace << "lcols: " << prelin.cols() << endl;
	}
	//int dim = prerays.cols() > prelin.cols() ? prerays.cols() : prelin.cols();
	int dim = prerays.rows() > 0? prerays.cols() : prelin.cols();
	dbgtrace << "dim: " << dim << endl;
		
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
	
	dbgtrace << "Adding directional rays" << endl;
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
	Set<int> newDirectional = sequence(newRays.rows()-1,product_directional.size() + complex_directional.size());
	
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
	if(premax.dim() == 0) { premax = premax | Set<int>();}
	for(int pmax = 0; pmax < maximalCones.dim(); pmax++) {
	    Set<int> product_cone = maximalCones[pmax];
	    if(!product_uses_homog && uses_homog) {
	      product_cone = product_cone + (-1);
	    }
	    for(int cmax = 0; cmax < premax.dim(); cmax++) {
	      Set<int> complex_cone = premax[cmax];
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
	
	//Compute the cross product of the local_restrictions
	Vector<Set<int> > new_local_restriction;
	if(local_restriction.dim() > 0 || pre_local_restriction.dim() > 0) {
	  //If one variety is not local, we take all its vertices for the product
	  Vector<Set<int> > product_locality(local_restriction);
	  if(product_locality.dim() == 0) {
	    for(Entire<Set<int> >::iterator aRay = entire(product_affine); !aRay.at_end(); aRay++) {
	      Set<int> single; single += *aRay;
	      product_locality |= single;
	    }
	  }
	  Vector<Set<int> > pre_locality(pre_local_restriction);
	  if(pre_locality.dim() == 0) {
	    for(Entire<Set<int> >::iterator aRay = entire(complex_affine); !aRay.at_end(); aRay++) {
	      Set<int> single; single += *aRay;
	      pre_locality |= single;
	    }
	  }
	  
	  for(int i = 0; i < product_locality.dim(); i++) {
	    Set<int> pAffine = product_locality[i] * product_affine;
	    Set<int> pDirectional = product_locality[i] * product_directional;
	    for(int j = 0; j < pre_locality.dim(); j++) {
	      Set<int> local_cone;
	      Set<int> cAffine = pre_locality[j] * complex_affine;
	      Set<int> cDirectional = pre_locality[j] * complex_directional;
	      //First add the affine rays: For each pair of affine rays add the corresponding index from
	      //affineIndices
	      for(Entire<Set<int> >::iterator pa = entire(pAffine); !pa.at_end(); pa++) {
		for(Entire<Set<int> >::iterator ca = entire(cAffine); !ca.at_end(); ca++) {
		    local_cone = local_cone + affineIndices[std::make_pair(*pa,*ca)];
		}		
	      }
	      //Now add the directional indices
	      for(Entire<Set<int> >::iterator pd = entire(pDirectional); !pd.at_end(); pd++) {
		local_cone = local_cone + pdirIndices[*pd];
	      }
	      for(Entire<Set<int> >::iterator cd = entire(cDirectional); !cd.at_end(); cd++) {
		local_cone = local_cone + cdirIndices[*cd];
	      }
	      new_local_restriction = new_local_restriction | local_cone;
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
	local_restriction = Vector<Set<int> > (new_local_restriction);
      }
    
      //Fill fan with result
      perl::Object result("WeightedComplex");
	result.take("RAYS") << rayMatrix;
	result.take("LINEALITY_SPACE") << linMatrix;
	result.take("MAXIMAL_CONES") << maximalCones;
	if(product_has_weights) result.take("TROPICAL_WEIGHTS") << weights;
	result.take("USES_HOMOGENEOUS_C") << product_uses_homog;
	result.take("LOCAL_RESTRICTION") << local_restriction;
	
      return result;
         
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  perl::Object facetRefinement(perl::Object fan, Matrix<Rational> facets) {
    for(int r = 0; r < facets.rows(); r++) {
      fan = intersect_complete_fan(fan,halfspace_complex(0,facets.row(r)));
    }
    return fan;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  perl::ListReturn fan_decomposition(perl::Object complex) {
    //Extract values
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    if(!uses_homog) {
      perl::ListReturn singlefan;
	singlefan << complex;
      return singlefan;
    }
    Matrix<Rational> rays = complex.give("RAYS");
    IncidenceMatrix<> cones = complex.give("MAXIMAL_CONES");
      IncidenceMatrix<> verticesInCones = T(cones);
    Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
    Set<int> affine = complex.give("VERTICES");
    Set<int> directional = complex.give("DIRECTIONAL_RAYS");
    int ambient_dim = complex.give("CMPLX_AMBIENT_DIM");
    Vector<Set<int> > local_restriction = complex.give("LOCAL_RESTRICTION");
    
    if(local_restriction.dim() > 0) {
      affine *= accumulate(local_restriction,operations::add());
    }
    
    
    //Now go through all affine rays
    perl::ListReturn result;
    for(Entire<Set<int> >::iterator v = entire(affine); !v.at_end(); v++) {
      //For each cone containing this vertex, we get a cone of the Star at v
      Set<int> conesatv = verticesInCones.row(*v);
      //Whether a ray has been added to the Star (for an affine ray, this actually means the distance
      Vector<bool> hasBeenAdded(rays.rows()); 
      //Contains the new row indices of the rays
      Map<int,int> newIndex;
      Matrix<Rational> starRays(0,ambient_dim+1);
      Vector<Set<int> > starCones;
      Vector<Integer> starWeights;
      //Now transform all the cones containing v
      for(Entire<Set<int> >::iterator catv = entire(conesatv); !catv.at_end(); catv++) {
	starWeights |= weights[*catv];
	Set<int> otherrays = cones.row(*catv) - *v;
	Set<int> newcone;
	for(Entire<Set<int> >::iterator r = entire(otherrays); !r.at_end(); r++) {
	    //If the ray is already in the matrix, just add its new index
	    if(hasBeenAdded[*r]) {
	      newcone += newIndex[*r];
	    }
	    else {
	      hasBeenAdded[*r] = true;
	      newIndex[*r] = starRays.rows();
	      newcone += starRays.rows();
	      //A directional ray is just copied, for an affine one we take the difference to v
	      if(directional.contains(*r)) {
		starRays /= rays.row(*r);
	      }
	      else {
		starRays /= (rays.row(*r) - rays.row(*v));
	      }
	    }
	}
	starCones |= newcone;
      } //End iterate all cones at v
      starRays = starRays.minor(All,~scalar2set(0));
      perl::Object fan("WeightedComplex");
	fan.take("INPUT_RAYS") << starRays;
	fan.take("INPUT_CONES") << starCones;
	fan.take("TROPICAL_WEIGHTS") << starWeights;
	fan.take("USES_HOMOGENEOUS_C") << false;
      result << fan;
    }//End iterate all vertices
    
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  perl::Object affineTransformation(perl::Object complex, Vector<Rational> translate, Matrix<Integer> matrix) {
    //Extract values
    Matrix<Rational> rays = complex.give("RAYS");
    Set<int> affine = complex.give("VERTICES");
    Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
    IncidenceMatrix<> cones = complex.give("MAXIMAL_CONES");
    Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    Vector<Set<int> > local_restriction = complex.give("LOCAL_RESTRICTION");
    
    //Transform rays and lin space
    linspace = linspace * matrix;
    rays = rays * matrix;
    for(Entire<Set<int> >::iterator r = entire(affine); !r.at_end(); r++) {
      rays.row(*r) = rays.row(*r) + translate;
    }
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("LINEALITY_SPACE") << linspace;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
      result.take("LOCAL_RESTRICTION") << local_restriction;
    return result;
    
  }

  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  perl::Object skeleton_complex(perl::Object complex, int k, bool preserve = false) {
    //Extract properties
    int cmplx_dim = complex.give("CMPLX_DIM");
    Matrix<Rational> rays = complex.give("RAYS");
    IncidenceMatrix<> maximalCones = complex.give("MAXIMAL_CONES");
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    Matrix<Rational> lineality = complex.give("LINEALITY_SPACE");
    int lineality_dim = complex.give("LINEALITY_DIM");
    Vector<Set<int> > local_restriction = complex.give("LOCAL_RESTRICTION");
    
    //If the skeleton dimension is too small, return the 0-cycle
    if(k < (uses_homog? 0 : 1) || k < lineality_dim) {
      return CallPolymakeFunction("zero_cycle");
    }
    
    //If the skeleton dimension is the fans dimension, return the fan
    if(k == cmplx_dim) {
      return complex;
    }
    
    Vector<Set<int> > new_local_restriction;
    
    //Now we compute the codimension one skeleton of the fan (cmplx_dim - k) times 
    IncidenceMatrix<> newMaximalCones = maximalCones;
    for(int i = 1; i <= (cmplx_dim - k); i++) {
      newMaximalCones = calculateCodimOneData(rays,newMaximalCones, uses_homog, lineality, local_restriction).codimOneCones;
    }
    
    //Now return the result - made irredundant, if preserve is false
    Matrix<Rational> newrays = rays;
    if(!preserve) {
      //Take the union of all cones to see what rays are uses
      Set<int> usedRays;
      for(int c = 0; c < newMaximalCones.rows(); c++) {
	usedRays += newMaximalCones.row(c);
      }
      
      if(local_restriction.dim() > 0) {
	  Map<int,int> index_map; //Maps indices of old rays to indices of new rays
	  int newIndex = 0;
	  for(Entire<Set<int> >::iterator uR = entire(usedRays); !uR.at_end(); uR++) {
	      index_map[*uR] = newIndex;
	      newIndex++;
	  }
	  for(int i = 0; i < local_restriction.dim(); i++) {
	      new_local_restriction |= attach_operation(local_restriction[i] * usedRays, pm::operations::associative_access<Map<int,int>,int>(&index_map));
	  }
      }
      
      
      
      newrays = newrays.minor(usedRays,All);
      newMaximalCones = newMaximalCones.minor(All,usedRays);
    }
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << newrays;
      result.take("MAXIMAL_CONES") << newMaximalCones;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
      result.take("LINEALITY_SPACE") << lineality;
      result.take("LOCAL_RESTRICTION") << new_local_restriction;
   
    return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&separateRayMatrix,"separateRayMatrix(Matrix<Rational>,$)");
  
  Function4perl(&compute_product_complex,"compute_product_complex(;@)");
  
  Function4perl(&facetRefinement,"facetRefinement(WeightedComplex,Matrix<Rational>)");
  
  Function4perl(&affineTransformation, "affineTransformation(WeightedComplex, Vector<Rational>, Matrix<Integer>)");
  
  Function4perl(&skeleton_complex,"calculate_skeleton_complex(WeightedComplex,$;$=1)");
  
  UserFunction4perl("# @category Tropical geometry"
		    "# Take a polyhedral complex and returns a list of all the local vertex fans, "
		    "# i.e. for each affine ray r, the list contains the fan Star_complex(r) "
		    "# (in non-homogeneous coordinates)"
		    "# @param WeightedComplex complex A tropical variety"
		    "# @return perl::ListReturn A list of WeightedComplex objects in "
		    "# non-homogeneous coordinates. The i-th complex corresponds to the i-th "
		    "# affine ray ( vertex). If the complex is not in homogeneous coordinates, "
		    "# the list contains just the complex itself ",
    &fan_decomposition,"fan_decomposition(WeightedComplex)");
  
  
}}