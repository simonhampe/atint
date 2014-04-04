/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 * 
 * ---
 * Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 * 
 * This file computes the stable intersection of two tropical varieties
 * using the Minkowski weight formula
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/PowerSet.h"
#include "polymake/Array.h"
#include "polymake/linalg.h"
#include "polymake/Smith_normal_form.h"
#include "polymake/RandomGenerators.h"
#include "polymake/atint/normalvector.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/cdd_helper_functions.h"
#include "polymake/atint/WeightedComplexRules.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
/*  
  //Documentation see perl wrapper
  Integer lattice_index(const Matrix<Integer> &lattice_rays) {
    if(lattice_rays.rows() < lattice_rays.cols()) {
      throw std::runtime_error("Cannot compute lattice index - not a full-dimensional lattice!");
    }
    
    //dbgtrace << "Computing minors" << endl;
    //dbgtrace << "Set is " << sequence(0, lattice_rays.rows()) << endl;
    //dbgtrace << "Number of elements is " << lattice_rays.cols() << endl;
    Array<Set<int> > minors = all_subsets_of_k(sequence(0,lattice_rays.rows()),lattice_rays.cols());
    Integer g = 0;
    //dbgtrace << "Computing gcd" << endl;
    for(int s = 0; s < minors.size(); s++) {
      g = gcd(g,abs(det(lattice_rays.minor(minors[s],All))));
    }
    
    return g;
    
  }*/

  ///////////////////////////////////////////////////////////////////////////////////////
 
  //Documentation see perl wrapper
  Integer lattice_index(const Matrix<Integer> &lattice_rays) {
    //Compute the Smith Normal form 
    SmithNormalForm<Integer> solution = smith_normal_form(lattice_rays);
    
    Integer result = 1;
    for(int i = 0; i < solution.rank; i++) {
      result *= solution.form(i,i);
    }
    
    return abs(result);
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the Star around a point in a given set of cones
   @param Vector<Rational> point The point around which the star is to be computed, given in homog. coordinates
   @param Matrix<Rational> rays The rays of the surrounding cones (needn't contain point) in homog. coordinates
   @param Vector<Set<int> > cones The surrounding cones (i.e. should all contain point). Should form a complex.
   @param Matrix<Rational> result_rays Will contain the rays of the result (in non-homog. coordinates). Might not be irredundant
   @param Vector<Set<int> > result_cones Will contain the cones of the result. A cone might be given by a redundant list of rays. The cones will be in the exact same order as in "cones", i.e. the i-th fan cone is the star of the i-th cone at point.
   */
  void computeStar(const Vector<Rational> &point, const Matrix<Rational> &rays, const IncidenceMatrix<> &cones,
		   Matrix<Rational> &result_rays, Vector<Set<int> > &result_cones) {
    //Prepare result variables
    result_rays = Matrix<Rational>(0,rays.cols()-1);
    result_cones = Vector<Set<int> >();
    
    
    Matrix<Rational> fan_rays = rays.minor(All,~scalar2set(0));
    Vector<Rational> fan_point = point.slice(~scalar2set(0));
    //Iterate all surrounding cones
    for(int sc = 0; sc < cones.rows(); sc++) {
      Set<int> cone;
      Set<int> surroundCone = cones.row(sc);
      for(Entire<Set<int> >::iterator r = entire(surroundCone); !r.at_end(); r++) {
	if(rays(*r,0) == 0) {
	    result_rays /= fan_rays.row(*r);
	}
	else {
	    result_rays /= (fan_rays.row(*r) - fan_point);
	}
	cone += (result_rays.rows()-1);
      }
      result_cones |= cone;
    }
    
    cdd_normalize_rays(result_rays,false);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the Minkowski multiplicity of two fans, i.e. find two cones that add up to full dimension. Then it chooses an interior vector in the difference and finds all cone differences containing it as an interior vector. For all these differences it adds the product of the weights times the lattice index of the sum of the lattices
   */
  Integer computeFanMultiplicity(const Matrix<Rational> &xrays, const Matrix<Rational> &xlin, 
				 const Vector<Set<int> > &xcones, const Vector<Integer> &xweights,
				 const Matrix<Rational> &yrays, const Matrix<Rational> &ylin, 
				 const Vector<Set<int> > &ycones, const Vector<Integer> &yweights) {
    Integer weight(0);
    solver<Rational> sv;
    
    // First, we compute all H-representations of xcone - ycone, keeping
    // only full-dimensional ones
    Vector< Matrix<Rational> > full_dimensional_cones;
    Vector<int> full_dimensional_xindex; //Keep track of associated x- and ycones
    Vector<int> full_dimensional_yindex;
    for(int xc = 0; xc < xcones.dim(); xc++) {
      for(int yc = 0; yc < ycones.dim(); yc++) {
	//dbgtrace << "Having cones " << xrays.minor(xcones[xc],All) << ", \n" << yrays.minor(ycones[yc],All) << endl;
	Matrix<Rational> x_sub_rays = xrays.minor(xcones[xc],All);
	Matrix<Rational> y_sub_rays = (- yrays.minor(ycones[yc],All));
	std::pair<Matrix<Rational>, Matrix<Rational> > eqs = 
	  sv.enumerate_facets(zero_vector<Rational>() | x_sub_rays / y_sub_rays, 
			      zero_vector<Rational>() | (xlin / ylin),true,false);
	if(eqs.second.rows() == 0){
	    //dbgtrace << "Is fulldimensional" << endl;
	    full_dimensional_cones |= eqs.first;
	    full_dimensional_xindex |= xc;
	    full_dimensional_yindex |= yc;
	}
      }
    }
    
    //If there are no full-dimensional cones, the result is 0
    if(full_dimensional_cones.dim() == 0) return weight;
    
    //Otherwise, we need to compute a generic vector. We compute a 
    // random vector and go through all cones. We add up appropriate
    // weights for those cones containing the point in their interior.
    // If we find one that contains the point in its boundary, we 
    // create another interior point and try again.
    bool point_found;
    UniformlyRandom<Rational> random_gen;
    Vector<Rational> interior_point(xrays.cols()+1);
    //dbgtrace << "Generating generic point" << endl;
    do {
      weight = Integer(0);
      copy(random_gen.begin(), entire(interior_point));
      interior_point[0] = 1;
      point_found = true;
      //dbgtrace << "Trying " << interior_point << endl;
      //Now go through all full-dimensional cones
      for(int fullcone = 0; fullcone < full_dimensional_cones.dim(); fullcone++) {
	//dbgtrace << "Checking fulldimension cone " << fullcone << endl;
	//dbgtrace << "Has facets " << full_dimensional_cones[fullcone] << endl;
	Vector<Rational> eq_check = full_dimensional_cones[fullcone] * interior_point;
	bool is_interior = true;
	bool is_in_boundary = false;
	for(int c = 0; c < eq_check.dim(); c++) {
	  if(eq_check[c] == 0) {
	    is_in_boundary = true; break;
	  }
	  if(eq_check[c] < 0) {
	    is_interior = false; break;
	  }
	}//END check for interiorness
	// If its in the boundary of something, try another point.
	if(is_in_boundary) {
	    //dbgtrace << "It is a boundary point. Trying another one..." << endl;
	    point_found = false; break;
	}
	//If its interior, add the appropriate weight.
	if(is_interior) {
	    //dbgtrace << "Is interior point of this cone, computing weight..." << endl;
	    //dbgtrace << "xweight: " << xweights[full_dimensional_xindex[fullcone]] << endl;
	    //dbgtrace << "yweight: " << yweights[full_dimensional_yindex[fullcone]] << endl;
	    Integer latticeIndex = lattice_index(
	      latticeBasisFromRays(
		xrays.minor(xcones[full_dimensional_xindex[fullcone]],All),xlin) / 
	      latticeBasisFromRays(
		yrays.minor(ycones[full_dimensional_yindex[fullcone]],All),ylin));
	    //dbgtrace << "lattice: " << latticeIndex<< endl;
	    weight += (xweights[full_dimensional_xindex[fullcone]] * yweights[full_dimensional_yindex[fullcone]] * latticeIndex);
	}
	
      }//END iterate full-dimensional cones
    } while(!point_found);
    
    
    return weight;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object minkowski_intersection(perl::Object X, perl::Object Y) {
    //Extract values
    int Xcodim = X.give("CMPLX_CODIMENSION");
    int Ycodim = Y.give("CMPLX_CODIMENSION");
    int Xambi  = X.give("CMPLX_AMBIENT_DIM");
    
    //dbgtrace << "Checking codimension" << endl;
    
    //If the codimensions of the varieties add up to something larger then CMPLX_AMBIENT_DIM, return the 0-cycle 
    if(Xcodim + Ycodim > Xambi) {
      return CallPolymakeFunction("zero_cycle");
    }
    
    //dbgtrace << "Homogenizing where necessary" << endl;
    
    //Make sure,both are homogeneous
    bool x_uses_homog = X.give("USES_HOMOGENEOUS_C");
    bool y_uses_homog = Y.give("USES_HOMOGENEOUS_C");
    
    if(!x_uses_homog) X = X.CallPolymakeMethod("homogenize");
    if(!y_uses_homog) Y = Y.CallPolymakeMethod("homogenize");
    
    //Extract values
    Matrix<Rational> xrays = X.give("RAYS");
    Matrix<Rational> xlin = X.give("LINEALITY_SPACE");
    IncidenceMatrix<> xcones = X.give("MAXIMAL_CONES");
    Vector<Integer> xweights = X.give("TROPICAL_WEIGHTS");
    int xambdim = X.give("CMPLX_AMBIENT_DIM");
    
    Matrix<Rational> yrays = Y.give("RAYS");
    Matrix<Rational> ylin = Y.give("LINEALITY_SPACE");
    IncidenceMatrix<> ycones = Y.give("MAXIMAL_CONES");
    Vector<Integer> yweights = Y.give("TROPICAL_WEIGHTS");
    int yambdim = Y.give("CMPLX_AMBIENT_DIM");
    
    if(xambdim != yambdim) {
      throw std::runtime_error("Cannot compute intersection product: Cycles live in different spaces.");
    }
   
    //Compute the expected dimension of the intersection product 
    int k = Xambi - (Xcodim + Ycodim);
    
    //Compute the intersection complex
    fan_intersection_result f = cdd_fan_intersection(xrays,xlin,xcones,yrays,ylin,ycones,true);
   
    //Now we compute the k-skeleton of the intersection complex together with the data of original maximal
    //cones containing these cones
    Matrix<Rational> interrays = f.rays;
    Matrix<Rational> interlin = f.lineality_space;
    int i_lineality_dim = rank(f.lineality_space);
    Vector<Set<int> > intercones;
    Vector<Set<int> > xcontainers;
    Vector<Set<int> > ycontainers;
    
    for(int ic = 0; ic < f.cones.rows(); ic++) {
      //Check that the cone dimension is at least the expected dimension
      int cone_dim = rank(interrays.minor(f.cones.row(ic),All)) + i_lineality_dim -1;
      if(cone_dim >= k) {
	//Now we compute the k-skeleton of the intersection cone
	Vector<Set<int> > singlecone; singlecone |= f.cones.row(ic);
	IncidenceMatrix<> k_skeleton_matrix(singlecone);
	for(int i = cone_dim; i > k; i--) {
	    k_skeleton_matrix = 
	      calculateCodimOneData(interrays, k_skeleton_matrix, true, interlin, IncidenceMatrix<>()).codimOneCones;
	}
	
	//Go through all cones and add them (if they haven't already been added)
	for(int kc = 0; kc < k_skeleton_matrix.rows(); kc++) {
	    int cone_index = -1;
	    for(int oc = 0; oc < intercones.dim(); oc++) {
	      //Since both cones have the same dimension, it suffices to check, whether the old cone
	      //is contained in the new cone
	      if((intercones[oc] * k_skeleton_matrix.row(kc)).size() == intercones[oc].size()) {
		cone_index = oc; break;
	      }
	    }
	    //If it doesn't exist yet, add it
	    if(cone_index == -1) {
	      intercones |= k_skeleton_matrix.row(kc);
	      xcontainers |= Set<int>();
	      ycontainers |= Set<int>();
	      cone_index = intercones.dim()-1;
	    }
	    
	    //Now add containers
	    xcontainers[cone_index] += f.xcontainers.row(ic);
	    ycontainers[cone_index] += f.ycontainers.row(ic);
	    
	}//END iterate all k-skeleton cones
	
      }//END if cone_dim >= k
    }//END iterate intersection cones
    
    //If no cones remain, return the zero cycle
    if(intercones.dim() == 0) {
      return CallPolymakeFunction("zero_cycle");
    }
     
    //dbgtrace << "Computing weights " << endl;
    
    //Now we compute weights
    Vector<Integer> weights(intercones.dim());
    Set<int> weight_zero_cones;
    
    Matrix<Rational> xlin_dehom = xlin.minor(All,~scalar2set(0));
    Matrix<Rational> ylin_dehom = ylin.minor(All,~scalar2set(0));
    
    for(int c = 0; c < intercones.dim(); c++) {
      //dbgtrace << "Computing on intersection cone " << c << endl;
      //Find interior point
      Vector<Rational> interior_point = accumulate(rows(interrays.minor(intercones[c],All)),operations::add());
      Rational count_vertices = accumulate(interrays.col(0).slice(intercones[c]),operations::add());
      if(count_vertices != 0) interior_point /= count_vertices;
      
      //dbgtrace << "Interior point is " << interior_point << endl;
      
      //Compute stars
      Matrix<Rational> xstar_rays, ystar_rays;
      Vector<Set<int> > xstar_cones, ystar_cones;
      
      //dbgtrace << "Computing stars " << endl;
      
      computeStar(interior_point, xrays, xcones.minor(xcontainers[c],All), xstar_rays, xstar_cones);
      computeStar(interior_point, yrays, ycones.minor(ycontainers[c],All), ystar_rays, ystar_cones);
      
      //dbgtrace << "X Star rays: " << xstar_rays << endl;
      //dbgtrace << "X Star cones: " << xstar_cones << endl;
      //dbgtrace << "Y Star rays: " << ystar_rays << endl;
      //dbgtrace << "Y Star cones: " << ystar_cones << endl;
      
      //dbgtrace << "Computing multiplicity " << endl;
      
      Integer w = computeFanMultiplicity(
	xstar_rays, xlin_dehom, xstar_cones, xweights.slice(xcontainers[c]),
	ystar_rays, ylin_dehom, ystar_cones, yweights.slice(ycontainers[c]));
      
      //dbgtrace << "Weight is " << w << endl;
      
      weights[c] = w;
      if(w == 0) weight_zero_cones += c;
      
    }
    
    //Check if any cones remain
    if(weight_zero_cones.size() == intercones.dim()) {
      return CallPolymakeFunction("zero_cycle");
    }
    
    //dbgtrace << "Done" << endl;
    
    //Clean up rays and cones
    
    IncidenceMatrix<> intercones_matrix(intercones);
      intercones_matrix = intercones_matrix.minor(~weight_zero_cones,All);
      Set<int> used_rays = accumulate(rows(intercones_matrix), operations::add());
      intercones_matrix = intercones_matrix.minor(All,used_rays);
      interrays = interrays.minor(used_rays,All);
      weights = weights.slice(~weight_zero_cones);
      
    //Finally create the result  
    perl::Object result("WeightedComplex");
      result.take("RAYS") << interrays;
      result.take("MAXIMAL_CONES") << intercones_matrix;
      result.take("USES_HOMOGENEOUS_C") << true;
      result.take("LINEALITY_SPACE") << interlin;
      result.take("TROPICAL_WEIGHTS") << weights;
    
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Integer minkowski_multiplicity(perl::Object X, perl::Object Y) {
    //Extract values
    bool x_uses_homog = X.give("USES_HOMOGENEOUS_C");
    bool y_uses_homog = Y.give("USES_HOMOGENEOUS_C");
    
    if(x_uses_homog || y_uses_homog) {
      throw std::runtime_error("Cannot compute Minkowski multiplicity. Need fans, not complexes");
    }
    
    Matrix<Rational> xrays = X.give("RAYS");
    Matrix<Rational> xlin = X.give("LINEALITY_SPACE");
    Vector<Set<int> > xcones = X.give("MAXIMAL_CONES");
    Vector<Integer> xweights = X.give("TROPICAL_WEIGHTS");
    int xambdim = X.give("CMPLX_AMBIENT_DIM");
    
    Matrix<Rational> yrays = Y.give("RAYS");
    Matrix<Rational> ylin = Y.give("LINEALITY_SPACE");
    Vector<Set<int> > ycones = Y.give("MAXIMAL_CONES");
    Vector<Integer> yweights = Y.give("TROPICAL_WEIGHTS");
    int yambdim = Y.give("CMPLX_AMBIENT_DIM");
    
    if(xambdim != yambdim) {
      throw std::runtime_error("Cannot compute Minkowski multiplicity: Fans live in different spaces.");
    }
    
    //Compute multiplicity
    return computeFanMultiplicity(xrays, xlin, xcones, xweights, 
				  yrays, ylin, ycones, yweights);
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Integer and lattice arithmetic"
		    "# This computes the index of a lattice in its saturation."
		    "# @param Matrix<Integer> m A list of (row) generators of the lattice."
		    "# @return Integer The index of the lattice in its saturation.",
		    &lattice_index,"lattice_index(Matrix<Integer>)");
  
  UserFunction4perl("# @category Intersection products"
		    "# Computes the intersection product of two tropical cycles in R^n"
		    "# @param WeightedComplex X A tropical cycle"
		    "# @param WeightedComplex Y A tropical cycle, living in the same space as X"
		    "# @return WeightedComplex The intersection product, always in homogeneous coordinates",
		    &minkowski_intersection,"intersect(WeightedComplex,WeightedComplex)");
  
  UserFunction4perl("# @category Intersection products"
		    "# Computes the minkowski multiplicity of two fans, i.e. it computes the unique weight of "
		    "# the Minkowski sum X + (-Y). More precisely: If the sum is not full-dimensional, it "
		    "# returns 0. Otherwise it runs over all pairs of cones whose Minkowski sum contains"
		    "# a generic vector and adds the product of their weight times the lattice index of the sum"
		    "# of the lattices."
		    "# @param WeightedComplex X A tropical fan"
		    "# @param WeightedComplex Y A tropical fan, living in the same space as X"
		    "# @return Integer The Minkowski multiplicity",
		    &minkowski_multiplicity, "minkowski_multiplicity(WeightedComplex, WeightedComplex)");
  
}}