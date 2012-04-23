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
  
  //Documentation see perl wrapper
  Integer lattice_index(const Matrix<Integer> &lattice_rays) {
    if(lattice_rays.rows() < lattice_rays.cols()) {
      throw std::runtime_error("Cannot compute lattice index - not a full-dimensional lattice!");
    }
    Array<Set<int> > minors = all_subsets_of_k(sequence(0,lattice_rays.rows()),lattice_rays.cols());
    Integer g = 0;
    for(int s = 0; s < minors.size(); s++) {
      g = gcd(g,abs(det(lattice_rays.minor(minors[s],All))));
    }
    
    return g;
    
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
    
    Vector<Rational> interior_point(xrays.cols());
    bool point_found = false;
    
    
    for(int xc = 0; xc < xcones.dim(); xc++) {
      for(int yc = 0; yc < ycones.dim(); yc++) {
	//dbgtrace << "Having cones " << xrays.minor(xcones[xc],All) << ", \n" << yrays.minor(ycones[yc],All) << endl;
	Matrix<Rational> x_sub_rays = xrays.minor(xcones[xc],All);
	Matrix<Rational> y_sub_rays = (- yrays.minor(ycones[yc],All));
	//Compute H-representation of xcone - ycone
	std::pair<Matrix<Rational>, Matrix<Rational> > eqs = 
	  sv.enumerate_facets(zero_vector<Rational>() | x_sub_rays / y_sub_rays, 
			      zero_vector<Rational>() | (xlin / ylin),true,false);
	  
	//We're only interested in full-dimensional cones
	if(eqs.second.rows() == 0) {
	    //dbgtrace << "Is fulldimensional" << endl;
	    //Compute an interior point, if necessary
	    if(!point_found) {
	      if(eqs.first.rows() == 0) {
		  interior_point = zero_vector<Rational>(xrays.cols() > xlin.cols()? xrays.cols() : xlin.cols());
	      }
	      else {
		Matrix<Rational> r = sv.enumerate_vertices(eqs.first,eqs.second,true,true).first;
		interior_point = accumulate(rows(r),operations::add());
		//dbgtrace << "Setting interior point to " << interior_point << endl;
		point_found = true;
	      }
	    }
	    //Otherwise check if this point is an interior point
	    else {
	      Vector<Rational> eq_check = eqs.first * interior_point;
	      bool is_interior = true;
	      for(int c = 0; c < eq_check.dim(); c++) {
		if(eq_check[c] <= 0) {
		  is_interior = false; break;
		}
	      }
	      if(!is_interior) continue;
	      //dbgtrace << "Contains interior point" << endl;
	    }//END check if is interior point
	    
	    //If we arrive here, compute weight
	    //dbgtrace << "xweight: " << xweights[xc] << endl;
	    //dbgtrace << "yweight: " << yweights[yc] << endl;
	    Integer latticeIndex = lattice_index(latticeBasisFromRays(x_sub_rays,xlin) / 
						  latticeBasisFromRays(y_sub_rays,ylin));
	    //dbgtrace << "lattice: " << latticeIndex<< endl;
	    weight += (xweights[xc] * yweights[yc] * latticeIndex);
	    
	}//END check if diff is full-dimensional
      }//END iterate ycones
    }//END iterate xcones
    
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
  
  UserFunction4perl("# @category Lattice arithmetic"
		    "# This computes the index of a lattice of rank n in Z^n"
		    "# @param Matrix<Integer> m A list of (row) generators of the lattice. The matrix must have"
		    "# full column rank, otherwise an error is thrown"
		    "# @return Integer The index of the lattice",
		    &lattice_index,"lattice_index(Matrix<Integer>)");
  
  UserFunction4perl("# @category Tropical geometry / Intersection theory"
		    "# Computes the intersection product of two tropical cycles in R^n"
		    "# @param WeightedComplex X A tropical cycle"
		    "# @param WeightedComplex Y A tropical cycle, living in the same space as X"
		    "# @return WeightedComplex The intersection product, always in homogeneous coordinates",
		    &minkowski_intersection,"intersect(WeightedComplex,WeightedComplex)");
  
  UserFunction4perl("# @category Tropical geometry / Intersection theory"
		    "# Computes the minkowski multiplicity of two fans, i.e. it computes the unique weight of "
		    "# the Minkowski sum X + (-Y) "
		    "# @param WeightedComplex X A tropical fan"
		    "# @param WeightedComplex Y A tropical fan, living in the same space as X"
		    "# @return Integer The Minkowski multiplicity",
		    &minkowski_multiplicity, "minkowski_multiplicity(WeightedComplex, WeightedComplex)");
  
}}