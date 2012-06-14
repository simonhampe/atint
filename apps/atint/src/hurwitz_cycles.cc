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
 * This file contains functions to compute Hurwitz cycles in M_0,n
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/psi_classes.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/morphism_special.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_pre_cycle(int k, Vector<int> degree, Vector<int> pullback_points = Vector<int>()) {
    
    //First, compute the psi-class product
    int n = degree.dim();
    
    if(pullback_points.dim() < n-3-k) {
      pullback_points |= zero_vector<int>(n-3-k - pullback_points.dim());
    }
    
    int big_moduli_dim = 2*n - k - 2;
    Vector<int> exponents = zero_vector<int>(n) | ones_vector<int>(n-2-k);
    //dbgtrace << "Computing psi product in M_N, N = " << big_moduli_dim << " with exponents " << exponents << endl;
    perl::Object P = psi_product(big_moduli_dim,exponents);
    
    //Then compute evaluation map pullbacks (in each case of the pullback point p_i cut out by max(x,p_i))
    Matrix<Rational> rat_degree(degree.dim(),0);
      rat_degree |= degree;
    
    for(int i = n+2; i <= 2*n - 2 - k; i++) {
      //dbgtrace << "Computing evaluation map pull back for i = " << i-n-1 << endl;
      perl::Object evi = evaluation_map(n-2-k, 1, rat_degree, i-n-1);
      Matrix<Rational> evi_matrix = evi.give("MATRIX");
      
      //Pulling back p_i = max(x,p_i) * R means we take the vector representing the morphism and 
      //attach a row below that has p_i at the end
      evi_matrix /= zero_vector<Rational>(evi_matrix.cols());
      evi_matrix(1, evi_matrix.cols()-1) = pullback_points[i-n-2];
      //Since we restrict ourselves to M_0,N x {0}, we actually ignore the last coefficient
      //of ev_i and replace it by the constant coefficient 0
      evi_matrix(0, evi_matrix.cols()-1) = 0;
      //dbgtrace << "Pullback evaluation matrix is " << evi_matrix << endl;
      perl::Object pb("MinMaxFunction");
	pb.take("FUNCTION_MATRIX") << Matrix<Rational>(evi_matrix);
	pb.take("USES_MIN") << false;

		
      P = divisor_minmax(P, pb);
    }
    
    return P;
    
  }//END function hurwitz_pre_cycle
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_curve(Vector<int> degree) {
    //First we compute the pre-cycle
    perl::Object precycle = hurwitz_pre_cycle(1, degree);
    
    //Then we extract values
    Matrix<Rational> rays = precycle.give("RAYS");
    IncidenceMatrix<> cones = precycle.give("MAXIMAL_CONES");
    Vector<Integer> weights = precycle.give("TROPICAL_WEIGHTS");
    bool uses_homog = precycle.give("USES_HOMOGENEOUS_C");
    
    //Result variables
    Matrix<Rational> curve_rays(0,rays.cols());
    Vector<Set<int> > curve_cones;
    Vector<Integer> curve_weights;
    
    //First of all, we find the homogenizing ray and remove it
    if(uses_homog) {
      int homog_index = -1;
      for(int r = 0; r < rays.rows(); r++) {
	if(rays(r,0) == 1) {
	  homog_index = r; break;
	}
      }
      rays = rays.minor(~scalar2set(homog_index),~scalar2set(0));
      cones = cones.minor(All,~scalar2set(homog_index));
    }
    
    //Then we apply the forgetful map to the rays
    perl::Object ffmap = CallPolymakeFunction("forgetful_map",2*degree.dim() -3, sequence(degree.dim()+1, degree.dim()-3));
    Matrix<Rational> ffmatrix = ffmap.give("MATRIX");
    rays = rays * T(ffmatrix);
    
    //Iterate all cones
    for(int c = 0; c < cones.rows(); c++) {
      //The cone now only consists of one ray
      int r_index = cones.row(c).front();
      //Check if this ray already exists
      int n_index = -1;
      for(int nray = 0; nray < curve_rays.rows(); nray++) {
	if(curve_rays.row(nray) == rays.row(r_index)) {
	  n_index = nray; break;
	}
      }
      
      //If it doesn't exist, add it
      if(n_index == -1) {
	curve_rays /= rays.row(r_index);
	Set<int> single_set;
	single_set += (curve_rays.rows()-1);
	curve_cones |= single_set;
	curve_weights |= weights[c];
      }
      //Otherwise add its weight to the appropriate cone 
      //(whose index is now equal to the ray index)
      else {
	curve_weights[n_index] += weights[c];
      }
      
    }//END iterate cones
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << curve_rays;
      result.take("MAXIMAL_CONES") << curve_cones;
      result.take("TROPICAL_WEIGHTS") << curve_weights;
      
    return result;
    
  }//END function hurwitz_curve
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Integer hurwitz_degree(Vector<int> degree) {
    //First we compute the pre-cycle
    perl::Object precycle = hurwitz_pre_cycle(1, degree);
    
    Vector<Integer> weights = precycle.give("TROPICAL_WEIGHTS");
    
    return accumulate(weights, operations::add()); 	
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
		    "# Computes the k-dimensional tropical Hurwitz cycle H_k(degree), except that it doesn't"
		    "# compute the push-forward to M_0,n"
		    "# @param int k The dimension of the Hurwitz cycle"
		    "# @param Vector<int> degree The degree of the covering. The sum over all entries should "
		    "# be 0 and if n := degree.dim, then 0 <= k <= n-3"
		    "# @param Vector<int> pullback_points The points p_i that should be pulled back to "
		    "# determine the Hurwitz cycle. Should have length n-3-k. If it is not given, all p_i"
		    "# are by default equal to 0 (same for missing points)"
		    "# @return perl::Object A WeightedComplex object representing the Hurwitz cycle H_k(degree) before push-forward",    
		    &hurwitz_pre_cycle, "hurwitz_pre_cycle($, Vector<Int>; Vector<Int> = new Vector<Int>())");
  
  UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
		    "# Computes the Hurwitz curve H_1(degree)"
		    "# @param Vector<int> degree The degree of the covering. The sum over all entries should "
		    "# be 0 and if n := degree.dim, then 0 <= 1 <= n-3"
		    "# @return perl::Object A WeightedComplex object representing the Hurwitz cycle H_1(degree). This will always be a fan cycle", 
		    &hurwitz_curve, "hurwitz_curve(Vector<Int>)");
  
  UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
		    "# Computes the Hurwitz degree H_0(degree)"
		    "# @param Vector<int> degree The degree of the covering. The sum over all entries should "
		    "# be 0 and if n := degree.dim, then 0 <= n-3"
		    "# @return Integer The Hurwitz degree H_0(degree)", 
		    &hurwitz_degree, "hurwitz_degree(Vector<Int>)");
  
}}