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
  perl::Object hurwitz_pre_cycle(int k, Vector<int> degree) {
    
    //First, compute the psi-class product
    int n = degree.dim();
    int big_moduli_dim = 2*n - k - 2;
    Vector<int> exponents = zero_vector<int>(n) | ones_vector<int>(n-2-k);
    //dbgtrace << "Computing psi product in M_N, N = " << big_moduli_dim << " with exponents " << exponents << endl;
    perl::Object P = psi_product(big_moduli_dim,exponents);
    
    //Then compute evaluation map pullbacks (in each case of the same point = 0 cut out by max(0,x))
    Matrix<Rational> rat_degree(degree.dim(),0);
      rat_degree |= degree;
    
    for(int i = n+2; i <= 2*n - 2 - k; i++) {
      //dbgtrace << "Computing evaluation map pull back for i = " << i-n-1 << endl;
      perl::Object evi = evaluation_map(n-2-k, 1, rat_degree, i-n-1);
      Matrix<Rational> evi_matrix = evi.give("MATRIX");
      
      //Pulling back 0 = max(x,0) * R means we take the vector representing the morphism and 
      //attach a zero row below
      evi_matrix /= zero_vector<Rational>(evi_matrix.cols());
      //Since we restrict ourselves to M_0,N x {0}, we actually ignore the last coefficient
      //and replace it by the constant coefficients (0,0)
      evi_matrix.col(evi_matrix.cols()-1) = zero_vector<Rational>(2);
      //dbgtrace << "Pullback evaluation matrix is " << evi_matrix << endl;
      perl::Object pb("MinMaxFunction");
	pb.take("FUNCTION_MATRIX") << Matrix<Rational>(evi_matrix);
	pb.take("USES_MIN") << false;

		
      P = divisor_minmax(P, pb);
    }
    
    return P;
    
  }//END function hurwitz_pre_cycle
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_curve(Vector<int> degree) {
    //First we compute the pre-cycle
    perl::Object precycle = hurwitz_pre_cycle(1, degree);
    
    //Then we extract values
    Matrix<Rational> rays = precycle.give("RAYS");
    IncidenceMatrix<> cones = precycle.give("MAXIMAL_CONES");
    Vector<Integer> weights = precycle.give("TROPICAL_WEIGHTS");
    
    //Result variables
    Matrix<Rational> curve_rays(0,rays.cols());
    Vector<Set<int> > curve_cones;
    Vector<Integer> curve_weights;
    
    //First of all, we find the homogenizing ray and remove it
    int homog_index = -1;
    for(int r = 0; r < rays.rows(); r++) {
      if(rays(r,0) == 1) {
	homog_index = r; break;
      }
    }
    rays = rays.minor(~scalar2set(homog_index),~scalar2set(0));
    cones = cones.minor(All,~scalar2set(homog_index));
    
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
      //Otherwise add its weigh to the appropriate cone 
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
		    "# @return perl::Object A WeightedComplex object representing the Hurwitz cycle H_k(degree) before push-forward",    
		    &hurwitz_pre_cycle, "hurwitz_pre_cycle($, Vector<Int>)");
  
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