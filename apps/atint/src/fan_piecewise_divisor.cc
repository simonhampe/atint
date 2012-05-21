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
 * This file computes divisors of piecewise polynomials on simplicial fans
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/basicoperations.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object piecewise_divisor(perl::Object fan, IncidenceMatrix<> cones, Vector<Integer> coefficients) {
    
    //Basic security checks
    if(cones.rows() != coefficients.dim()) 
	throw std::runtime_error("Cannot compute divisor of piecewise polynomial: Number of cones does not match number of coefficients");
    
    //Compute fan dimension
    int fan_dim = fan.give("CMPLX_DIM");
    
    //First we compute the appropriate skeleton of fan
    if(cones.rows() == 0) return fan;
    int result_dim = fan_dim - cones.row(0).size();
    perl::Object skeleton = skeleton_complex(fan,result_dim,true);
    
    //Extract values of skeleton
    Matrix<Rational> sk_rays = skeleton.give("RAYS");
    IncidenceMatrix<> sk_cones = skeleton.give("MAXIMAL_CONES");
    
    //This will contain the weights of the cones in the linear combination
    Vector<Integer> result_weights = zero_vector<Integer>(sk_cones.rows());
    
    //Now go through the divisors psi_tau for all cones tau in cones
    for(int tau = 0; tau < cones.rows(); tau++) {
      if(coefficients[tau] != 0) {
	//Create function matrix
	Matrix<Rational> psi_tau(0,sk_rays.rows());
	Set<int> tau_set = cones.row(tau);
	if(tau_set.size() != fan_dim - result_dim) 
	  throw std::runtime_error("Cannot compute divisor of piecewise polynomials: Cones have different dimension.");
	for(Entire<Set<int> >::iterator ts = entire(tau_set); !ts.at_end(); ts++) {
	  psi_tau /= unit_vector<Rational>(sk_rays.rows(),*ts);
	}
	
	//Compute divisor
	perl::Object divisor = divisorByValueMatrix(fan,psi_tau);
	
	//Extract cones, rays and weights
	Matrix<Rational> div_rays = divisor.give("RAYS");
	IncidenceMatrix<> div_cones = divisor.give("MAXIMAL_CONES");
	Vector<Integer> div_weights = divisor.give("TROPICAL_WEIGHTS");
	
	//Associate to each ray its original index 
	Map<int,int> div_ray_to_old;
	for(int dr = 0; dr < div_rays.rows(); dr++) {
	  for(int oray = 0; oray < sk_rays.rows(); oray++) {
	      if(sk_rays.row(oray) == div_rays.row(dr)) {
		div_ray_to_old[dr] = oray; break;
	      }
	  }
	}//END translate ray indices
	
	//Now go through all d-dimensional cones in the divisor and insert their weight at the appropriate point
	for(int rho = 0; rho < div_cones.rows(); rho++) {
	  //Map rho rays to old rays
	  Set<int> rho_old = 
	      attach_operation(div_cones.row(rho), pm::operations::associative_access<Map<int,int>, int>(&div_ray_to_old));
	  //Find the original cone equal to that
	  for(int oc = 0; oc < sk_cones.rows(); oc++) {
	    if( (sk_cones.row(oc) * rho_old).size() == result_dim) {
	      result_weights[oc] += (coefficients[tau] * div_weights[rho]); 
	      break;
	    }
	  }
	}//END iterate divisor cones
      }//END if coeff !=0
    }//END iterate cones 
    
    //Clean up by removing weight zero cones
    Set<int> used_cones;
    for(int c = 0; c < result_weights.dim(); c++) {
      if(result_weights[c] != 0) used_cones += c;
    }
    
    Set<int> used_rays = accumulate(rows(sk_cones.minor(used_cones,All)),operations::add());
    sk_rays = sk_rays.minor(used_rays,All);
    sk_cones = sk_cones.minor(used_cones,used_rays);
    result_weights = result_weights.slice(used_cones);
    
    //Return result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << sk_rays;
      result.take("MAXIMAL_CONES") << sk_cones;
      result.take("TROPICAL_WEIGHTS") << result_weights;
      
   return result;
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Intersection theory "
		    "# Computes a divisor of a linear sum of certain piecewise polynomials on a simplicial "
		    "# fan. "
		    "# @param WeightedComplex F A simplicial fan without lineality space in non-homog."
		    "# coordinates"
		    "# @param IncidenceMatrix cones A list of cones of F (not maximal, but all of the same "
		    "# dimension). Each cone t corresponds to a piecewise polynomial psi_t, defined by "
		    "# subsequently applying the rational functions that are 1 one exactly one ray of t and "
		    "# 0 elsewhere. "
		    "# @param Vector<Integer> coefficients A list of coefficients a_tcorresponding to the "
		    "# cones. "
		    "# @return WeightedComplex The divisor sum_t a_t psi_t * F",
		    &piecewise_divisor, "piecewise_divisor(WeightedComplex, IncidenceMatrix, Vector<Integer>)");
  
}}