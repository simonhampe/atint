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
 * This file contains some basic functions related to simplicial fans and their diagonals
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
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object simplicial_with_diagonal(perl::Object fan) {
    
    //Extract values
    Matrix<Rational> rays = fan.give("RAYS");
    IncidenceMatrix<> cones = fan.give("MAXIMAL_CONES");
    Vector<Integer> weights = fan.give("TROPICAL_WEIGHTS");
    
    //Compute new rays (first the first component, then the second, then the diagonal rays)
    Matrix<Rational> p_rays(3*rays.rows(), 2* rays.cols());
    p_rays.minor(sequence(0,rays.rows()), sequence(0, rays.cols())) = 
      p_rays.minor(sequence(rays.rows(), rays.rows()), sequence(rays.cols(), rays.cols())) = 
      p_rays.minor(sequence(2*rays.rows(), rays.rows()), sequence(0,rays.cols())) =
      p_rays.minor(sequence(2*rays.rows(), rays.rows()), sequence(rays.cols(), rays.cols())) = rays;
    
    //Compute translation maps for ray indices
    Map<int,int> second_comp;
    Map<int,int> diagon_comp;
    for(int r = 0; r < rays.rows(); r++) {
      second_comp[r] = r+rays.rows();
      diagon_comp[r] = r+ 2* rays.rows();
    }
      
    Vector<Set<int> > p_cones;
    Vector<Integer> p_weights;
    //Go through pairs of cones
    for(int c1 = 0; c1 < cones.rows(); c1++) {
      Set<int> cone1 = cones.row(c1);
      for(int c2 = 0; c2 < cones.rows(); c2++) {
	  Set<int> cone2 = cones.row(c2);
	  Set<int> inter = cone1 * cone2;
	  //dbgtrace << "Cones " << cone1 << ", " << cone2 << " intersect in " << inter.size() << " elements " << endl;
	  //If the cones share no ray, their product is not subdivided
	  if(inter.size() == 0) {
	      Set<int> cone2_translate = 
		attach_operation(cone2, pm::operations::associative_access<Map<int,int>, int>(&second_comp));
	      p_cones |= (cone1 + cone2_translate);
	      p_weights |= weights[c1] * weights[c2];
	  }
	  else {
	    //Otherwise we get a subdivision cone for each subset of the intersection
	    Array<Set<int> > subdivisions = pm::all_subsets(inter);
	    //dbgtrace << "Replaced by " << subdivisions.size() << " cones " << endl;
	    for(int sd = 0; sd < subdivisions.size(); sd++) {
		Set<int> complement = inter - subdivisions[sd];
		//We replace the elements from the subset of the intersection in the FIRST cone
		//by the corresponding diagonal rays and replace the elements of the rest of the
		//intersection in the SECOND cone by the corresponding diagonal rays
		Set<int> first_cone_remaining = cone1 - subdivisions[sd];
		Set<int> second_cone_remaining = cone2 - complement;
		  second_cone_remaining = 
		    attach_operation(second_cone_remaining,
				     pm::operations::associative_access<Map<int,int>, int>(&second_comp));
		
		Set<int> first_cone_diagonal = 
		    attach_operation(subdivisions[sd], 
				     pm::operations::associative_access<Map<int,int>, int>(&diagon_comp));
		Set<int> second_cone_diagonal = 
		    attach_operation(complement,
				     pm::operations::associative_access<Map<int,int>, int>(&diagon_comp));
		    
		Set<int> final_cone = first_cone_remaining + second_cone_remaining + 
				      first_cone_diagonal + second_cone_diagonal;
		
		p_cones |= final_cone;
		p_weights |= weights[c1] * weights[c2];
		    
	    }//END iterate subdivisions
	  }//END compute subdivision cones
	  
      }//END iterate second cone
    }//END iterate first cone
    
    //Return result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << p_rays;
      result.take("MAXIMAL_CONES") << p_cones;
      result.take("TROPICAL_WEIGHTS") << p_weights;
    return result;
  }//END function simplicial_with_diagonal
  
   ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Matrix<Rational> simplicial_piecewise_system(perl::Object fan) {
    
    //Extract values
    int fan_dim = fan.give("CMPLX_DIM");
    
    //Compute diagonal subdivision and skeleton
    perl::Object diag_fan = simplicial_with_diagonal(fan);
    perl::Object skeleton = skeleton_complex(diag_fan,fan_dim,true);
    
    //Extract values
    Matrix<Rational> old_rays = skeleton.give("RAYS");
    IncidenceMatrix<> skeleton_cones = skeleton.give("MAXIMAL_CONES");
    
    //Prepare solution matrix
    Matrix<Rational> result(skeleton_cones.rows(),0);
    
    //Iterate through all cones of the skeleton
    for(int tc = 0; tc < skeleton_cones.rows(); tc++) {
      //Create the function(s) Psi_tau
      Vector<int> tc_rays(skeleton_cones.row(tc));
      Matrix<Rational> psi_tau(0,old_rays.rows());
      for(int r = 0; r < tc_rays.dim(); r++) {
	psi_tau /= (unit_vector<Rational>(old_rays.rows(), tc_rays[r]));
      }
            
      //Compute the divisor
      perl::Object divisor = divisorByValueMatrix(diag_fan, psi_tau);
      
      //Extract cones, rays and weights
      Matrix<Rational> div_rays = divisor.give("RAYS");
      IncidenceMatrix<> div_cones = divisor.give("MAXIMAL_CONES");
      Vector<Integer> div_weights = divisor.give("TROPICAL_WEIGHTS");
      
      //Associate to each ray its original index 
      Map<int,int> div_ray_to_old;
      for(int dr = 0; dr < div_rays.rows(); dr++) {
	for(int oray = 0; oray < old_rays.rows(); oray++) {
	    if(old_rays.row(oray) == div_rays.row(dr)) {
	      div_ray_to_old[dr] = oray; break;
	    }
	}
      }//END translate ray indices
      
      //Now go through all d-dimensional cones in the divisor and insert their weight at the appropriate point
      Vector<Integer> tau_column(skeleton_cones.rows());
      for(int rho = 0; rho < div_cones.rows(); rho++) {
	//Map rho rays to old rays
	Set<int> rho_old = 
	    attach_operation(div_cones.row(rho), pm::operations::associative_access<Map<int,int>, int>(&div_ray_to_old));
	//Find the original cone equal to that
	for(int oc = 0; oc < skeleton_cones.rows(); oc++) {
	  if( (skeleton_cones.row(oc) * rho_old).size() == fan_dim) {
	    tau_column[oc] = div_weights[rho]; break;
	  }
	}
      }//END iterate divisor cones
      result |= tau_column;
      
    }//END iterate skeleton cones
    
    return result;
    
  }//END simplicial_piecewise_kernel
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Matrix<Rational> simplicial_diagonal_system(perl::Object fan) {
    //Extract values
    int fan_dim = fan.give("CMPLX_DIM");
    Matrix<Rational> fan_rays = fan.give("RAYS");
    
    //Compute diagonal subdivision and skeleton
    perl::Object diag_fan = simplicial_with_diagonal(fan);
    perl::Object skeleton = skeleton_complex(diag_fan,fan_dim,true);
    
    Matrix<Rational> nspace = simplicial_piecewise_system(fan);
    
    //Identify diagonal cones
    IncidenceMatrix<> sk_cones = skeleton.give("MAXIMAL_CONES");
    Set<int> diag_rays = sequence(2*fan_rays.rows(), fan_rays.rows());
    Set<int> diag_cones;
    for(int c = 0; c < sk_cones.rows(); c++) {
      if((diag_rays * sk_cones.row(c)).size() == sk_cones.row(c).size()) {
	diag_cones += c;
      }
    }
    
    Vector<Rational> v(nspace.rows());
    v.slice(diag_cones) = ones_vector<Rational>(diag_cones.size());
    
    nspace |= v;
    return nspace;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Intersection theory / Diagonal-related"
		    "# This function takes a simplicial fan F (in non-homog. coordinates and without "
		    "# lineality space) and computes the coarsest subdivision of F x F containing all "
		    "# diagonal rays (r,r)"
		    "# @param WeightedComplex F A simplicial fan without lineality space in non-homog."
		    "# coordinates"
		    "# @return WeightedComplex The product complex FxF subdivided such that it contains "
		    "# all diagonal rays",
		    &simplicial_with_diagonal, "simplicial_with_diagonal(WeightedComplex)");
  
  UserFunction4perl("# @category Intersection theory / Diagonal-related"
		    "# This function takes a d-dimensional simplicial fan F and computes the linear system "
		    "# defined in the following way: For each d-dimensional cone t in the diagonal subdivision of FxF, let psi_t be the "
		    "# piecewise polynomial defined by subsequently applying the rational functions that "
		    "# are 1 one exactly one ray of t and 0 elsewhere. Now for which coefficients a_t"
		    "# is sum_t a_t psi_t * (FxF) = 0?"
		    "# @param WeightedComplex F  A simplicial fan without lineality space in non-homog."
		    "# coordinates"
		    "# @return Matrix<Rational> The above mentioned linear system. The rows "
		    "# are equations, the columns correspond to d-dimensional cones of FxF in the order given "
		    "# by skeleton_complex(simplicial_with_diagonal(F), d,TRUE)",
		    &simplicial_piecewise_system, "simplicial_piecewise_system(WeightedComplex)");
  
  UserFunction4perl("# @category Intersection theory / Diagonal-related"
		    "# This function computes the inhomogeneous version of simplicial_piecewise_system"
		    "# in the sense that it computes the result of the above mentioned function (i.e. "
		    "# which coefficients for the piecewise polynomials yield the zero divisor)"
		    "# and adds another column at the end where only the entries corresponding to the "
		    "# diagonal cones are 1, the rest is zero. This can be seen as appending a non-zero"
		    "# solution to the system",
		    &simplicial_diagonal_system, "simplicial_diagonal_system(WeightedComplex)");
  
}}