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
 * Contains basic operations on morphisms
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a morphism object and homogenizes it
   @param perl::Object morphism A Morphism object
   @return perl::Object A Morphism object, the same morphism as given as input, but on a homogenized domain
   */
  perl::Object morphism_homogenize(perl::Object morphism) {
    //Extract values
    perl::Object domain = morphism.give("DOMAIN");
    bool is_homog = domain.give("USES_HOMOGENEOUS_C");
    if(is_homog) return morphism;
    Matrix<Rational> ray_values = morphism.give("RAY_VALUES");
    Matrix<Rational> lin_values = morphism.give("LIN_VALUES");
    
    //Homogenize
    perl::Object homog_domain = domain.CallPolymakeMethod("homogenize");
    Matrix<Rational> cmplx_rays = homog_domain.give("CMPLX_RAYS");
    
    //Find the new ray with x0 = 1
    for(int r = 0; r < cmplx_rays.rows(); r++) {
      if(cmplx_rays(r,0) == 1) {
	Matrix<Rational> new_ray_values = ray_values.minor(sequence(0,r),All) / zero_vector<Rational>();
	new_ray_values /= ray_values.minor(sequence(r,ray_values.rows()-r),All);
	perl::Object result("Morphism");
	  result.take("DOMAIN") << homog_domain;
	  result.take("RAY_VALUES") << new_ray_values;
	  result.take("LIN_VALUES") << lin_values;
	return result;
      }
    }
    
    //To avoid compiler warnings:
    return perl::Object("Morphism");
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&morphism_homogenize,"morphism_homogenize(Morphism)");
  
}}