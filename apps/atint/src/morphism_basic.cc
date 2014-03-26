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
#include "polymake/atint/refine.h"
#include "polymake/atint/morphism_composition.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the sum of two morphisms (which should be defined on the same support)
   @param perl::Object f A Morphism object
   @param perl::Object g A Morphism object, whose DOMAIN has the same support as f's DOMAIN and whose image has the same ambient dimension as f's image
   @return perl::Object A Morphism object representing f+g
   */
  perl::Object add_morphisms(perl::Object f, perl::Object g) {
    //First we treat the special case where both are global
    bool f_global = f.give("IS_GLOBAL");
    bool g_global = g.give("IS_GLOBAL");
    if(f_global && g_global) {
      //dbgtrace << "Both are global - computing global sum" << endl;
      Matrix<Rational> fmatrix = f.give("MATRIX");
      Vector<Rational> ftranslate = f.give("TRANSLATE");
      Matrix<Rational> gmatrix = g.give("MATRIX");
      Vector<Rational> gtranslate = g.give("TRANSLATE");
      
      perl::Object addition("Morphism");
	addition.take("MATRIX") << fmatrix + gmatrix;
	addition.take("TRANSLATE") << ftranslate + gtranslate;
	
      return addition;
    }
    
    //dbgtrace << "Homogenizing where necessary" << endl;
    
    //First we homogenize where necessary
    perl::Object fDomain = f.give("DOMAIN");
    perl::Object gDomain = g.give("DOMAIN");
    bool fhomog = fDomain.give("USES_HOMOGENEOUS_C");
    bool ghomog = gDomain.give("USES_HOMOGENEOUS_C");
    if(fhomog || ghomog) {
      if(!fhomog) {
	f = f.CallPolymakeMethod("homogenize");
	fDomain = f.give("DOMAIN");
      }
      if(!ghomog) {
	g = g.CallPolymakeMethod("homogenize");
	gDomain = g.give("DOMAIN");
      }
    }
    
    //dbgtrace << "Computing refinement " << endl;
    
    //Then compute the common refinement of the domains
    RefinementResult r = refinement(fDomain,gDomain,true,true,false,true);
      perl::Object nDomain = r.complex;
	Matrix<Rational> x_rayrep = r.rayRepFromX;
	Matrix<Rational> y_rayrep = r.rayRepFromY;
	Matrix<Rational> x_linrep = r.linRepFromX;
	Matrix<Rational> y_linrep = r.linRepFromY;
	
	Matrix<Rational> f_rayval = f.give("RAY_VALUES");
	Matrix<Rational> g_rayval = g.give("RAY_VALUES");
	Matrix<Rational> f_linval = f.give("LIN_VALUES");
	  f_linval = T(f_linval);
	Matrix<Rational> g_linval = g.give("LIN_VALUES");
	  g_linval = T(g_linval);
	
	Matrix<Rational> fval = T(f_rayval) | f_linval;  
	Matrix<Rational> gval = T(g_rayval) | g_linval;
	
      Matrix<Rational> rays = nDomain.give("CMPLX_RAYS");
      Matrix<Rational> linspace = nDomain.give("LINEALITY_SPACE");
      
      //dbgtrace << "Computing sums of values " << endl;
      
      //dbgtrace << "Ray values " << endl;
      //Now compute ray values
      Matrix<Rational> rValues(0, f_rayval.rows());
      for(int r = 0; r < rays.rows(); r++) {
	rValues /= (fval * x_rayrep.row(r)) + (gval * y_rayrep.row(r) );
      }
      
      //dbgtrace << "Lin values " << endl;
      //Now compute lin values
      Matrix<Rational> lValues(0, f_linval.rows());
      for(int l = 0; l < linspace.rows(); l++) {
	lValues /= (f_linval * x_linrep.row(l)) + (g_linval * y_linrep.row(l) );
      }
      
      perl::Object result("Morphism");
	result.take("DOMAIN") << nDomain;
	result.take("RAY_VALUES") << rValues;
	result.take("LIN_VALUES") << lValues;
      
      return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&morphism_homogenize,"morphism_homogenize(Morphism)");
  Function4perl(&add_morphisms,"add_morphisms(Morphism, Morphism)");
  
  
}}