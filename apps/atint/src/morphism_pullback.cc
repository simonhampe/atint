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
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 This file contains functions to compute pull-backs of rational functions
 along morphisms
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/linalg.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/morphism_pullback.h"
#include "polymake/atint/morphism_composition.h"

namespace polymake { namespace atint { 
    
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  void computeConeFunction(const Matrix<Rational> &rays, const Matrix<Rational> &linspace, bool uses_homog, const Matrix<Rational> &ray_values, const Matrix<Rational> &lin_values, Vector<Rational> &translate, Matrix<Rational> &matrix) {
    
    //First we need to compute a cone basis
    
    //Convert vertices to differences of vertices in homogeneous case
    Vector<Rational> basepoint = zero_vector<Rational>(rays.cols());
    Vector<Rational> basepoint_value = zero_vector<Rational>(ray_values.cols());
    Matrix<Rational> converted_rays;
    Matrix<Rational> converted_values;
    if(!uses_homog) {
      converted_rays = rays;
      converted_values = ray_values;
    }
    else {
      bool basepoint_found = false;
      for(int r = 0; r < rays.rows(); r++) {
	if(rays(r,0) == 1) {
	    if(basepoint_found) {
	      converted_rays /= (rays.row(r) - basepoint);
	      converted_values /= (ray_values.row(r) - basepoint_value);
	    }
	    else {
	      basepoint = rays.row(r);
	      basepoint_value = ray_values.row(r);
	      basepoint_found = true;
	    }
	}
	else {
	    converted_rays /= rays.row(r);
	    converted_values /= ray_values.row(r);
	}
      }
    }//END convert rays
    converted_rays /= linspace;
    converted_values /= lin_values;
    //Remove homog. coordinate if necessary
    if(uses_homog) {
      converted_rays = converted_rays.minor(All,~scalar2set(0));
      basepoint = basepoint.slice(~scalar2set(0));
    }
    
    //dbgtrace << "Basepoint: " << basepoint << endl;
    //dbgtrace << "Basepoint value: " << basepoint_value << endl;
    
    //Compute basis of rays
    Set<int> ray_basis = basis_rows(converted_rays);
    converted_rays = converted_rays.minor(ray_basis,All);
    
    //dbgtrace << "Converted rays: " << converted_rays << endl;
    
    //Now compute a column basis for the computation of the transformation matrix
    Set<int> I = basis_cols(converted_rays);
    Matrix<Rational> inverse = inv(T( converted_rays.minor(All,I)));
    Matrix<Rational> trafo(converted_rays.rows(), converted_rays.cols());
    trafo.minor(All,I) = inverse;
    
    //dbgtrace << "Trafo: " << trafo << endl;
    
    //Compute function matrix:
    Matrix<Rational> values = converted_values.minor(ray_basis,All);
    matrix = T(values) * trafo;
    
    //dbgtrace << "Matrix: " << matrix << endl;
    
    //Finally, compute the translate
    translate = basepoint_value - matrix * basepoint;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  void computeConeFunction(const Matrix<Rational> &rays, const Matrix<Rational> &linspace, bool uses_homog, const Vector<Rational> &ray_values, const Vector<Rational> &lin_values, Rational &translate, Vector<Rational> &functional) {
    //Convert input values
    Matrix<Rational> convert_ray_values(0,ray_values.dim());
      convert_ray_values /= ray_values;
    Matrix<Rational> convert_lin_values(0,lin_values.dim());
      convert_lin_values /= lin_values;
    Vector<Rational> convert_translate;
    Matrix<Rational> convert_functional;
    
    //Compute result
    computeConeFunction(rays, linspace, uses_homog, convert_ray_values, convert_lin_values, convert_translate, convert_functional);
    
    //Convert result
    translate = convert_translate[0];
    functional = convert_functional.row(0);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  perl::Object pb_minmax_global(perl::Object morphism, perl::Object function) {
    //Extract values
    Matrix<Rational> m_matrix = morphism.give("MATRIX");
    Vector<Rational> m_translate = morphism.give("TRANSLATE");
    Matrix<Rational> f_linear = function.give("LINEAR_COEFFICIENTS");
    Vector<Rational> f_constant = function.give("CONSTANT_COEFFICIENTS");
    bool uses_min = function.give("USES_MIN");
    
    Vector<Rational> new_constant = f_linear*m_translate + f_constant;    
    Matrix<Rational> new_linear = f_linear*m_matrix;
    
    perl::Object result("MinMaxFunction");
      result.take("LINEAR_COEFFICIENTS") << new_linear;
      result.take("CONSTANT_COEFFICIENTS") << new_constant;
      result.take("USES_MIN") << uses_min;
      
    return result;    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the pull-back of a RationalFunction along a Morphism. The function assumes that the morphism is either global and surjective or that its image is contained in the domain of the function.
   @param perl::Object morphism A Morphism object
   @param perl::Object function A RationalFunction object, If the morphism is not global and surjective, its domain should contain the image of the morphism.
   @return perl::Object A RationalFunction object, the pull-back of function along morphism (always given on a homogeneous domain)
   */
  perl::Object pb_general(perl::Object morphism, perl::Object function) {
    //First we need to convert the function to a "morphism"
    perl::Object fdomain = function.give("DOMAIN");
    Vector<Rational> frayvalues = function.give("RAY_VALUES");
    Vector<Rational> flinvalues = function.give("LIN_VALUES");
    bool uses_min = function.give("USES_MIN");
    
    Matrix<Rational> frayconv(frayvalues.dim(), 0);
      frayconv |= frayvalues;
    Matrix<Rational> flinconv(flinvalues.dim(), 0);
      flinconv |= flinvalues;
      
    perl::Object fmorphism("Morphism");
      fmorphism.take("DOMAIN") << fdomain;
      fmorphism.take("RAY_VALUES") << frayconv;
      fmorphism.take("LIN_VALUES") << flinconv;
    
    //Now compute composition (as a "morphism")
    
    perl::Object comp = morphism_composition(morphism, fmorphism);
    
    //Finally convert back to rational function
    
    perl::Object cdomain = comp.give("DOMAIN");
    Matrix<Rational> crayval = comp.give("RAY_VALUES");
    Matrix<Rational> clinval = comp.give("LIN_VALUES");
    
    perl::Object result("RationalFunction");
      result.take("DOMAIN") << cdomain;
      if(crayval.cols() > 0) result.take("RAY_VALUES") << crayval.col(0);
      if(clinval.cols() > 0) result.take("LIN_VALUES") << clinval.col(0);
      result.take("USES_MIN") << uses_min;
    
    return result;
    
  }//END pb_general
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&pb_minmax_global,"pb_minmax_global(Morphism,MinMaxFunction)");
  Function4perl(&pb_general, "pb_general(Morphism, RationalFunction)");
  
}}