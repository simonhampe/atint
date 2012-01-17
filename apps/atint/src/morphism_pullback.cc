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
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/refine.h"

namespace polymake { namespace atint { 
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Computes the pull-back of a min/max function along a global affine linear function
    @param perl::Object morphism A Morphism object that represents a global affine linear function
    @param perl::Object function A MinMaxFunction object
    @return perl::Object The pull-back function
  */
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
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&pb_minmax_global,"pb_minmax_global(Morphism,MinMaxFunction)");
  
}}