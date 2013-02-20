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
 * Computes basic stuff on MinMaxFunctions
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/minmax_functions.h"

namespace polymake { namespace atint { 
  
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object add_minmax_functions(perl::Object f, perl::Object g, bool reduce) {
    //Extract relevant values
    Matrix<Rational> fmatrix = f.give("FUNCTION_MATRIX");
    bool f_uses_min = f.give("USES_MIN");
    Matrix<Rational> gmatrix = g.give("FUNCTION_MATRIX");
    bool g_uses_min = g.give("USES_MIN");
    
    if(f_uses_min != g_uses_min) {
      throw std::runtime_error("Cannot compute sum of tropical polynomials as polynomial. Different min/max-attitude");
    }
    
    if(fmatrix.cols() != gmatrix.cols()) {
      throw std::runtime_error("Cannot compute sum of tropical polynomials. Different domain.");
    }
    
    Matrix<Rational> rmatrix(0,fmatrix.cols());
    
    for(int i = 0; i < fmatrix.rows(); i++) {
      for(int j = 0; j < gmatrix.rows(); j++) {
	rmatrix /= (fmatrix.row(i) + gmatrix.row(j));
      }
    }
    
    if(reduce) {
      //Move constant coefficients to the front and add a homogenizing one
      Matrix<Rational> shifted_matrix = 
	rmatrix.col(rmatrix.cols()-1) | rmatrix.minor(All,sequence(0,rmatrix.cols()-1));
      shifted_matrix = ones_vector<Rational>(shifted_matrix.rows()) | shifted_matrix;
      
      //Canonicalize vertices
      solver<Rational> sv;
      std::pair<Matrix<Rational>, Matrix<Rational> > facets = 
	sv.enumerate_facets(shifted_matrix, Matrix<Rational>(0,shifted_matrix.cols()), false,false);
      rmatrix = sv.enumerate_vertices(facets.first,facets.second,false,true).first;
      
      //Put matrix back into proper form
      rmatrix = rmatrix.minor(All,sequence(2,rmatrix.cols()-2)) | rmatrix.col(1);      
      
      
    }
    
    perl::Object result("MinMaxFunction");
      result.take("FUNCTION_MATRIX") << rmatrix;
      result.take("USES_MIN") << f_uses_min;
      
    return result;
    
  }
  
  perl::Object test(int n, int j) {
    
    Array<Set<int> > jsubsets = pm::all_subsets_of_k(sequence(0,n),j);
    
    perl::Object result("MinMaxFunction");
    for(int s = 0; s < jsubsets.size(); s++) {
      perl::Object f("MinMaxFunction");
      Matrix<Rational> xmatrix(2*j, n);
      Matrix<Rational> ymatrix(2*j, n);
      xmatrix.minor(sequence(0,j),jsubsets[s]) = unit_matrix<Rational>(j);
      ymatrix.minor(sequence(j,j), jsubsets[s]) = unit_matrix<Rational>(j);
      Matrix<Rational> fmatrix = xmatrix | ymatrix | zero_vector<Rational>();
      f.take("FUNCTION_MATRIX") << fmatrix;
      f.take("USES_MIN") << false;
      if(s == 0) result = f;
      else result = add_minmax_functions(result,f);
    }
    
    return result;
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object scale_minmax_function(perl::Object f, Rational a) {
    //Extract values
    Matrix<Rational> fmatrix = f.give("FUNCTION_MATRIX");
    bool uses_min = f.give("USES_MIN");
    
    fmatrix = a * fmatrix;
    
    perl::Object result("MinMaxFunction");
      result.take("FUNCTION_MATRIX") << fmatrix;
      result.take("USES_MIN") << uses_min;
      
    return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Morphisms and functions"
		    "# Computes the sum of two MinMaxFunctions that are supposed to have the same USES_MIN property. You can also simply use '+' instead."
		    "# @param MinMaxFunction f An arbitrary MinMaxFunction"
		    "# @param MinMaxFunction g A MinMaxFunction such that f->[[USES_MIN]] == g->[[USES_MIN]] and is"
		    "# defined on the same domain."
		    "# @param Bool reduce Optional. False by default. If true, the function reduces the amount of terms by computing vertices of the newton polytope of the sum"
		    "# @return MinMaxFunction The sum of both functions. ",
		    &add_minmax_functions, "add_minmax_functions(MinMaxFunction, MinMaxFunction;$= 0)");
  
  UserFunction4perl("# @category Morphisms and functions"
		    "# Scales a MinMaxFunction by a given Rational a. You can also use the '*' operator."
		    "# @param MinMaxFunction f An arbitrary MinMaxFunction"
		    "# @param Rational a A scalar values"
		    "# @return MinMaxFunction The scaled function",
		    &scale_minmax_function, "scale_minmax_function(MinMaxFunction, Rational)");
  
  Function4perl(&test,"mmtest($,$)");
  
}}