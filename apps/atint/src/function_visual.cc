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
 * This file computes visualization data for rational functions and morphisms
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/morphism_pullback.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the labels for each cell of a function domain, given its ray and lineality values
   @param perl::Object domain A WeightedComplex object, representing the domain of the function
   @param Matrix<Rationl> ray_values Values of the function on the rays, given as row vectors in non-homog. coordinates
   @param Matrix<Rational> lin_values Values of the function on the lineality space, given as row vectors in non-homog. coordinates.
   @return A list of std::strings
   */
  perl::ListReturn computeFunctionLabels(perl::Object domain, Matrix<Rational> ray_values, Matrix<Rational> lin_values) {
    //Extract values
    Matrix<Rational> rays = domain.give("CMPLX_RAYS");
    IncidenceMatrix<> cones = domain.give("CMPLX_MAXIMAL_CONES");
    Matrix<Rational> lineality = domain.give("LINEALITY_SPACE");
    bool uses_homog = domain.give("USES_HOMOGENEOUS_C");
    
    perl::ListReturn result;
    
    
    
    for(int mc = 0; mc < cones.rows(); mc++) {
      dbgtrace << "Computing representation of cone " << mc << endl;
      Matrix<Rational> matrix;
      Vector<Rational> translate;
      computeConeFunction(rays.minor(cones.row(mc),All), lineality, uses_homog, ray_values.minor(cones.row(mc),All), lin_values, translate, matrix);
      
      std::ostringstream sstream;
      pm::PlainPrinter<> rep(sstream);
      if(matrix.rows() > 1) {
	rep << "(" << translate << ")" << " + ";
	for(int i = 0; i < matrix.rows(); i++) {
	  rep << "[" << matrix.row(i) << "]";
	}
      }
      //We have a special representation format for functions to R
      else {
	bool hadnonzeroterm = false;
	for(int i = 0; i < matrix.cols(); i++) {
	    if(matrix(0,i) != 0) {
	      if(hadnonzeroterm) rep << " + ";
	      hadnonzeroterm = true;
	      if(matrix(0,i) < 0) rep << "(" << matrix(0,i) << ")";
	      else rep << matrix(0,i);
	      rep << "*x_" << (i+1);
	    }
	}
	if(translate[0] < 0 && hadnonzeroterm) rep << " - " << (-translate[0]);
	if(translate[0] > 0 && hadnonzeroterm) rep << " + " << translate[0];
	if(!hadnonzeroterm) rep << translate[0];
      }
      result << sstream.str();
      
    }
    
    return result;
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&computeFunctionLabels,"computeFunctionLabels(WeightedComplex, Matrix<Rational>, Matrix<Rational>)");
  
}}