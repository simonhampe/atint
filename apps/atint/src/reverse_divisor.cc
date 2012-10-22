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
 * 
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Matrix<Rational> cutting_functions(perl::Object fan, Vector<Integer> weight_aim) {
    
    //Extract values
    Map<int, Map<int, Vector<Rational> > > summap = fan.give("LATTICE_NORMAL_FCT_VECTOR");
    Matrix<Rational> summatrix = fan.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
    IncidenceMatrix<> codim_in_cones = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    Vector<Integer> weights = fan.give("TROPICAL_WEIGHTS");
    
    //Compute equation matrix
    Matrix<Rational> equations(0,summatrix.cols());
    
    //Compute equation for each codimension one cone
    for(int c = 0; c < codim_in_cones.rows(); c++) {
      //Coefficients are the sum of the lattice normal function vectors minus
      //the sum function vector
      Vector<Rational> ceq(summatrix.cols());
      Set<int> adjacent_cones = codim_in_cones.row(c);
      for(Entire<Set<int> >::iterator ac = entire(adjacent_cones); !ac.at_end(); ac++) {
	ceq += weights[*ac] * (summap[c])[*ac];
      }//END iterate adjacent maximal cells
      ceq -= summatrix.row(c);
      equations /= ceq;
    }//END iterate codim 1 cells
    
    //Finally add desired weights as additional coefficients
    equations |= weight_aim;
//     dbgtrace << "Equations: " << equations << endl;
    return null_space(equations);
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Rational Functions / Inverse problem"
		    "# Takes a weighted complex and a list of desired weights on its codimension one"
		    "# faces and computes all possible rational functions on (this subdivision of )"
		    "# the complex"
		    "# @param WeightedComplex F A tropical variety"
		    "# @param Vector<Integer> weight_aim A list of weights, whose length should be equal to the number of CODIM_1_FACES. Gives the desired weight on each "
		    "# codimension one face"
		    "# @return Matrix<Rational> The space of rational functions defined on this "
		    "# particular subdivision. Each row is a generator. The columns correspond to "
		    "# values on CMPLX_RAYS and LINEALITY_SPACE, except the last one, which is either 0 (then this "
		    "# function cuts out zero and can be added to any solution) or non-zero (then "
		    "# normalizing this entry to -1 gives a function cutting out the desired weights "
		    "# on the codimension one skeleton",
		    &cutting_functions, "cutting_functions(WeightedComplex, Vector<Integer>)");
  
}}