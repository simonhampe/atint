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
#include "polymake/atint/basicoperations.h"
#include "polymake/atint/cdd_helper_functions.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  perl::Object weighted_curve_space(int f, int t) {
    
    //First create product of moduli spaces
    dbgtrace << "Creating moduli space" << endl;
    std::vector<perl::Object> moduli;
    perl::Object mf = CallPolymakeFunction("tropical_m0n",f);
      int mf_amb = f*(f-3)/2;
    perl::Object mfp = CallPolymakeFunction("tropical_m0n",f+1);
      int mfp_amb = (f+1)*(f-2)/2;
    moduli.push_back(mf);
    for(int i = 1; i <= t; i++) { moduli.push_back(mfp);}
    
    dbgtrace << "Computing product" << endl;
    
    perl::Object prod = compute_product_complex(moduli);
    
    //Compute the product of the forgetful map
    
    dbgtrace << "Creating maps " << endl;
    
    perl::Object ffmap = CallPolymakeFunction("forgetful_map",f+1,Set<int>(f+1));
    Matrix<Rational> ffmatrix = ffmap.give("MATRIX");
    
    dbgtrace << "Gluing matrix together" << endl;
    
    Matrix<Rational> fibre_matrix(t*mf_amb, mf_amb + t*mfp_amb);
    //Insert unit matrices for projection onto first coordinate
    dbgtrace << "inserting unit matrix" << endl;
    for(int i = 0; i < t; i++) {
      fibre_matrix.minor(sequence(i*mf_amb, mf_amb),sequence(0,mf_amb)) = - unit_matrix<Rational>(mf_amb);
    }//END insert unit matrices
    
    dbgtrace << "Inserting forgetful map" << endl;
    
    //Insert forgetful maps
    for(int i = 0; i < t; i++) {
      fibre_matrix.minor(sequence(i*mf_amb,mf_amb), sequence(mf_amb + i*mfp_amb, mfp_amb)) =
			      ffmatrix;
    }//END insert forgetful maps
    
    Matrix<Rational> fibre_space = null_space(fibre_matrix);
    
    //Now intersect each cone with the equation that ft_{f+1} on the last coordinates equals the first 
    //coordinate, i.e. fibre_matrix * vector = 0
    
    dbgtrace << "Computing intersections" << endl;
    
    Matrix<Rational> prays = prod.give("RAYS");
    IncidenceMatrix<> pcones = prod.give("MAXIMAL_CONES");
    
    Matrix<Rational> rrays(0, prays.cols());
    Vector<Set<int> > rcones;
    
    for(int pc = 0; pc < pcones.rows(); pc++) {
      
      Matrix<Rational> inter = cdd_cone_intersection(prays.minor(pcones.row(pc),All),
						      Matrix<Rational>(0,prays.cols()),
						      Matrix<Rational>(0,prays.cols()),
						      fibre_space,false).first;
      if(rank(inter) == f+t-3) {
	rcones |= Set<int>(insert_rays(rrays,inter,false,false));
      }      
      
    }//END iterate product cones
    
    //Return result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rrays;
      result.take("MAXIMAL_CONES") << rcones;
      result.take("TROPICAL_WEIGHTS") << ones_vector<Integer>(rcones.dim());
    return result;
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&weighted_curve_space,"wcs($,$)");
  
}}