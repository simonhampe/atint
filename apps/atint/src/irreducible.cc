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
 * This file contains functions to check irreducibility of a tropical cycle
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"
#include "polymake/atint/normalvector.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief This function computes a system of linear equations that determine whether a given tropical cycle is irreducible.
   @param perl::Object C A WeightedComplex object that is assumed to be balanced
   @param Map<int,int> A reference that will be set to a map describing how to translate an element in the kernel of this matrix. More precisely: For maximal cone i map(i) describes the index of the solution (corresponding to a column index of the result of this function) where you should look for the weight of the cone. I.e. if you assign to all cones the weight that you find in an integer solution at the positions described by this map, you obtain a balanced cycle.
   @return Matrix<Integer> A system of linear equations: For each codimension one cell t it contains a block describing the balancing condition around this cell (i.e. of the form (u_1,... u_k,v_1,..v_l), where the u_i are the lattice normal vectors at t and v_1,...v_l is a lattice basis for t. Additionally the matrix contains equations to make sure that a single maximal cone is assigned the same weight in each of the above systems
   */
  Matrix<Integer> cycle_irreducibility_system(perl::Object C, Map<int,int> &index_translator) {
    
    //dbgtrace << "Extracting values..." << endl;
    
    //Extract values
    int ambient_dim = C.give("CMPLX_AMBIENT_DIM");
    int dim = C.give("CMPLX_DIM");
    Matrix<Rational> rays = C.give("RAYS");
    Matrix<Rational> lineality = C.give("LINEALITY_SPACE");
    IncidenceMatrix<> maximal_cones = C.give("MAXIMAL_CONES");
    IncidenceMatrix<> codim_1_faces = C.give("CODIM_1_FACES");
    IncidenceMatrix<> codim_in_maximal = C.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximal_to_codim = T(codim_in_maximal);
    Map<int,Map<int, Vector<Integer> > > lattice_normals = C.give("LATTICE_NORMALS");
    bool uses_homog = C.give("USES_HOMOGENEOUS_C");
    
    //Compute future dimension of (the upper part of ) the matrix
    int f_rows = codim_1_faces.rows() * ambient_dim;
    int f_cols = codim_1_faces.rows() * (dim-1);
    for(int s = 0; s < codim_1_faces.rows(); s++) {
      f_cols += codim_in_maximal.row(s).size();
    }
    
    Matrix<Integer> F(f_rows,f_cols);
    
    Vector<bool> cone_translated(maximal_cones.rows()); //Whether the column index for cone i is computed
    index_translator = Map<int,int>();
    
    //Generate balancing equation blocks
    int rowindex = 0;
    int colindex = 0;
    Map<int,Map<int,int> > normal_col_index; //Saves for each normal vector its column index
    for(int tau = 0; tau < codim_1_faces.rows(); tau++) {
      //dbgtrace << "Computing equation block for codim one cell " << tau << endl;
      normal_col_index[tau] = Map<int,int>();
      
      Matrix<Integer> M_tau(ambient_dim,0);
      
      //dbgtrace << "Appending normal vectors " << endl;
      
      //Add all normal vectors
      Set<int> adjacent_maximal = codim_in_maximal.row(tau);
      for(Entire<Set<int> >::iterator sigma = entire(adjacent_maximal); !sigma.at_end(); sigma++) {
	  Vector<Integer> lnormal = ((lattice_normals[tau])[*sigma]);
	    if(uses_homog) lnormal = lnormal.slice(~scalar2set(0));
	  M_tau |= lnormal;
	  (normal_col_index[tau])[*sigma] = colindex + M_tau.cols()-1;
	  if(!cone_translated[*sigma]) {
	      cone_translated[*sigma] = true;
	      index_translator[*sigma] = (normal_col_index[tau])[*sigma];
	  }
      }
      
      //dbgtrace << "Appending lattice basis" << endl;
      
      //Compute a lattice basis for this cell
      if(dim > 1) {
	Matrix<Integer> lbasis = latticeBasisFromRays(rays.minor(codim_1_faces.row(tau),All),lineality,uses_homog);
	  if(uses_homog) lbasis = lbasis.minor(All,~scalar2set(0));
	M_tau |= T(lbasis);
      }
      
      //dbgtrace << "Inserting equation block " <<  M_tau << endl;
      
      //Insert M_tau
      F.minor(sequence(rowindex,M_tau.rows()), sequence(colindex, M_tau.cols())) = M_tau;
      rowindex += M_tau.rows();
      colindex += M_tau.cols();
      
    }//END generate balancing equation blocks
    
    
    //Now generate maximal cone glue equations
    for(int sigma = 0; sigma < maximal_cones.rows(); sigma++) {
	//dbgtrace << "Computing glueing equations for maximal cone " << sigma << endl;
	Vector<int> faces(maximal_to_codim.row(sigma));
	if(faces.dim() <= 1) continue;
	//For all pairs of codim one cells in sigma, make sure the corresponding normal vectors
	//are assigned the same weight, i.e. the corresponding entries in the solution must be the same
	//More precisely, we just have to demand that the solution for the system of the "first" codimension
	//one cell of sigma is the same as the solution for any other system
	Matrix<Integer> N_sigma(0,f_cols);
	
	for(int tau = 1; tau < faces.dim(); tau++) {
	  N_sigma /= ( unit_vector<Integer>(f_cols, (normal_col_index[faces[0]])[sigma]) - 
			 unit_vector<Integer>(f_cols, (normal_col_index[faces[tau]])[sigma]));
	}
	
	
	F /= N_sigma;
    }//END generate maximal cone glue equations
    
    return F;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Matrix<Integer> cycle_weight_space(perl::Object C) {
    
    IncidenceMatrix<> maximal_cones = C.give("MAXIMAL_CONES");
    
    //First we compute the irreducibility system
    Map<int,int> index_translator;
    Matrix<Integer> irred_system = cycle_irreducibility_system(C,index_translator);
    
    //Then compute the solution system
    Matrix<Integer> ker = null_space(irred_system);
    
    Matrix<Rational> solutions(ker.rows(), maximal_cones.rows());
    //Go through all kernel solutions
    for(int k = 0; k < ker.rows(); k++) {
      //Translate into weight solution
      for(int m = 0; m< maximal_cones.rows(); m++) {
	solutions(k,m) = ker(k, index_translator[m]);
      }
    }
    
    //Now compute a lattice basis
    return latticeBasisFromRays(solutions, Matrix<Rational>(0,solutions.cols()), false);
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  bool is_irreducible(perl::Object C) {
    //First we compute the gcd of the weights of C
    Vector<Integer> weights = C.give("TROPICAL_WEIGHTS");
    if(weights.dim() == 0) return true;
    Integer g = weights[0];
    for(int w = 1; w < weights.dim(); w++) {
      g = gcd(g,weights[w]);
    }
    if(g != 1) {
      return false;
    }
    
    Map<int,int> dummy;
    Matrix<Integer> F = cycle_irreducibility_system(C, dummy);
    return rank(null_space(F)) == 1;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
//   Function4perl(&cycle_irreducibility_system, "cycle_irreducibility_system(WeightedComplex)");
  
  UserFunction4perl("# @category Tropical geometry / Irreducibility "
		    "# This function computes whether a given tropical cycle is irreducible"
		    "# @param WeightedComplex C A tropical cycle"
		    "# @return Bool Whether the cycle is irreducible",
		    &is_irreducible, "is_irreducible(WeightedComplex)");
  
  UserFunction4perl("# @category Tropical geometry / Irreducibility"
		    "# This function computes a Z-basis for the space of weight distributions on a "
		    "# tropical cycle making it balanced (i.e. the cycle is irreducible, if and only if"
		    "# the basis consists of one element"
		    "# @param WeightedComplex C A tropical cycle"
		    "# @return Matrix<Integer> A Z-basis of the weight space (given as row vectors)",
		    &cycle_weight_space, "cycle_weight_space(WeightedComplex)");
	
  
}}