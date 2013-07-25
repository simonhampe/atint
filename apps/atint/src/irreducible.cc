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
  
//   //Documentation see perl wrapper
//   Matrix<Integer> cycle_weight_system(perl::Object C) {
//     
//     //Extract values
//     int ambient_dim = C.give("CMPLX_AMBIENT_DIM");
//     int dim = C.give("CMPLX_DIM");
//     Matrix<Rational> rays = C.give("RAYS");
//     Matrix<Rational> lineality = C.give("LINEALITY_SPACE");
//     IncidenceMatrix<> maximal_cones = C.give("MAXIMAL_CONES");
//     IncidenceMatrix<> codim_1_faces = C.give("CODIM_1_FACES");
//     IncidenceMatrix<> codim_in_maximal = C.give("CODIM_1_IN_MAXIMAL_CONES");
//     IncidenceMatrix<> maximal_to_codim = T(codim_in_maximal);
//     Map<int,Map<int, Vector<Integer> > > lattice_normals = C.give("LATTICE_NORMALS");
//     bool uses_homog = C.give("USES_HOMOGENEOUS_C");
//     
//     //Prepare matrix
//     Matrix<Rational> result(0, maximal_cones.rows());
//     
//     //Iterate all codimension one cells
//     for(int tau = 0; tau < codim_1_faces.rows(); tau++) {
//       //Create the linear equation system for the balancing condition at tau
//       Matrix<Rational> M_tau(ambient_dim,0);
//       
//       //dbgtrace << "Appending normal vectors " << endl;
//       
//       //Add all normal vectors
//       Set<int> adjacent_maximal = codim_in_maximal.row(tau);
//       for(Entire<Set<int> >::iterator sigma = entire(adjacent_maximal); !sigma.at_end(); sigma++) {
// 	  Vector<Integer> lnormal = ((lattice_normals[tau])[*sigma]);
// 	    if(uses_homog) lnormal = lnormal.slice(~scalar2set(0));
// 	  M_tau |= lnormal;
//       }//END add normal vectors
//       
//       //Compute a lattice basis for this cell
//       if(dim > 1) {
// 	Matrix<Integer> lbasis = latticeBasisFromRays(rays.minor(codim_1_faces.row(tau),All),lineality,uses_homog);
// 	  if(uses_homog) lbasis = lbasis.minor(All,~scalar2set(0));
// 	M_tau |= T(lbasis);
//       }//END compute lattice basis
//       
//       //Now compute the kernel of this system
//       Matrix<Rational> K_tau = null_space(M_tau);
//       
//       //dbgtrace << "Kernel of matrix is " << K_tau << endl;
//       
//       //Compute the projection to the weight space:
//       Matrix<Rational> P_tau(K_tau.rows(),maximal_cones.rows());
//       P_tau.minor(All,adjacent_maximal) = K_tau.minor(All, sequence(0, maximal_cones.rows()));
//       //The remaining cones can have arbitrary weights, so we append the appropriate unit vectors
//       P_tau /= (unit_matrix<Rational>(maximal_cones.rows()).minor(sequence(0,maximal_cones.rows()) - adjacent_maximal,All));
//       
//       
//       //dbgtrace << "Projected to " << P_tau << endl;
//       
//       //Now we compute the equations in the weight space,add them to the current equations
//       //and reduce
//       result /= null_space(P_tau);
//       result = result.minor(basis_rows(result),All);
//       
//     }//END iterate codimension one cells
//     
//     return makeInteger(result);
//     
//   }//END old_cycle_weight_system
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Matrix<Rational> cycle_weight_space(perl::Object C) {
    
    //Extract values
    int ambient_dim = C.give("CMPLX_AMBIENT_DIM");
    int dim = C.give("CMPLX_DIM");
    Matrix<Rational> rays = C.give("RAYS");
    Matrix<Rational> linspace = C.give("LINEALITY_SPACE");
    IncidenceMatrix<> maximal_cones = C.give("MAXIMAL_CONES");
    IncidenceMatrix<> codim_1_faces = C.give("CODIM_1_FACES");
    IncidenceMatrix<> codim_in_maximal = C.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximal_to_codim = T(codim_in_maximal);
    Map<int,Map<int, Vector<Integer> > > lattice_normals = C.give("LATTICE_NORMALS");
    bool uses_homog = C.give("USES_HOMOGENEOUS_C");
        
    if(dim == 0) {
      return unit_matrix<Rational>(maximal_cones.rows());
    }
    
    //dbgtrace << "Computing equivalence classes" << endl;
    
    //Compute equivalence classes of maximal cones
    //and while doing so identify those classes which can not be balanced (i.e. which have to have
    // weight 0).
    Vector<Set<int> > subdivision_classes;
    Set<int> zero_classes; //List of classes that must have weight 0.
    Set<int> hasBeenAdded; //List of cones that know their class
    Map<int,int> class_index; //Maps each cone to the index of its class in subdivision_classes
    for(int mc = 0; mc < maximal_cones.rows(); mc++) {
      if(!hasBeenAdded.contains(mc)) {
	Set<int> mc_class; mc_class += mc;
	class_index[mc] = subdivision_classes.dim();
	bool is_zero_class = false;
	//Elements in this queue are already in the class but their neighbours might not:
	std::list<int> queue;
	  queue.push_back(mc);
	while(queue.size() > 0) {
	  int node = queue.front(); queue.pop_front();
	  Set<int> node_codim = maximal_to_codim.row(node);
	  for(Entire<Set<int> >::iterator nc = entire(node_codim); !nc.at_end(); nc++) {
	    if(codim_in_maximal.row(*nc).size() == 2) {
	      int other_cone = *((codim_in_maximal.row(*nc) - node).begin());
	      if(!hasBeenAdded.contains(other_cone)) {
		hasBeenAdded += other_cone;
		mc_class += other_cone;
		class_index[other_cone] = subdivision_classes.dim();
		//If it is not a zero class yet, check if it is now
		if(!is_zero_class) {
		  //The lattice normals must add up to 0 mod V_{nc}
		  Vector<Rational> lsum ((lattice_normals[*nc])[node] + (lattice_normals[*nc])[other_cone]);
		  Matrix<Rational> vtau = rays.minor(codim_1_faces.row(*nc),All) / linspace;
		  if(rank(vtau/lsum) > rank(vtau)) is_zero_class = true;
		}//END check if zero class
	      }
	    }//END if two-valent
	  }//END iterate all codim-1-faces
	}//END iterate queue
	if(is_zero_class) zero_classes += subdivision_classes.dim();
	subdivision_classes |= mc_class;	
      }//END if not in a class yet
    }//END iterate maximal cones to find equiv. classes
    
    //dbgtrace << "Zero classes: " << zero_classes << endl;
    //dbgtrace << "Clas indices: " << class_index << endl;
    
    //Linear system matrix
    Matrix<Rational> system_matrix(0,subdivision_classes.dim());
    
    //dbgtrace << "Computing system matrices" << endl;
    
    //Iterate codimension one faces to compute local signature systems
    for(int tau = 0; tau < codim_1_faces.rows(); tau++) {
      
      //dbgtrace << "Computing signature neighbours" << endl;
      //Create local system matrix for tau
      Matrix<Rational> Mtau(ambient_dim,0);
      //Find signature neighbours
      Set<int> nonsig_neighbours;
      Set<int> all_neighbours = codim_in_maximal.row(tau);
      for(Entire<Set<int> >::iterator mc = entire(all_neighbours); !mc.at_end(); mc++) {
	//See if there is another one with the same class
	for(Entire<Set<int> >::iterator other_mc = entire(all_neighbours); !other_mc.at_end(); other_mc++) {
	    if(*other_mc != *mc && class_index[*mc] == class_index[*other_mc]) {
	      nonsig_neighbours += *mc; break;
	    }
	}//END iterate rest of neighbours
      }//END iterate all neighbours to find signature ones
      Vector<int> sig_neighbours(all_neighbours - nonsig_neighbours);
      
      //dbgtrace << "Signature neighbours " << sig_neighbours << endl;
      
      for(int smc = 0; smc < sig_neighbours.dim(); smc++) {
	Mtau |= (lattice_normals[tau])[sig_neighbours[smc]];
      }//END iterate signature neighbours
      
      //dbgtrace << "Computing lattice basis" << endl;
      
      if(dim > 1) {
	Matrix<Integer> lbasis = 	
	  latticeBasisFromRays(rays.minor(codim_1_faces.row(tau),All),linspace,uses_homog);
	Mtau |= T(lbasis);
      }//END compute lattice basis
      
      if(uses_homog) Mtau = Mtau.minor(~scalar2set(0),All);
      
      //dbgtrace << "Appending result " << endl;
      
      //Compute kernel
      Matrix<Rational> Ktau = null_space(Mtau);
      //dbgtrace << "Kernel: " << Ktau << endl;
      //Compute the conversion to the total weight space:
      Set<int> remaining_classes = sequence(0,subdivision_classes.dim());
      Matrix<Rational> Ptau(Ktau.rows(), subdivision_classes.dim());
      for(int sig = 0; sig < sig_neighbours.dim(); sig++) {
	Ptau.col(class_index[sig_neighbours[sig]]) = Ktau.col(sig);
	remaining_classes -= (class_index[sig_neighbours[sig]]);
      }//END copy results
            
      //The remaining classes can have arbitrary weights here, so we append the appropriate unit vectors
      for(Entire<Set<int> >::iterator rc = entire(remaining_classes); !rc.at_end(); rc++) {
	Ptau /= unit_vector<Rational>(subdivision_classes.dim(),*rc);
      }//END add unit vectors for remaining classes
      
      //dbgtrace << "Transform: " << Ptau << endl;
      
      //Compute equations, attach to total system and reduce
      system_matrix /= null_space(Ptau);
      system_matrix = system_matrix.minor(basis_rows(system_matrix),All);
      
	
    }//END iterate codimension one faces
    
    //dbgtrace << "Transforming..." << endl;
    
    //To compute the final subdivision weight space, we add the equation that all zero classes
    // must have weight 0
    for(Entire<Set<int> >::iterator zc = entire(zero_classes); !zc.at_end(); zc++) {
      system_matrix /= unit_vector<Rational>(system_matrix.cols(), *zc);
    }//END add zero class equations
    
    //Compute the solution
    Matrix<Rational> subdiv_space = null_space(system_matrix);
    
    //Transform to weight space of actual complex
    Matrix<Rational> result(subdiv_space.rows(), maximal_cones.rows());
    for(int mc = 0; mc < maximal_cones.rows(); mc++) {
      result.col(mc) = subdiv_space.col(class_index[mc]);
    }//END transform result
    
    return result;
    
  }//END weight_properties
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
//   //Documentation see perl wrapper
//   Matrix<Integer> cycle_weight_space(perl::Object C) {
//     
//     //First we compute the weight system
//     Matrix<Integer> weight_system = cycle_weight_system(C);
//     //The weight solutions are the kernel...
//     Matrix<Rational> ker(null_space(weight_system));
//     //... but we have to make them a lattice basis
//     return latticeBasisFromRays(ker, Matrix<Rational>(0, ker.cols()), false);
//     /*
//     
//     IncidenceMatrix<> maximal_cones = C.give("MAXIMAL_CONES");
//     
//     //First we compute the irreducibility system
//     Map<int,int> index_translator;
//     Matrix<Integer> irred_system = cycle_irreducibility_system(C,index_translator);
//     
//     //Then compute the solution system
//     Matrix<Integer> ker = null_space(irred_system);
//     
//     Matrix<Rational> solutions(ker.rows(), maximal_cones.rows());
//     //Go through all kernel solutions
//     for(int k = 0; k < ker.rows(); k++) {
//       //Translate into weight solution
//       for(int m = 0; m< maximal_cones.rows(); m++) {
// 	solutions(k,m) = ker(k, index_translator[m]);
//       }
//     }
//     
//     //Now compute a lattice basis
//     return latticeBasisFromRays(solutions, Matrix<Rational>(0,solutions.cols()), false);*/
//     
//   }
  
  ///////////////////////////////////////////////////////////////////////////////////////
 /* 
  //Documentation see perl wrapper
  perl::Object cycle_weight_cone(perl::Object C) {
    
    //First we compute the weight system
    Matrix<Rational> wsys = C.give("WEIGHT_SYSTEM");
    
    //Create inequalities of the positive orthant
    Matrix<Rational> positive = unit_matrix<Rational>(wsys.cols());
    
    //Now create the cone
    perl::Object c("polytope::Cone");
      c.take("INEQUALITIES") << positive;
      c.take("EQUATIONS") << wsys;
      
    return c;
    
  }*/
  
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
    
//     Map<int,int> dummy;
    Matrix<Integer> F = C.give("WEIGHT_SPACE");
    return rank(F) == 1;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------

    Function4perl(&is_irreducible,"is_irreducible(WeightedComplex)");
//     Function4perl(&cycle_weight_space,"cycle_weight_space(WeightedComplex)");
//     Function4perl(&cycle_weight_system,"cycle_weight_system(WeightedComplex)");
//     Function4perl(&cycle_weight_cone,"cycle_weight_cone(WeightedComplex)");
    Function4perl(&cycle_weight_space,"cycle_weight_space(WeightedComplex)");

}}