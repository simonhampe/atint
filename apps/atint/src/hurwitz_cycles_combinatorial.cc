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
 * Compute Hurwitz cycles combinatorially
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/PowerSet.h"
#include "polymake/permutations.h"
#include "polymake/atint/morphism_special.h"
#include "polymake/atint/specialvarieties.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/cdd_helper_functions.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a RationalCurve and a list of node indices. Then inserts additional leaves (starting from N_LEAVES+1) at these nodes and returns the resulting RationalCurve object
   @param perl::Object curve A RationalCurve object
   @param Vector<int> nodes A list of node indices of the curve
   */
  perl::Object insert_leaves(perl::Object curve, Vector<int> nodes) {
    //Extract values
    int max_leaf = curve.give("N_LEAVES");
    Vector<Set<int> > sets = curve.give("SETS");
    Vector<Rational> coeffs = curve.give("COEFFS");
    IncidenceMatrix<> nodes_by_sets = curve.give("NODES_BY_SETS");
    Vector<Set<int> > nodes_by_leaves = curve.give("NODES_BY_LEAVES");
    
    for(int n_el = 0; n_el < nodes.dim(); n_el++) {
      int n = nodes[n_el];
      max_leaf++;
      //Easy case: There is already a leaf at this vertex
      if(nodes_by_leaves[n].size() > 0) {
	int ref_leaf = *( nodes_by_leaves[n].begin());
	for(int s = 0; s < sets.dim(); s++) {
	    if(sets[s].contains(ref_leaf)) {
	      sets[s] += max_leaf;
	    }
	}
      }//END if there is already a leaf
      else {
	//Normalize the sets at the node to point away from it. Order them in an arbitrary way
	//Intersect the first two, I_1 and I_2: I_1 points away from the node, if and only if
	//I_1 cap I_2 = empty or I_1. The subsequent sets I_k only have to fulfill 
	//I_k cap I_1 = empty
	Vector<int> adjacent_sets(nodes_by_sets.row(n));
	Vector<Set<int> > normalized_sets;
	Set<int> first_inter = sets[adjacent_sets[0]] * sets[adjacent_sets[1]];
	normalized_sets |= 
	    (first_inter.size() == 0 || first_inter.size() == sets[adjacent_sets[0]].size()) ? sets[adjacent_sets[0]] : (sequence(1, max_leaf-1) - sets[adjacent_sets[0]]);
	for(int as = 1; as < adjacent_sets.dim(); as++) {
	    Set<int> subseq_inter = normalized_sets[0] * sets[adjacent_sets[as]];
	    normalized_sets |=
	    (subseq_inter.size() == 0)? sets[adjacent_sets[as]] : (sequence(1,max_leaf-1) - sets[adjacent_sets[as]]);
	}
	//Now for each set count, how many of the adjacent sets intersect it nontrivially
	//It points away from the node, if and only if this number is 1 (otherwise its all)
	for(int s = 0; s < sets.dim(); s++) {
	    int inter_count = 0;
	    for(int ns = 0; ns < normalized_sets.dim(); ns++) {
	      if( (sets[s] * normalized_sets[ns]).size() > 0) inter_count++;
	      if(inter_count > 1) break;
	    }
	    if(inter_count > 1) {
	     sets[s] += max_leaf; 
	    }	    
	}
      }//END if there are only bounded edges
      
      nodes_by_leaves[n] += max_leaf;
    }//END iterate nodes
    
    //Create result
    perl::Object result("RationalCurve");
      result.take("SETS") << sets;
      result.take("COEFFS") << coeffs;
      result.take("N_LEAVES") << max_leaf;
    return result;
    
  }//END insert_leaves
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a RationalCurve object and returns the rays v_I corresponding to its bounded edges as a matrix of row vectors in matroid coordinates.
   @param perl::Object curve A RationalCurve object
   @return Matrix<Rational> The rays corresponding to the bounded edges
   */
  Matrix<Rational> edge_rays(perl::Object curve) {
    IncidenceMatrix<> sets = curve.give("SETS");
    int n = curve.give("N_LEAVES");
    Matrix<Rational> result(0, n*(n-1)/2);
    for(int s = 0; s < sets.rows(); s++) {
      perl::Object rcurve("RationalCurve");
	rcurve.take("SETS") << sets.minor(scalar2set(s),All);
	rcurve.take("N_LEAVES") << n;
	rcurve.take("COEFFS") << ones_vector<Rational>(1);
      Vector<Rational> rray = rcurve.CallPolymakeMethod("matroid_vector");
      result /= rray;
    }
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes all ordered lists of k elements in {0,...,n-1}
   @param int n The size of the set {0,...,n-1}
   @param int k The size of the ordered subsets
   @return Matrix<int> All ordered k-sets as a list of row vectors, i.e. a matrix of dimension ( (n choose k) * k!) x k
   */
  Matrix<int> ordered_k_choices(int n, int k) {
    
    Matrix<int> result(0,k);
    
    //Compute all k-choices of the set {0,..,n-1}
    Array<Set<int> > kchoices = all_subsets_of_k(sequence(0,n),k);
    
    //Compute all permutations on a k-set
    AllPermutations<> kperm = all_permutations(k);
    
    for(int ch = 0; ch < kchoices.size(); ch++) {
      Vector<int> kvec(kchoices[ch]);
      for(Entire<AllPermutations<> >::iterator p = entire(kperm);!p.at_end(); p++) {
	result /= (permuted(kvec,*p));
      }
    }
    
    return result;
  }//END ordered_k_choices
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_subdivision(int k, Vector<int> degree, Vector<Rational> points = Vector<Rational>()) {
    
    //Make points default points if not given
    if(points.dim() < degree.dim() - 3- k) {
      points = points | zero_vector<Rational>( degree.dim()-3-k-points.dim());
    }
    
    int n = degree.dim();
    int big_moduli_dim = 2*n - k - 2;
    
    //Compute M_0,n and extract cones
    perl::Object m0n = CallPolymakeFunction("tropical_m0n",n);
    Matrix<Rational> mn_rays = m0n.give("RAYS");
    IncidenceMatrix<> mn_cones = m0n.give("MAXIMAL_CONES");
    
    //Create evaluation maps
    Matrix<Rational> ev_maps(0, big_moduli_dim);
    Matrix<Rational> rat_degree(n,0);
       rat_degree |= degree;
    for(int i = n+2; i <= 2*n-2-k; i++) {
      perl::Object evi = evaluation_map(n-2-k, 1, rat_degree, i-n-1);
      Matrix<Rational> evimatrix = evi.give("MATRIX");
      ev_maps /= (evimatrix.row(0).slice(sequence(0, evimatrix.cols()-1)));
    }
    
    //Will contain the rays/cones of the subdivided M_0,n in homog. coordinates
    Matrix<Rational> subdiv_rays(0,mn_rays.cols()+1);
    Vector<Set<int> > subdiv_cones;
    
    //Iterate all cones of M_0,n and compute their refinement as well as the 
    //Hurwitz cycle cones within
    for(int mc = 0; mc < mn_cones.rows(); mc++) {
	dbglog << "Refining cone " << mc << " of " << mn_cones.rows() << endl;
	//This will be a model of the subdivided cone of M_0,n
	Matrix<Rational> model_rays = unit_matrix<Rational>(degree.dim()-3);
	Vector<Set<int> > model_cones; model_cones |= sequence(0, degree.dim()-3);
	perl::Object model_complex("WeightedComplex");
	  model_complex.take("RAYS") << model_rays;
	  model_complex.take("MAXIMAL_CONES") << model_cones;
	model_complex = model_complex.CallPolymakeMethod("homogenize");  
	  
	//Extract the combinatorial type to compute evaluation maps
	perl::Object mc_type = CallPolymakeFunction("rational_curve_from_cone",m0n,degree.dim(),mc);
	  Vector<int> node_degrees = mc_type.give("NODE_DEGREES");
	  
	//Iterate over all possible ordered choices of (#vert-k)-subsets of nodes  
	Matrix<int> fix_node_sets = ordered_k_choices(node_degrees.dim(), node_degrees.dim()-k);
	dbgtrace << "Iterating all choices of fixed vertices " << endl;
	
	for(int nchoice = 0 ; nchoice < fix_node_sets.rows(); nchoice++) {
	  //Compute combinatorial type obtained by adding further ends to chosen nodes
	  perl::Object higher_type = insert_leaves(mc_type, fix_node_sets.row(nchoice));	  
	  //Convert evaluation maps to local basis of rays
	  Matrix<Rational> local_basis = edge_rays(higher_type);
	  Matrix<Rational> converted_maps = ev_maps * T(local_basis);
	  
	  //Now refine along each evaluation map
	  for(int evmap = 0; evmap < converted_maps.rows(); evmap++) {
	    //Compute half-space fan induced by the equation ev_i >= / = / <= p_i
	    //Then refine the model cone along this halfspace
	    perl::Object evi_halfspace = halfspace_complex(points[evmap],converted_maps.row(evmap));
	      evi_halfspace = evi_halfspace.CallPolymakeMethod("homogenize");
	    model_complex = refinement(model_complex, evi_halfspace, false,false,false,true,false).complex;
	  }//END iterate evaluation maps
	  
	}//END iterate choices of fixed nodes
	
	dbgtrace << "Re-converting refined cone " << endl;
	
	//Finally convert the model cone back to M_0,n-coordinates
	Matrix<Rational> model_subdiv_rays = model_complex.give("RAYS");
	IncidenceMatrix<> model_subdiv_cones = model_complex.give("MAXIMAL_CONES");
	//First convert the rays back
	Matrix<Rational> model_conv_rays = 
	  model_subdiv_rays.col(0) | (model_subdiv_rays.minor(All,~scalar2set(0)) * edge_rays(mc_type));
	Vector<int> model_rays_perm  = insert_rays(subdiv_rays, model_conv_rays,false,true);
	Map<int,int> ray_index_map;
	for(int mrp = 0; mrp < model_rays_perm.dim(); mrp++) {
	  ray_index_map[mrp] = model_rays_perm[mrp];
	}
	//Then use the above index list to transform the cones
	for(int msc = 0; msc < model_subdiv_cones.rows(); msc++) {
	  subdiv_cones |= attach_operation(model_subdiv_cones.row(msc), pm::operations::associative_access<Map<int,int>,int>(&ray_index_map));
	}

    }//END iterate cones of M_0,n
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << subdiv_rays;
      result.take("MAXIMAL_CONES") << subdiv_cones;
      result.take("TROPICAL_WEIGHTS") << ones_vector<Integer>(subdiv_cones.dim());
      result.take("USES_HOMOGENEOUS_C") << true;
      
    return result;
  }//END hurwitz_subdivision
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&insert_leaves,"insert_leaves(RationalCurve, Set<Int>)");
  Function4perl(&edge_rays,"edge_rays(RationalCurve)");
  
  UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
		    "# This function computes a subdivision of M_0,n containing the Hurwitz cycle"
		    "# H_k(x), x = (x_1,...,x_n) as a subfan. If k = n-4, this subdivision is the unique"
		    "# coarsest subdivision fulfilling this property"
		    "# @param Int k The dimension of the Hurwitz cycle, i.e. the number of moving vertices"
		    "# @param Vector<Int> degree The degree x. Should add up to 0"
		    "# @param Vector<Rational> points Optional. Should have length n-3-k. Gives the images of "
		    "# the fixed vertices (besides 0). If not given all fixed vertices are mapped to 0"
		    "# and the function computes the subdivision of M_0,n containing the recession fan of H_k(x)"
		    "# @return WeightedComplex A subdivision of M_0,n, in homogeneous coordinates",
		    &hurwitz_subdivision,"hurwitz_subdivision($,Vector<Int>;Vector<Rational> = new Vector<Rational>())");
}}