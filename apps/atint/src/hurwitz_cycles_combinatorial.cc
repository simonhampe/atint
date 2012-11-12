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
#include "polymake/linalg.h"
#include "polymake/atint/morphism_special.h"
#include "polymake/atint/specialvarieties.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/cdd_helper_functions.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
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
  
  /**
    @brief Takes a matrix and computes the gcd of the absolute values of the maximal minors
   */
  Integer gcd_maxminor(Matrix<Rational> map) {
    Array<Set<int> > r_choices = all_subsets_of_k(sequence(0,map.cols()),map.rows());
    int g = 0;
    for(int r = 0; r < r_choices.size(); r++) {
      g = gcd(g,det(map.minor(All,r_choices[r])).to_int());
    }
    return Integer(abs(g));
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief This takes a list of cones (in terms of ray indices) and their weights and inserts another cone wit given weight. If the new cone already exists in the list, nothing happens
   */
  void insert_cone(Vector<Set<int> > &cones, Vector<Integer> &weights, Set<int> ncone, Integer nweight) {
    int ncone_index = -1;
    for(int c = 0; c < cones.dim(); c++) {
      Set<int> inter = ncone * cones[c];
      if( inter.size() == ncone.size() && inter.size() == cones[c].size()) {
	ncone_index = c; break;
      }
    }
    if(ncone_index == -1) {
      cones |= ncone;
      weights |= nweight;
    }
    else {
      weights[ncone_index] += nweight;
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Checks whether the cone of the positive orthant intersects all of the hyperplanes defined by ev_maps.row(i) = p_i in its interior (i.e. is not contained in one of the half spaces >= pi / <= pi), where p_i = 1,2,... This can be used to check whether a choice of fixed vertices is valid (because the corresponding refinement has to have at least one interior cone for a generic choice of points p_i). In fact, this is equivalent to saying that each row of the matrix has at least one positive entry and that the hyperplanes intersect the orthant in a cone, s.t. for each entry _j there is at least one ray or vertex v s.t. v_j > 0.
   */
  bool is_valid_choice(Matrix<Rational> ev_maps) {
    //Check if all maps have a positive entry.
    for(int r = 0; r < ev_maps.rows(); r++) {
      bool found_positive = false;
      for(int c = 0; c < ev_maps.cols(); c++) {
	if(ev_maps(r,c) > 0) {
	    found_positive = true; break;
	}
      }
      if(!found_positive) return false;
    }
    
    //Now compute the intersection of all hyperplanes with the positive orthant
    Matrix<Rational> ineq = unit_matrix<Rational>(ev_maps.cols());
      ineq = zero_vector<Rational>() | ineq;
    Matrix<Rational> eq = ev_maps;
      eq = (- Vector<Rational>(sequence(1,ev_maps.rows()))) | eq;
    try {
      Matrix<Rational> rays = solver<Rational>().enumerate_vertices(ineq,eq, false,true).first;
      bool found_positive = false;
      for(int c = 1; c < rays.cols(); c++) {
	for(int r = 0; r < rays.rows(); r++) {
	    if(rays(r,c) > 0) {
	      found_positive = true; break;
	    }
	}
	if(!found_positive) return false;
      }
    }
    catch(...) {
      //This might happen, if the system has no solution
      return false;
    }
    return true;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes either a subdivision of M_0,n containing H_k(degree) or the cycle itself. The function 
   */
  perl::Object hurwitz_computation(int k, Vector<int> degree, Vector<Rational> points, bool compute_only_subdivision) {
    
    //Make points default points ( = 0) if not given and find out if the points are generic
    bool is_generic = true;
    if(points.dim() < degree.dim() - 3- k) {
      points = points | zero_vector<Rational>( degree.dim()-3-k-points.dim());
      is_generic = false;
    }
    if( (Set<Rational>(points)).size() < points.dim()) is_generic = false;
    
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
    Vector<Integer> subdiv_weights;
    
    //Iterate all cones of M_0,n and compute their refinement as well as the 
    //Hurwitz cycle cones within
    for(int mc = 0; mc < mn_cones.rows(); mc++) {
	pm::cout << "Refining cone " << mc << " of " << mn_cones.rows() << endl;
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
	//dbgtrace << "Iterating all choices of fixed vertices " << endl;
	
	//We save the evaluation matrices for use in the computation
	//of tropical weights in the Hurwitz cycle
	Vector<Matrix<Rational> > evmap_list;
	
	for(int nchoice = 0 ; nchoice < fix_node_sets.rows(); nchoice++) {
	  //Compute combinatorial type obtained by adding further ends to chosen nodes
	  perl::Object higher_type = insert_leaves(mc_type, fix_node_sets.row(nchoice));
	  std::string s = higher_type.CallPolymakeMethod("to_string");
	  //dbgtrace << "Intersecting with type " << s << endl;
	  //Convert evaluation maps to local basis of rays
	  Matrix<Rational> local_basis = edge_rays(higher_type);
	  Matrix<Rational> converted_maps = ev_maps * T(local_basis);
	  
	  //Check if this choice of fixed vertices induces some valid Hurwitz type
	  //(in the generic case)
	  if(!is_valid_choice(converted_maps)) continue;
	  
	  evmap_list |= converted_maps;
	  
	  //Now refine along each evaluation map
	  for(int evmap = 0; evmap < converted_maps.rows(); evmap++) {
	    //Compute half-space fan induced by the equation ev_i >= / = / <= p_i
	    //Then refine the model cone along this halfspace
	    perl::Object evi_halfspace = halfspace_complex(points[evmap],converted_maps.row(evmap));
	      evi_halfspace = evi_halfspace.CallPolymakeMethod("homogenize");
	    model_complex = refinement(model_complex, evi_halfspace, false,false,false,true,false).complex;
	  }//END iterate evaluation maps
	  
	}//END iterate choices of fixed nodes
	
	//If we are actually computing the Hurwitz cycle, we take the evalutation map
	//corresponding to each choice of vertex placement. We then compute the k-cones
	//in the subdivision on which this map is (p_1,..,p_r) and add the gcd of the maximal minors of the
	//matrix as tropical weigh to these cones
	IncidenceMatrix<> model_hurwitz_cones;
	Vector<Integer> model_hurwitz_weights;
	Matrix<Rational> model_hurwitz_rays;
	if(!compute_only_subdivision) {
	  perl::Object skeleton = CallPolymakeFunction("skeleton_complex",model_complex,k,false);
	  //We go through all dimension - k cones of the subdivision
	  Matrix<Rational> k_rays = skeleton.give("RAYS");
	  Vector<Set<int > > k_cones = skeleton.give("MAXIMAL_CONES");
	    model_hurwitz_weights = zero_vector<Integer>(k_cones.dim());
	    //dbgtrace << "Rays: " << k_rays << endl;
	    //dbgtrace << "Cones: " << k_cones << endl;
	  Set<int> non_zero_cones;
	  for(int m = 0; m < evmap_list.dim(); m++) {
	    //Compute gcd of max-minors
	    //dbgtrace << "Matrix is " << evmap_list[m] << endl;
	    Integer g = gcd_maxminor(evmap_list[m]);
	    for(int c = 0; c < k_cones.dim(); c++) {
	      //Check if ev maps to the p_i on this cone (rays have to be mapped to 0!)
	      Set<int> cset = k_cones[c];
	      bool maps_to_pi = true;
	      for(Entire<Set<int> >::iterator cit = entire(cset);!cit.at_end(); cit++) {
		Vector<Rational> cmp_vector = 
		  ( k_rays(*cit,0) == 1? points : zero_vector<Rational>(points.dim())) ;
		if( evmap_list[m] * k_rays.row(*cit).slice(~scalar2set(0)) != cmp_vector ){
		    maps_to_pi = false;break;
		}
	      }
	      if(maps_to_pi) {
		model_hurwitz_weights[c] += g;
		if(g != 0) non_zero_cones += c;
		//dbgtrace << "Adding weight to cone " << c << endl;
	      }
	    }
	  }//END iterate evaluation maps
	  //Identify cones with non-zero weight and the rays they use
	  Set<int> used_rays = accumulate(k_cones.slice(non_zero_cones),operations::add());
	  model_hurwitz_rays = k_rays.minor(used_rays,All);
	  model_hurwitz_cones = IncidenceMatrix<>(k_cones).minor(non_zero_cones,used_rays);
	  model_hurwitz_weights = model_hurwitz_weights.slice(non_zero_cones);
	  
	  //dbgtrace << "Model weight: " << model_hurwitz_weights << endl;
	  
	}//END compute hurwitz weights
	
	//dbgtrace << "Re-converting refined cone " << endl;
	
	//Finally convert the model cone back to M_0,n-coordinates
	Matrix<Rational> model_subdiv_rays = model_complex.give("RAYS");
	IncidenceMatrix<> model_subdiv_cones = model_complex.give("MAXIMAL_CONES");
	Vector<Integer> model_subdiv_weights = ones_vector<Integer>(model_subdiv_cones.rows());
	//If we want to compute the Hurwitz cycle, we take the k-cones instead
	if(!compute_only_subdivision) {
	  model_subdiv_rays = model_hurwitz_rays;
	  model_subdiv_cones = model_hurwitz_cones;
	  model_subdiv_weights = model_hurwitz_weights;
	}
	//First convert the rays back
	//dbgtrace << "Final rays are: " << model_subdiv_rays << endl;
	Matrix<Rational> model_conv_rays = 
	  model_subdiv_rays.col(0) | (model_subdiv_rays.minor(All,~scalar2set(0)) * edge_rays(mc_type));
	  //dbgtrace << "Converted rays are: " << model_conv_rays << endl;
	Vector<int> model_rays_perm  = insert_rays(subdiv_rays, model_conv_rays,false,true);
	Map<int,int> ray_index_map;
	for(int mrp = 0; mrp < model_rays_perm.dim(); mrp++) {
	  ray_index_map[mrp] = model_rays_perm[mrp];
	}
	//Then use the above index list to transform the cones
	for(int msc = 0; msc < model_subdiv_cones.rows(); msc++) {
	  insert_cone(subdiv_cones, subdiv_weights, 
		      attach_operation(model_subdiv_cones.row(msc), pm::operations::associative_access<Map<int,int>,int>(&ray_index_map)), 
		      model_subdiv_weights[msc]);
	}

    }//END iterate cones of M_0,n
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << subdiv_rays;
      result.take("MAXIMAL_CONES") << subdiv_cones;
      result.take("TROPICAL_WEIGHTS") << subdiv_weights;
      result.take("USES_HOMOGENEOUS_C") << true;
      
    return result;
  }//END hurwitz_computation
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_subdivision(int k, Vector<int> degree, Vector<Rational> points = Vector<Rational>()) {
    return hurwitz_computation(k,degree,points,true);
  }
  
  //Documentation see perl wrapper
  perl::Object hurwitz_cycle(int k,  Vector<int> degree, Vector<Rational> points = Vector<Rational>()) {
    return hurwitz_computation(k,degree,points,false);
  }
  
  
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
  
  UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
		    "# This function computes the Hurwitz cycle H_k(x), x = (x_1,...,x_n)"
		    "# @param Int k The dimension of the Hurwitz cycle, i.e. the number of moving vertices"
		    "# @param Vector<Int> degree The degree x. Should add up to 0"
		    "# @param Vector<Rational> points Optional. Should have length n-3-k. Gives the images of "
		    "# the fixed vertices (besides 0). If not given all fixed vertices are mapped to 0"
		    "# and the function computes the recession fan of H_k(x)"
		    "# @return WeightedComplex H_k(x), in homogeneous coordinates",
		    &hurwitz_cycle, "hurwitz_cycle($,Vector<Int>;Vector<Rational> = new Vector<Rational>())");
		    
}}