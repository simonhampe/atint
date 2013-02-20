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
 * This file contains functions to compute Hurwitz cycles in M_0,n
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/Graph.h"
#include "polymake/atint/psi_classes.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/morphism_special.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
//   using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_pre_cycle(int k, Vector<int> degree, Vector<Rational> pullback_points = Vector<Rational>()) {
    
    solver<Rational> sv;
    
    //First, compute the psi-class product
    int n = degree.dim();
    
    if(pullback_points.dim() < n-3-k) {
      pullback_points |= zero_vector<Rational>(n-3-k - pullback_points.dim());
    }
    
    int big_moduli_dim = 2*n - k - 2;
    Vector<int> exponents = zero_vector<int>(n) | ones_vector<int>(n-2-k);
    //dbgtrace << "Computing psi product in M_N, N = " << big_moduli_dim << " with exponents " << exponents << endl;
    perl::Object P = psi_product(big_moduli_dim,exponents);
    
    if(n == 4) return P;
    
    //Compute evalutation maps and pullbacks
    Matrix<Rational> ev_maps(0, big_moduli_dim);
    Vector<perl::Object> pb_functions;
    Matrix<Rational> rat_degree(degree.dim(),0);
      rat_degree |= degree;
    for(int i = n+2; i <= 2*n-2-k; i++) {
      //dbgtrace << "Computing evaluation map pull back for i = " << i-n-1 << endl;
      perl::Object evi = evaluation_map(n-2-k, 1, rat_degree, i-n-1);
      Matrix<Rational> evi_matrix = evi.give("MATRIX");
      
      //We save the linear map to R defined by ev_i in a separate matrix for sorting out cones
      ev_maps /= (evi_matrix.row(0).slice(sequence(0, evi_matrix.cols()-1)));
      
      //Pulling back p_i = max(x,p_i) * R means we take the vector representing the morphism and 
      //attach a row below that has p_i at the end
      evi_matrix /= zero_vector<Rational>(evi_matrix.cols());
      evi_matrix(1, evi_matrix.cols()-1) = pullback_points[i-n-2];
      //Since we restrict ourselves to M_0,N x {0}, we actually ignore the last coefficient
      //of ev_i and replace it by the constant coefficient 0 (for the min-max-function)
      evi_matrix(0, evi_matrix.cols()-1) = 0;
      //dbgtrace << "Pullback evaluation matrix is " << evi_matrix << endl;
      perl::Object pb("MinMaxFunction");
	pb.take("FUNCTION_MATRIX") << Matrix<Rational>(evi_matrix);
	pb.take("USES_MIN") << false;

      pb_functions |= pb;
    }//END compute pullback functions
    
    //dbgtrace << "Throwing out cones" << endl;
    //Before we actually compute the divisor of the pullback of points along evaluation maps, we sort out 
    //those cones, that do not even map to these points
    Matrix<Rational> psi_rays = P.give("RAYS");
    IncidenceMatrix<> psi_cones = P.give("MAXIMAL_CONES");
    Vector<Integer> psi_weights = P.give("TROPICAL_WEIGHTS");
    
    ev_maps = T(ev_maps);
    Matrix<Rational> inequalities = unit_matrix<Rational>(big_moduli_dim-3 - (n-2-k));
	inequalities = zero_vector<Rational>() | inequalities;

    //dbgtrace << "Evaluation maps: " << ev_maps << "\npos. orthant: " << inequalities << endl;
	
    Set<int> kept_cones = sequence(0, psi_cones.rows());
    //FIXME: This doesn't seem to work for non-generic points?? (or k = 0?)
//     if(k <= 1) {
//       kept_cones = Set<int>();
//       for(int c = 0; c < psi_cones.rows(); c++) {
// 	dbgtrace << "Computing preimage of points in cone " << c << " of " << psi_cones.rows() << endl;
// 	//Compute the intersection of the preimages of all points p_i (as a cone in the current cone)
// 	//If it is empty, throw out this cone
// 	dbgtrace << "Rays: " << psi_rays.minor(psi_cones.row(c),All) << endl;
// 	Matrix<Rational> equalities = T( psi_rays.minor(psi_cones.row(c),All) * ev_maps);
// 	equalities = (-Vector<Rational>(pullback_points)) | equalities;
// 	//We intersect with the positive orthant, since we only allow positive scalar combinations of rays
// 	dbgtrace << "Eq: " << equalities << endl;
// 
// 	try {
// 	  Matrix<Rational> sol_rays = (sv.enumerate_vertices(inequalities, equalities, true,true)).first;
// 	  if(sol_rays.rows() > 0) {
// 	    dbgtrace << "Solution: " << sol_rays << endl;
// 	    kept_cones += c;
// 	  }
// 	}
// 	catch(...){ //Catch any "infeasible-system"-exception
// 	  //Go to the next cone
// 	  continue;
// 	}
// 	
// 	
//       }//END sort out cones
//     }
    //dbgtrace << "Cones remaining: " << kept_cones.size() << " out of " << psi_cones.rows() << endl;
    
    Set<int> used_rays = accumulate(rows(psi_cones.minor(kept_cones,All)), operations::add());
    P = perl::Object("WeightedComplex");
      P.take("RAYS") << (psi_rays.minor(used_rays,All));
      P.take("MAXIMAL_CONES") << (psi_cones.minor(kept_cones,used_rays));
      P.take("TROPICAL_WEIGHTS") << (psi_weights.slice(kept_cones));
      P.take("IS_UNIMODULAR") << true;
   
    //Now compute the divisor
    for(int i = 0; i < pb_functions.dim(); i++) {
      //dbgtrace << "Computing divisor of function " << (i+1) << " of " << pb_functions.dim() << endl;
      P = divisor_minmax(P, pb_functions[i]);
    }//END compute divisors
    
    return P;
    
  }//END function hurwitz_pre_cycle
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_push_cycle(int k, Vector<int> degree, Vector<Rational> points) {
    int n = degree.dim();
    int big_n = 2*n - k - 2;
    
    //Check points
    if(points.dim() < n-k-3) {
      throw std::runtime_error("Cannot compute Hurwitz cycle. Not enough points given.");
    }
    for(int i = 0; i < points.dim(); i++) {
      for(int j = i+1; j < points.dim(); j++) {
	if(points[i] == points[j] || points[i] == 0) {
	  throw std::runtime_error("Cannot compute Hurwitz cycle. Need distinct points != 0 to allow easy push-forward. Try different points or hurwitz_pre_cycle(...)");
	}
      }
    }
    
    //First we compute the pre-cycle    
    perl::Object precycle =  hurwitz_pre_cycle(1, degree, points);
    Matrix<Rational> rays = precycle.give("RAYS");
    IncidenceMatrix<> cones = precycle.give("MAXIMAL_CONES");
    Vector<Integer> weights = precycle.give("TROPICAL_WEIGHTS");
    bool uses_homog = precycle.give("USES_HOMOGENEOUS_C");
    
    //Then apply the forgetful map to the rays
    perl::Object ffmap = CallPolymakeFunction("forgetful_map",big_n, sequence(degree.dim()+1, n-k-2));
    Matrix<Rational> ffmatrix = ffmap.give("MATRIX");
    //If we use homog. coordiantes, keep the homogenizing coord.
    Vector<Rational> homog_coord = rays.col(0);
    if(uses_homog) rays = rays.minor(All,~scalar2set(0));
    
    rays = T(ffmatrix * T(rays));
    if(uses_homog) rays = homog_coord | rays;
    
    //Check for double rays
    Map<int,int> ray_map;
    Set<int> used_rays;
    for(int r = 0; r < rays.rows(); r++) {
      ray_map[r] = r;
      for(int s = 0; s < r; s++) {
	if(rays.row(r) == rays.row(s)) {
	    ray_map[r] = s; break;
	}
      }
      if(ray_map[r] == r) used_rays += r;
    }//END check doubles
    
    //Translate cones to new ray indices
    Vector<Set<int> > translated_cones;
    for(int oc = 0; oc < cones.rows(); oc++) {
      Set<int> tcone = 
	attach_operation(cones.row(oc), pm::operations::associative_access<Map<int,int>, int>(&ray_map));
      translated_cones |= tcone;
    }
    
    //Now we check if two cones intersect transversally and replace them accordingly
    solver<Rational> sv;
    bool added_something = true;
    while(added_something) {
      added_something = false;
      //Go through all cone pairs of cones that don't intersect in a vertex already
      for(int first_cone = 0; first_cone < translated_cones.dim() && !added_something; first_cone++) {
	Matrix<Rational> fconerays = rays.minor(translated_cones[first_cone],All);
	for(int sec_cone = first_cone+1; sec_cone < translated_cones.dim(); sec_cone++) {
	  if( (translated_cones[first_cone] * translated_cones[sec_cone]).size() == 0) {
	    //Compute H-reps
	    std::pair<Matrix<Rational>, Matrix<Rational> > frep = 
	      sv.enumerate_facets( zero_vector<Rational>() | fconerays, Matrix<Rational>(0, rays.cols()+1), true,false);
	    std::pair<Matrix<Rational>, Matrix<Rational> > srep =
	      sv.enumerate_facets( zero_vector<Rational>() | rays.minor(translated_cones[sec_cone],All),
				   Matrix<Rational>(0,rays.cols()+1), true,false);
	    //Compute intersection
	    Matrix<Rational> vertex_matrix = sv.enumerate_vertices( frep.first / srep.first, 
							     frep.second / srep.second,
							    true,true).first.minor(All,~scalar2set(0));
							    
	    if(vertex_matrix.rows() == 0 || vertex_matrix(0,0) == 0) continue;
	    
	    //Normalize
	    Vector<Rational> vertex = vertex_matrix.row(0) / vertex_matrix(0,0);
	   
	    rays = rays / vertex;
	    int vertex_index = rays.rows() -1;
	    used_rays += vertex_index;
	    //Replace cones
	    added_something = true;
	    
	    Integer fweight = weights[first_cone];
	    Integer sweight = weights[sec_cone];
	    Vector<int> fset(translated_cones[first_cone]);
	    Vector<int> sset(translated_cones[sec_cone]);
	    //Remove old cones
	    Set<int> rem_cones; rem_cones += first_cone; rem_cones += sec_cone;
	    translated_cones = translated_cones.slice(~rem_cones);
	    weights = weights.slice(~rem_cones);
	    //Add new ones
	    Set<int> a1,a2,b1,b2;
	    a1 += fset[0]; a1 += vertex_index;
	    a2 += fset[1]; a2 += vertex_index;
	    b1 += sset[0]; b1 += vertex_index;
	    b2 += sset[1]; b2 += vertex_index;
	    translated_cones |= a1; translated_cones |= a2;
	    translated_cones |= b1; translated_cones |= b2;
	    weights |= fweight; weights |= fweight;
	    weights |= sweight; weights |= sweight;
	    
	    break;
// 	    }
// 	    catch(...) {//Catch any "infeasible system" exceptions
// 	      //Do nothing
// 	    }
	    
	  }//END if intersection is empty
	}//END iterate second cone
      }//END iterate first cone
    }//END check transversal intersections
    

    IncidenceMatrix<> cone_matrix(translated_cones);
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays.minor(used_rays,All);
      result.take("MAXIMAL_CONES") << cone_matrix.minor(All, used_rays);
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
      
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_curve(Vector<int> degree) {
    //First we compute the pre-cycle
    perl::Object precycle =  hurwitz_pre_cycle(1, degree, Vector<Rational>(sequence(1,degree.dim()-4)));
    
    //Then we extract values
    Matrix<Rational> rays = precycle.give("RAYS");
    IncidenceMatrix<> cones = precycle.give("MAXIMAL_CONES");
    Vector<Integer> weights = precycle.give("TROPICAL_WEIGHTS");
    bool uses_homog = precycle.give("USES_HOMOGENEOUS_C");
    Set<int> directional = precycle.give("DIRECTIONAL_RAYS");
    
    //Result variables
    Matrix<Rational> curve_rays(0,rays.cols());
    Vector<Set<int> > curve_cones;
    Vector<Integer> curve_weights;
    
    //First of all, we remove the homogenizing coordinate
    if(uses_homog) {
//       int homog_index = -1;n
//       for(int r = 0; r < rays.rows(); r++) {
// 	if(rays(r,0) == 1) {
// 	  homog_index = r; break;
// 	}
//       }
//       rays = rays.minor(~scalar2set(homog_index),~scalar2set(0));
//       cones = cones.minor(All,~scalar2set(homog_index));
      rays = rays.minor(All,~scalar2set(0));
    }
    
    //Then we apply the forgetful map to the rays
    perl::Object ffmap = CallPolymakeFunction("forgetful_map",2*degree.dim() -3, sequence(degree.dim()+1, degree.dim()-3));
    Matrix<Rational> ffmatrix = ffmap.give("MATRIX");
    rays = rays * T(ffmatrix);
    
    //Normalize rays
    for(int r = 0; r < rays.rows(); r++) {
      for(int c = 0; c < rays.cols(); c++) {
	if(rays(r,c) != 0) {
	    rays.row(r) /= abs(rays(r,c));
	    break;
	}
      }
    }
    
    //Iterate all cones
    for(int c = 0; c < cones.rows(); c++) {
      Set<int> cdir = cones.row(c) * directional;
      if(cdir.size() == 0) continue;
      //The cone now only consists of one ray
      int r_index = cdir.front();
      //Check if this ray already exists
      int n_index = -1;
      for(int nray = 0; nray < curve_rays.rows(); nray++) {
	if(curve_rays.row(nray) == rays.row(r_index)) {
	  n_index = nray; break;
	}
      }
      
      //If it doesn't exist, add it
      if(n_index == -1) {
	curve_rays /= rays.row(r_index);
	Set<int> single_set;
	single_set += (curve_rays.rows()-1);
	curve_cones |= single_set;
	curve_weights |= weights[c];
      }
      //Otherwise add its weight to the appropriate cone 
      //(whose index is now equal to the ray index)
      else {
	curve_weights[n_index] += weights[c];
      }
      
    }//END iterate cones
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << curve_rays;
      result.take("MAXIMAL_CONES") << curve_cones;
      result.take("TROPICAL_WEIGHTS") << curve_weights;
      
    return result;
    
  }//END function hurwitz_curve
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  Integer hurwitz_degree(Vector<int> degree) {
    //First we compute the pre-cycle
    perl::Object precycle = hurwitz_pre_cycle(0, degree);
    
    Vector<Integer> weights = precycle.give("TROPICAL_WEIGHTS");
    
    return accumulate(weights, operations::add()); 	
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object hurwitz_graph(perl::Object cycle) {
    //Homogenize
    cycle = cycle.CallPolymakeMethod("homogenize");
    //Extract values
    Matrix<Rational> rays = cycle.give("RAYS");
      rays = rays.minor(All,~scalar2set(0));
    IncidenceMatrix<> cones = cycle.give("MAXIMAL_CONES");
    IncidenceMatrix<> rays_in_cones = T(cones);
    Set<int> vertices = cycle.give("VERTICES");
    Set<int> dirrays = cycle.give("DIRECTIONAL_RAYS");
    
    Vector<std::string> labels;
    
    //Compute how many vertices the graph has: #vertices + sum of directional rays at each vertex
    //Also create vertex labels
    int nodes = vertices.size();
    Vector<Set<int> > dir_at_vertices;
    for(Entire<Set<int> >::iterator v = entire(vertices); !v.at_end(); v++) {
      Set<int> rays_at_v = accumulate(rows(cones.minor(rays_in_cones.row(*v),All)),operations::add());
	rays_at_v *= dirrays;
      dir_at_vertices |= rays_at_v;
      nodes += rays_at_v.size();
      
      perl::Object vcurve = CallPolymakeFunction("rational_curve_from_moduli",rays.row(*v));
      std::string vstring = vcurve.CallPolymakeMethod("to_string");
      labels |= vstring;
    }
    
    //Graph object. The first #vertices nodes correspond to vertices, the rest are directional rays
    Graph<> G(nodes);
    
    
    //First create the bounded part
    IncidenceMatrix<> bounded_part = T(rays_in_cones.minor(vertices,All));
    for(int r = 0; r < bounded_part.rows(); r++) {
      if(bounded_part.row(r).size() > 1) {
	Vector<int> bp(bounded_part.row(r));
	G.edge( bp[0], bp[1] );
      }
    }//END create bounded part
    
    //Attach directional rays and create remaining labels
    int next_node_index = vertices.size();
    for(int v = 0; v < vertices.size(); v++) {
      Vector<int> raylist(dir_at_vertices[v]);
      for(int i = 0; i < raylist.dim(); i++) {
	G.edge(v,next_node_index);
	next_node_index++;
	
	perl::Object icurve = CallPolymakeFunction("rational_curve_from_moduli",rays.row(raylist[i]));
	std::string istring = icurve.CallPolymakeMethod("to_string");
	labels |= istring;
      }
    }//END attach directional rays
    
    perl::Object result("graph::Graph");
      result.take("N_NODES") << nodes;
      result.take("ADJACENCY") << G;
      result.take("NODE_LABELS") << labels;
      
    return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
//   UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
// 		    "# Computes the k-dimensional tropical Hurwitz cycle H_k(degree), except that it doesn't"
// 		    "# compute the push-forward to M_0,n"
// 		    "# @param int k The dimension of the Hurwitz cycle"
// 		    "# @param Vector<Int> degree The degree of the covering. The sum over all entries should "
// 		    "# be 0 and if n := degree.dim, then 0 <= k <= n-3"
// 		    "# @param Vector<Rational> pullback_points The points p_i that should be pulled back to "
// 		    "# determine the Hurwitz cycle (in addition to 0). Should have length n-3-k. If it is not given, all p_i"
// 		    "# are by default equal to 0 (same for missing points)"
// 		    "# @return perl::Object A WeightedComplex object representing the Hurwitz cycle H_k(degree) before push-forward",    
// 		    &hurwitz_pre_cycle, "hurwitz_pre_cycle($, Vector<Int>; Vector<Rational> = new Vector<Rational>())");
// 
// //   UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
// // 		    "# Computes the k-dimensional tropical Hurwitz cycle H_k(degree), including"
// // 		    "# push-forward to M_0,n"
// // 		    "# @param Int k The dimension of the Hurwitz cycle"
// // 		    "# @param Vector<Int> degree The degree of the covering. The sum over all entries should "
// // 		    "# be 0 and if n := degree.dim, then 0 <= k <= n-3"
// // 		    "# @param Vector<Rational> pullback_points The points p_i that should be pulled back to "
// // 		    "# determine the Hurwitz cycle. Should have length n-3-k and BE ALL DISTINCT AND NONZERO."
// // 		    "# If points are missing or some are equal, an error is thrown"
// // 		    "# @return perl::Object A WeightedComplex object representing the Hurwitz cycle H_k(degree)",    
// // 		    &hurwitz_cycle, "hurwitz_cycle($, Vector<Int>, Vector<Rational>)");
//   
//   UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
// 		    "# Computes the Hurwitz curve H_1(degree)"
// 		    "# @param Vector<int> degree The degree of the covering. The sum over all entries should "
// 		    "# be 0 and if n := degree.dim, then 0 <= 1 <= n-3"
// 		    "# @return perl::Object A WeightedComplex object representing the Hurwitz cycle H_1(degree). This will always be a fan cycle", 
// 		    &hurwitz_curve, "hurwitz_curve(Vector<Int>)");
/*  
  UserFunction4perl("# @category Tropical geometry / Hurwitz cycles"
		    "# Computes the Hurwitz degree H_0(degree)"
		    "# @param Vector<int> degree The degree of the covering. The sum over all entries should "
		    "# be 0 and if n := degree.dim, then 0 <= n-3"
		    "# @return Integer The Hurwitz degree H_0(degree)", 
		    &hurwitz_degree, "hurwitz_degree(Vector<Int>)");*/
  
  UserFunction4perl("# @category Hurwitz cycles"
		    "# Takes as input a hurwitz curve and computes the corresponding graph object"
		    "# Directional rays are modeled as terminal vertices. Each vertex (including directional rays)"
		    "# is labeled with its combinatorial type"
		    "# @param WeightedComplex cycle A Hurwitz curve object (or any one-dimensional cycle in M_0,n)"
		    "# @return graph::Graph",
		    hurwitz_graph,"hurwitz_graph(WeightedComplex)");
  
  Function4perl(&hurwitz_push_cycle, "hpc($, Vector<Int>, Vector<Rational>)");
  
}}