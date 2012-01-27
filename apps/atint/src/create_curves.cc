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
 * This file contains functions to create tropical curves
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object rational_curve_embedding(const Matrix<Rational> &delta, const perl::Object &type) {
    //Extract values
    IncidenceMatrix<> nodes_by_sets = type.give("NODES_BY_SETS");
    IncidenceMatrix<> nodes_by_leaves = type.give("NODES_BY_LEAVES");
    IncidenceMatrix<> sets = type.give("SETS");
    Vector<int> coeffs = type.give("COEFFS");
    int n_leaves = type.give("N_LEAVES");	
    if(n_leaves != delta.rows()) {
      throw std::runtime_error("Cannot create curve embedding. Degree does not match number of leaves");
    }
    
    int ambient_dim = delta.cols();
    
    
    //Result variables
    //The first k rows contain the nodes in order of their appearance in
    //NODES_BY_*
    Matrix<Rational> rays(nodes_by_leaves.rows(),ambient_dim+1); 
    Vector<Set<int> > cones;
    
    //This vector tells us whether the position of a vertex has already been computed
    Array<bool> computed(nodes_by_leaves.rows());
    
    //We arbitrarily place the first node at (0,...,0)
    rays.row(0) = unit_vector<Rational>(rays.cols(),0);
    computed[0] = true;
    
    //Put all remaining nodes in a queue
    //Then iterate: When a node doesn't have a neighbour whose position has been computed yet,
    //move it to the back
    std::list<int> queue;
    for(int v = 1; v < nodes_by_sets.rows(); v++) {
      queue.push_back(v);
    }
    
    while(queue.size() > 0) {
      dbgtrace << "Computed: " << computed << endl;
      int nextv = queue.front();
	queue.pop_front();
      dbgtrace << "Looking for neighbours of " << nextv << endl;
      int neighbour = -1;
      for(int nb = 0; nb < computed.size(); nb++) {
	if(nb != nextv && computed[nb] && (nodes_by_sets.row(nb) * nodes_by_sets.row(nextv)).size() > 0) {
	    neighbour = nb; break;
	}
      }//END search for neighbours
      if(neighbour == -1) {
	queue.push_back(nextv); continue;
      }
      
      dbgtrace << "Neighbour is " << neighbour << endl;
      
      //Compute orientation of edge: Take any other edge at nextv. It the intersection with
      //the partition of the connecting edge is empty or the connecting set, then 
      //the edge points away from nextv, otherwise towards it
      int edge_index = *((nodes_by_sets.row(neighbour) * nodes_by_sets.row(nextv)).begin());
      Set<int> edge_set = sets.row(edge_index);
      Set<int> compare_set;
      if(nodes_by_leaves.row(nextv).size() > 0) {
	compare_set += *(nodes_by_leaves.row(nextv).begin());
      }
      else {
	int otherset = *((nodes_by_sets.row(nextv) - edge_index).begin());
	compare_set = sets.row(otherset);
      }
      int sign = +1;
      Set<int> inter = edge_set * compare_set;
      if(inter.size() == 0 || inter.size() == edge_set.size()) {
	sign = -1;
      }
      
      //Now compute the vertex
      Vector<Rational> sum_of_leaves(rays.cols()-1);
      for(Entire<Set<int> >::iterator d = entire(edge_set); !d.at_end(); d++) {
	sum_of_leaves += delta.row(*d-1);
      }
      sum_of_leaves = 0 | sum_of_leaves;
      rays.row(nextv) = rays.row(neighbour) + coeffs[edge_index] * sign * sum_of_leaves;
      computed[nextv] = true;
      
      //Create the cone
      Set<int> cone_set;
	cone_set += nextv; cone_set += neighbour;
      cones |= cone_set;
    }//END compute vertices
    
    //Finally we attach the leaves - but there may be doubles!
    Matrix<Rational> leaves = zero_vector<Rational>() | delta;
    Map<int, int> leaf_ray_index;
    for(int l = 0; l < leaves.rows(); l++) {
      int index = -1;
      for(int r = nodes_by_leaves.rows(); r < rays.rows(); r++) {
	if(rays.row(r) == leaves.row(l)) {
	  index = r; break;
	}
      }
      if(index == -1) {
	rays /= leaves.row(l);
	leaf_ray_index[l] = rays.rows()-1;
      }
      else {
	leaf_ray_index[l] = index;
      }
    }
    
    for(int v = 0; v < nodes_by_leaves.rows(); v++) {
      Set<int> leaves_here = nodes_by_leaves.row(v);
      if(leaves_here.size() > 0) {
	for(Entire<Set<int> >::iterator l = entire(leaves_here); !l.at_end(); l++) {
	    Set<int> cone_set;
	    cone_set += v;
	    cone_set += leaf_ray_index[*l-1];
	    cones |= cone_set;
	}
      }
    }
    
    Vector<Integer> weights = ones_vector<Integer>(cones.dim());
    
    dbgtrace << "Rays: " << rays << endl;
    dbgtrace << "Cones: "<< cones << endl;
    
    
    //Create result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << Matrix<Rational>(rays);
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("USES_HOMOGENEOUS_C") << true;
    
    return result;
    
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry / Rational curves"
		    "# This function creates an embedding of a rational tropical curve using"
		    "# a given abstract curve and degree"
		    "# @param Matrix<Rational> delta The degree of the curve in non-homogeneous "
		    "# coordinates. The number of rows"
		    "# should correspond to the number of leaves of type and the number of columns"
		    "# is the dimension of the space in which the curve should be realized"
		    "# @param RationalCurve type An abstract rational curve"
		    "# @return WeightedComplex The corresponding embedded complex, weighted"
		    "# with 1. The position of the curve is determined by the first node, "
		    "# which is always placed at the origin",
		    &rational_curve_embedding, "rational_curve_embedding(Matrix<Rational>, RationalCurve)");
  
}}