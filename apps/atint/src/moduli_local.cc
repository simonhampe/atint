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
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/moduli.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes all paired ordered Prüfer sequences of order n, i.e. all sequences of length 2n-4 where each element from n,...,2n-3 occurs twice and where removing every second occurence of an element yields an ascending sequence)
   @param int n Should be greater equal 3
   @return Vector<Vector<int> > A list of index sequences
   */
  Vector<Vector<int> > computePrueferSequences(int n) {
    if(n < 3) throw std::runtime_error("Cannot compute M_n cones for n < 3");
    
    //Compute number of sequences
    int noOfCones = count_mn_cones(n);
    
    Vector<Vector<int> > result;
    
    //The idea of the algorithm is the following: Indices describes the Prüfer sequence in the 
    //following way:
    //If the entry at position i is j, then the second occurence of i is placed at the (j+1)-st
    // *possible* (i.e. not yet filled) position. Note that the first occurence of each i is alredy
    // determined by the placement of all entries lower than i.
    Vector<int> indices = ones_vector<int>(n-2);
    for(int iteration = 0; iteration < noOfCones; iteration++) {
      Vector<int> baseSequence = zero_vector<int>(2*n -4);
      for(int i = 0; i < n-1; i++) {
	//Go through the zero entries of baseSequence. If it is the first or the indices[i]+1-th, 
	//insert an n+i
	int zero_count = -1;
	for(int entry = 0; entry < baseSequence.dim(); entry++) {
	  if(baseSequence[entry] == 0) {
	      zero_count++;
	      if(zero_count == 0) {
		baseSequence[entry] = n+i;
	      }
	      if(zero_count == indices[i]) {
		baseSequence[entry] = n+i;
		break;
	      }
	  }
	}
      }
      result |= baseSequence;
      //Increase the indices vector by "1"---------------------------------------------------    
      if(iteration < noOfCones-1) {
	int counterindex = n-3;
	while(indices[counterindex] == 2*(n-counterindex)-5) {
	  indices[counterindex] = 1;
	  counterindex--;
	}
	indices[counterindex]++;
      }
    }
    
    return result;
  } //END computePrueferSequences
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a ordered paired Prüfer sequence (as for example computed by computePrueferSequences and decodes it into a list of edge partitions of the combinatorial type
   @param Vector<int> pseq An ordered paired Prüfer sequence
   @return Vector<Set<int> > A list of partitions of the corresponding combinatorial type
   */
  Vector<Set<int> > decodePrueferSequence(const Vector<int> &pseq) {
    //Construct vertex set
    int n = pseq[0]; //The first element is always the number of leaves
    Set<int> V = sequence(0,2*n-2);
    Vector<Set<int> > adjacencies(n-2); //Which leaves lie "behind" which interior vertex?
    Vector<Set<int> > result;
    
    dbgtrace << "Connecting leaves" << endl;
    int firstindex = 0; //We pretend that pseq starts at this index
    //Connect leaves
    for(int i = 0; i < n; i++) {
      adjacencies[pseq[firstindex]-n] += i;
      V = V - i;
      firstindex++;
    }//END add leaves
    
    dbgtrace << "Connecting edges" << endl;
    dbgtrace << "V: " << V << endl;
    dbgtrace << "Adjacencies: " << adjacencies << endl;
    
    //Now create edges
    for(int i = 1; i <= n-3; i++) {
      //If there are only two vertices left, connect them
      if(i == n-3) {
	Vector<int> lasttwo(V);
	result |= adjacencies[lasttwo[0]-n];
      }
      else {
	//Find the minimal element in V that is not in the sequence (starting at firstindex)
	Set<int> pset(pseq.slice(~sequence(0,firstindex)));
	int smallest = -1;
	for(Entire<Set<int> >::iterator vit = entire(V); !vit.at_end(); vit++) {
	  if(!pset.contains(*vit)) {
	      smallest = *vit;break;
	  }
	}//END look for smallest in V\P
	Set<int> Av = adjacencies[smallest-n];
	result |= Av;
	adjacencies[pseq[firstindex]-n] += Av;
	V = V - smallest;
	firstindex++;
      }
    }//END create edges
    
    return result;
  }//END decodePrueferSequence
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object local_mn(std::vector<perl::Object> curves) {
    
    //Compute the set of all vertex valences > 3
    Set<int> valences;
    for(unsigned int c = 0; c < curves.size(); c++) {
	Vector<int> degrees = curves[c].give("NODE_DEGREES");
	valences += (Set<int>(degrees) - 3);
    }
    dbgtrace << "Valences: " << valences << endl;
    
    dbgtrace << "Computing Prüfer sequences" << endl;
    //Now we compute combinatorially all M_0,v for all possible valences v
    Map<int, Vector<Vector<Set<int> > > > combinatorial_mns;
    for(Entire<Set<int> >::iterator v = entire(valences); !v.at_end();v++) {
      Vector<Vector<int> > pseq = computePrueferSequences(*v);
      Vector<Vector<Set<int> > > cmns;
      for(int p = 0; p < pseq.dim(); p++) cmns |= decodePrueferSequence(pseq[p]);
      combinatorial_mns[*v] = cmns;
    }
    dbgtrace << "Done." << endl;
    //Now iterate through all curves and compute their adjacent maximal cones
    
    Vector<Set<int> > rays; //Ray list in terms of partitions
    Vector<Set<int> > cones; //Cone list in terms of ray indices
    Vector<Set<int> > local_cones; 
    
    int n_leaves = curves[0].give("N_LEAVES");
    Set<int> all_leaves = sequence(1,n_leaves);
    
    //First we add all the curve rays and construct the corresponding local cones
    for(unsigned int cu = 0; cu < curves.size(); cu++) {
      IncidenceMatrix<> set_list = curves[cu].give("SETS");
      Set<int> l_cone;
      //For each ray, check if it exists already
      for(int cu_ray = 0; cu_ray < set_list.rows(); cu_ray++) {
	//Normalize
	Set<int> cu_set = set_list.row(cu_ray);
	if(cu_set.contains(n_leaves)) cu_set = all_leaves - cu_set;
	int ray_index = -1;
	for(int oray = 0; oray < rays.dim(); oray++) {
	  if(rays[oray] == cu_set) {
	    ray_index = oray; break;
	  }
	}
	if(ray_index == -1) {
	    rays |= cu_set;
	    ray_index = rays.dim()-1;
	}
	l_cone += ray_index;
      }
      local_cones |= l_cone;
    }//END create curve rays and local cones
    

    
    //Then we construct the actual cones
    for(unsigned int cu = 0; cu < curves.size(); cu++) {
      dbgtrace << "Computing on curve " << cu+1 << endl;
      //We iteratively compute the cartesian product of the M_0,ns at the different vertices
      Vector<Set<int> > cones_so_far; 
	cones_so_far |= local_cones[cu];
      IncidenceMatrix<> nodes_by_sets = curves[cu].give("NODES_BY_SETS");
      IncidenceMatrix<> nodes_by_leaves = curves[cu].give("NODES_BY_LEAVES");
      IncidenceMatrix<> set_list = curves[cu].give("SETS");
      Vector<int> degrees = curves[cu].give("NODE_DEGREES");
      for(int node = 0; node < nodes_by_sets.rows(); node++) {
	dbgtrace << "Computing on node " << node << endl;
	if(degrees[node] > 3) {
	  //We have to translate the cones of the M_0,val(node) first
	  Vector<Set<int> > translated_cones;
	  Vector<Vector<Set<int> > > valence_cones = combinatorial_mns[degrees[node]];
	  //The semantics are the following: The edges adjacent to node are numbered
	  //in this order: First the leaves ordered from lowest to highest number, then the 
	  //bounded edges, ordered according to the order of their appearance in NODES_BY_SETS
	  Vector<Set<int> > adjacent_edges;
	    Set<int> leaves = nodes_by_leaves.row(node);
	    for(Entire<Set<int> >::iterator l = entire(leaves); !l.at_end(); l++) {
	      Set<int> singleset; singleset += *l;
	      adjacent_edges |= singleset;
	    }
	    Vector<int> sets(nodes_by_sets.row(node));
	    for(int s = 0; s < sets.dim(); s++) {
	      //We normalize the sets to 'point away' from the node, i.e. they should have
	      //empty intersection with the first leaf (or already normalized set)
	      //If there is no leaf or normalized set, the intersection with the next set 
	      //must be either emtpy or this set. Otherwise we take the complement.
	      Set<int> compare_set;
	      if(adjacent_edges.dim() > 0) compare_set = adjacent_edges[0];
	      else compare_set = set_list.row(sets[1]);
	      int isize = (set_list.row(sets[s]) * compare_set).size();
	      if(isize == 0 || isize == set_list.row(sets[s]).size()) {
		adjacent_edges |= set_list.row(sets[s]);
	      }
	      else {
		adjacent_edges |= all_leaves - set_list.row(sets[s]);
	      }
	    }
	  //Now go through the cones and translate
	  Vector<Set<int> > new_cones;
	  for(int vc = 0; vc < valence_cones.dim(); vc++) {
	    dbgtrace << "Computing on cone " << vc << endl;
	    Set<int> ray_indices;
	    //Go through all rays of the cone
	    for(int r = 0; r < valence_cones[vc].dim(); r++) {
	      //Translate the ray
	      dbgtrace << "Adjacent: " << adjacent_edges << endl;
	      dbgtrace << "Cone: " << valence_cones[vc] << endl;
	      Set<int> ray_partition = accumulate(adjacent_edges.slice(valence_cones[vc][r]),operations::add());
	      //To avoid doubles via complements, we make sure that N_LEAVES is always in the
	      //complement
	      if( ray_partition.contains(n_leaves)) ray_partition = all_leaves - ray_partition;
	      dbgtrace << "Ray is " << ray_partition << endl;
	      //Check if this ray already exists
	      int ray_index = -1;
	      for(int oray = 0; oray < rays.dim(); oray++) {
		if(rays[oray] == ray_partition) {
		  ray_index = oray; break;
		}
	      }
	      if(ray_index == -1) {
		rays |= ray_partition;
		ray_index = rays.dim()-1;
	      }
	      dbgtrace << "Index is: " << ray_index << endl;
	      ray_indices += ray_index;
	    }//END iterate cone rays
	    new_cones |= ray_indices;
	  }//END iterate local mn cones
	  dbgtrace << "New cones: " << new_cones << endl;
	  //Finally take the cartesian products of the new cones with all the old ones
	  Vector<Set<int> > replace_cones_so_far;
	  for(int nc = 0; nc < new_cones.dim(); nc++) {
	    for(int oc = 0; oc < cones_so_far.dim(); oc++) {
	      replace_cones_so_far |= (new_cones[nc] + cones_so_far[oc]);
	    }
	  }
	  cones_so_far = replace_cones_so_far;
	}//END if deg(v) > 3
      }//END iterate nodes
      //Now we have to check for cones that already exist
      Set<int> double_cones;
      for(int nc = 0; nc < cones_so_far.dim(); nc++) {
	for(int oc = 0; oc < cones.dim(); oc++) {
	    if(cones[oc] == cones_so_far[nc]) double_cones += nc;
	}
      }
      cones |= cones_so_far.slice(~double_cones);
//       dbgtrace << "Cones: " << cones << endl;
//       dbgtrace << "Rays: " << rays << endl;
    }//END iterate curves
    
    //Finally we convert the rays to matroid coordinates
    Matrix<int> E = pair_index_map(n_leaves-1);
    Matrix<Rational> bergman_rays(rays.dim(),(n_leaves*(n_leaves-1))/2 - n_leaves);
    Vector<Rational> onlyones = ones_vector<Rational>(bergman_rays.cols());
    for(int r = 0; r < rays.dim(); r++) {
      Vector<int> raylist(rays[r]);
//       dbgtrace << "Converting ray " << raylist << endl;
      Vector<Rational> newray(bergman_rays.cols());
      for(int k = 0; k < raylist.dim()-1; k++) {
	  for(int l = k+1; l < raylist.dim(); l++) {
	    int newrayindex = E(raylist[k]-1,raylist[l]-1);
// 	    dbgtrace << "nri: " << newrayindex << endl;
	    //If the newrayindex is one higher than the ray dimension, 
	    //this means it is first of all the last pair. Also, we don't
	    //add -e_n but e_1 + ... + e_{n-1} (as we mod out lineality)
	    if(newrayindex < newray.dim()) {
		newray[newrayindex] = -1;
	    }
	    else {
		newray = newray + onlyones;
	    }
	  }
      }
//       dbgtrace << "Result: " << newray << endl;
      bergman_rays.row(r) = newray; 
    }
    
    Vector<Integer> weights = ones_vector<Integer>(cones.dim());
    
    dbgtrace << "Rays " << rays << endl;
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << bergman_rays;
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("LOCAL_RESTRICTION") << local_cones;
      return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry / Moduli spaces / Local geometry"
		    "# Computes the moduli space M_0,n locally around a given list of combinatorial"
		    "# types. More precisely: It computes the weighted complex consisting of all"
		    "# maximal cones containing any of the given combinatorial types and localizes "
		    "# at these types "
		    "# This should only be used for curves of small codimension. What the function "
		    "# actually does, is that it combinatorially computes the cartesian products "
		    "# of M_0,v's, where v runs over the possible valences of vertices in the curves"
		    "# For max(v) <= 8 this should terminate in a reasonable time (depending on the "
		    "# number of curves)"
		    "# The coordinates are the same that would be produced by the function "
		    "# tropical_m0n"
		    "# @param Int n Should be >= 3"
		    "# @return WeightedComplex The local complex",
		    &local_mn,"local_m0n(;@)");
  
  Function4perl(&decodePrueferSequence,"dcp(Vector<Int>)");
}}