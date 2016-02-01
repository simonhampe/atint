/*
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA  02110-1301, USA.

   ---
   Copyright (C) 2011 - 2015, Simon Hampe <simon.hampe@googlemail.com>

   Computes the lattice of chains of a lattice.
   */

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Map.h"
#include "polymake/graph/HasseDiagram.h"
#include "polymake/tropical/cyclic_chains.h"

namespace polymake { namespace tropical {

   using polymake::graph::HasseDiagram;

   /*
    * @brief Given a chain of sets in terms of rows of an incidence matrix and another set S, also
    * as a row index, determines whether that set S fits into the existing chain to form a larger chain.
    */
   bool is_chain(const IncidenceMatrix<> &faces, const Set<int> &chain, int additional_face) {
      Set<int> face_set = faces.row(additional_face);
      for(Entire<Rows<pm::MatrixMinor<const IncidenceMatrix<>, const Set<int>&, const pm::all_selector& > > >::iterator cs = entire(rows(faces.minor(chain,All))); 
            !cs.at_end(); cs++){
         Set<int> inter = (*cs) * face_set;
         if(inter.size() != (*cs).size() && inter.size() != face_set.size()) return false;
      }
      return true;
   }

   /*
    * Takes a face lattice and computes the lattice of all chains, which contain the
    * top and bottom node. This adds an artificial top node. 
    */
   HasseDiagram chain_lattice(perl::Object facelattice) {

      IncidenceMatrix<> faces = facelattice.give("FACES");
      int top_node = facelattice.CallPolymakeMethod("top_node");
      int bottom_node = facelattice.CallPolymakeMethod("bottom_node");

      HasseDiagram HD;
      HasseDiagram::_filler chain_filler = filler(HD,true);

      //The bottom node consists of top/bottom node of the original lattice 
      Set<int> bottom_set = scalar2set(top_node) + scalar2set(bottom_node);
      int new_bottom_index = chain_filler.add_node(bottom_set);

      //A list of nodes of the currently maximal dimension (i.e. size)
      // We add to these nodes in the next step
      Set<int> max_dim_nodes = scalar2set(new_bottom_index);
      //A list of nodes which are maximal, i.e. where we cant add any sets
      Set<int> final_nodes;

      while(max_dim_nodes.size() > 0) {
         chain_filler.increase_dim();
         Set<int> new_max_dim_nodes;
         Map<Set<int>, int> new_face_indices; 
         for(Entire<Set<int> >::iterator mdn = entire(max_dim_nodes); !mdn.at_end(); mdn++) {
            bool added_some_set = false;
            Set<int> mdn_face = HD.face(*mdn);

            Set<int> flats_to_test = sequence(0, faces.rows()) - mdn_face;
            for(Entire<Set<int> >::iterator ftt = entire(flats_to_test);  !ftt.at_end(); ftt++) {
               if(is_chain(faces, mdn_face, *ftt)) {
                  Set<int> newset = mdn_face + (*ftt);
                  Map<Set<int>,int>::iterator ns_index_it = new_face_indices.find(newset);
                  int ns_index;
                  //If its a new node, add it
                  if(ns_index_it == new_face_indices.end()){
                     ns_index = chain_filler.add_node( newset);
                     new_face_indices[newset] = ns_index;
                     new_max_dim_nodes += ns_index;
                  }
                  else {
                     ns_index = (*ns_index_it).second;
                  }
                  //In any case, connect it
                  chain_filler.add_edge( *mdn, ns_index);
                  added_some_set = true;
               }
            }//END iterate over flats to test.

            if(!added_some_set) {
               final_nodes += *mdn;
            }
         }//END iterate over maximal dim nodes

         max_dim_nodes = new_max_dim_nodes;
      }

      //Add the final, artificial top node
      int final_node = chain_filler.add_node(Set<int>());
      for(Entire<Set<int> >::iterator fn = entire(final_nodes); !fn.at_end(); fn++) {
         chain_filler.add_edge(*fn, final_node);
      }

      return HD;
   }


   //Computes the set of indices of nodes that lie above a given node
   Set<int> nodes_above(const HasseDiagram& HD, int node) {
      Set<int> result(HD.out_adjacent_nodes( node));
      std::list<int> queue;
      for(Entire<Set<int> >::iterator oa = entire(result); !oa.at_end(); oa++) queue.push_back(*oa);

      while(queue.size() > 0) {
         int next = queue.front();
         queue.pop_front();
         Set<int> nbrs = HD.out_adjacent_nodes( next);
         for(Entire<Set<int> >::iterator adj = entire(nbrs); !adj.at_end(); adj++) {
            result += *adj;
            queue.push_back(*adj);
         }
      }
      return result;
   }

   /*
    * @brief Takes a Hasse diagram and computes for each node n the value of the moebius function
    * mu(n,1), where 1 is the maximal element.
    * @return Vector<int> Each entry corresponds to the node of the same index.
    */
   Vector<int> top_moebius_function(HasseDiagram HD) {
      Vector<int> result(HD.nodes());

      result[ HD.top_node() ] = 1;
      int HD_dim = HD.dim();

      for(int r = HD_dim-1; r >= 0; r--) {
         Set<int> n_of_dim = HD.nodes_of_dim(r);
         for(Entire<Set<int> >::iterator nr = entire(n_of_dim); !nr.at_end(); nr++) {
            int value = 0;
            Set<int> above = nodes_above(HD, *nr);
            for(Entire<Set<int> >::iterator ab = entire(above); !ab.at_end(); ab++) {
               value -= result[ *ab];
            }
            result[*nr] = value;
         }
      }

      //The bottom node needs to be set manually - its minus the sum over all other values
      int bottom_value = accumulate( result, operations::add());
      result[ HD.bottom_node()] = -bottom_value;

      return result;
   }


   /*perl::Object chain_lattice_wrapper(perl::Object fl) {
      HasseDiagram HD = chain_lattice(fl);
      perl::Object result("FaceLattice");
         result.take("ADJACENCY") << HD.graph();
         result.take("FACES") << HD.faces();
         result.take("DIMS") << HD.dims();
      cout << top_moebius_function(HD) << endl;
      return result;
   }


   UserFunction4perl("", &chain_lattice_wrapper, "chain_lattice(FaceLattice)");*/

}}
