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
#include "polymake/graph/maximal_chains.h"
#include "polymake/graph/hasse_diagram_complex.h"
#include "polymake/tropical/cyclic_chains.h"
#include "polymake/Bitset.h"
#include "polymake/PowerSet.h"

namespace polymake { namespace tropical {

   //Computes the set of indices of nodes that lie above a given node
   Bitset nodes_above(const Lattice<BasicDecoration>& HD, int node) {
      Bitset result(HD.out_adjacent_nodes( node));
      std::list<int> queue;
      for(const auto& oa : result) queue.push_back(oa);

      while(queue.size() > 0) {
         int next = queue.front();
         queue.pop_front();
         Set<int> nbrs = HD.out_adjacent_nodes( next);
         for(auto &adj : nbrs) {
            result += adj;
            queue.push_back(adj);
         }
      }
      return result;
   }

   /*
    * @brief Takes a Hasse diagram and computes for each node n the value of the moebius function
    * mu(n,1), where 1 is the maximal element.
    * @return Vector<int> Each entry corresponds to the node of the same index.
    */
   Vector<int> top_moebius_function(const Lattice<BasicDecoration>& HD) {
      Vector<int> result(HD.nodes());

      result[ HD.top_node() ] = 1;
      int HD_dim = HD.rank();

      for(int r = HD_dim-1; r >= 0; r--) {
         auto n_of_dim = HD.nodes_of_rank(r);
         for(auto& nr : n_of_dim) {
            int value = 0;
            Bitset above = nodes_above(HD, nr);
            for(const auto& ab : above) {
               value -= result[ ab];
            }
            result[nr] = value;
         }
      }

      //The bottom node needs to be set manually - its minus the sum over all other values
      int bottom_value = accumulate( result, operations::add());
      result[ HD.bottom_node()] = -bottom_value;

      return result;
   }

   //Computes the lattice of chains which contain bottom and top node
   Lattice<BasicDecoration> cyclic_chains(const Lattice<BasicDecoration>& lattice) {
      IncidenceMatrix<> max_chains(maximal_chains(lattice,false,false));
      const int bottom_index = lattice.bottom_node();
      const int top_index = lattice.top_node();
      Array<IncidenceMatrix<> > max_coatoms(max_chains.rows());
      Array<int> maximal_dims(max_chains.rows());
      if(lattice.graph().nodes() <= 2){
         return graph::hasse_diagram_general(max_chains, max_coatoms, 0, maximal_dims, graph::lattice::RankRestriction(),
               graph::lattice::TopologicalType(true,true), Set<int>());
      }
      auto mv_it = entire(max_coatoms);
      auto mc_it = entire(rows(max_chains));
      auto md_it = entire(maximal_dims);
      int dim = 0;
      const Set<int> extreme_nodes = scalar2set(bottom_index) + scalar2set(top_index);
      for(; !mc_it.at_end(); ++mc_it, ++mv_it, ++md_it) {
         int size = mc_it->size();
         *md_it = size-1;
         dim = std::max(dim, size-1);
         *mv_it = IncidenceMatrix<>(size-2,lattice.graph().nodes());
         Set<int> non_extremal = *mc_it - extreme_nodes;
         auto removals = entire(all_subsets_of_k( non_extremal, size-3));
         for(auto inc_rows = entire(rows(*mv_it)); !inc_rows.at_end(); ++inc_rows, ++removals) {
            *inc_rows = extreme_nodes + *removals;
         }
      }
      return graph::hasse_diagram_general(max_chains, max_coatoms, dim, maximal_dims, graph::lattice::RankRestriction(), graph::lattice::TopologicalType(), Set<int>());

   }


}}
