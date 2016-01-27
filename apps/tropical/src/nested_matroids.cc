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

   Computations on nested matroids
   */

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/Map.h"
#include "polymake/graph/HasseDiagram.h"
#include "polymake/tropical/cyclic_chains.h"

namespace polymake { namespace tropical {

   using polymake::graph::HasseDiagram;

   /*
    * @brief Computes the maximal transversal presentation of a nested matroid from its chain of
    * cyclic flats
    * @param int n The size of the ground set
    * @param Array<Set<int> > flats The cyclic flats, ordered from smallest to largest
    * @param Array<int> ranks The ranks of the corresponding flats
    */
   IncidenceMatrix<> presentation_from_chain(int n, const IncidenceMatrix<>& flats, const Array<int> ranks) {
      Set<int> coloops = sequence(0, n) - flats[flats.rows()-1];
      int total_rank = coloops.size() + ranks[ranks.size()-1];
      IncidenceMatrix<> result(total_rank, n);

      //First: coloops as complements of largest cyclic flat
      int current_index = 0;
      for(int i = 0; i < coloops.size(); i++,current_index++) {
         result[i] = coloops;
      }

      //Move backwards in the list of flats
      for(int j =  flats.rows()-2; j >= 0; j--) {
         Set<int> complement = sequence(0,n) - flats[j];
         int occ = ranks[j+1] - ranks[j];
         for(int k = 0; k < occ; k++, current_index++) {
            result[current_index] = complement;
         }
      }

      return result;
   }


   //Compute a representation of a matroid in the basis of nested matroids
   //in the matroid intersection ring.
   perl::ListReturn matroid_nested_decomposition(perl::Object matroid) {
      int n = matroid.give("N_ELEMENTS");
      perl::Object flats = matroid.give("LATTICE_OF_CYCLIC_FLATS");
      IncidenceMatrix<> cyclic_flats = flats.give("FACES");
      //Determine orientation of flats 
      bool flipped = cyclic_flats.row(0).size() > cyclic_flats.row(cyclic_flats.rows()-1).size();

      HasseDiagram chains = chain_lattice(flats); 
      IncidenceMatrix<> chain_faces = chains.faces();
      Vector<int> coefficients = top_moebius_function(chains);
      Set<int> supp = support( coefficients);

      Array<IncidenceMatrix<> > nested_presentations(supp.size());
      Array<int> final_coefficients(supp.size());

      int current_index = 0;
      for(Entire<Set<int> >::iterator s = entire(supp); !s.at_end(); s++, current_index++) {
         final_coefficients[current_index] = coefficients[*s];
         IncidenceMatrix<> current_faces = cyclic_flats.minor(chains.face( *s),All);
         IncidenceMatrix<> ordered_faces = flipped? 
            IncidenceMatrix<>(current_faces.rows(), current_faces.cols(), rentire(rows(current_faces))) : 
            current_faces;
         


      }

      

      perl::ListReturn result;
         result << nested_presentations;
         result << final_coefficients;

      return result;
   }



   Function4perl(&presentation_from_chain, "presentation_from_chain($, $,$)");

}}
