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
 * This file computes matroid fans via their flats (i.e. very inefficiently)
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/Map.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the set of flats and all maximal chains of flats of a matroid
   @param perl::Object matroid A matroid::Matroid object
   @param Vector<Vector<Set<int> > > flats A reference. Will be filled with the flats. At position i = 0,..,rank(matroid)-1 it contains a list of all flats with rank (i+1) (i.e. the empty set and the complete set are excluded)
   @param Vector<Vector<int> > chains A reference. Will be filled with the maximal chains of flats: Each element in the list is again a list of indices. chains[i][j] is the index of the rank-(j+1)-flat in flats[j] for the i-th chain.
   */
  void chains_of_flats(perl::Object matroid, Vector<Vector<Set<int> > > &flats, Vector<Vector<int> > &chains) {
    //First we need to compute all flats of the matroid
    Array<Set<int> > bases = matroid.give("BASES");
    int rank = matroid.give("RANK");
    int n = matroid.give("N_ELEMENTS");
    
    //Contains at position i the rank i+1-flats
    flats = Vector<Vector<Set<int> > >(rank-1); 
    
    //Go through all subsets
    Array<Set<int> > ssets= all_subsets(sequence(0,n));
    for(int s = 0; s < ssets.size(); s++) {
      if(ssets[s].size() == 0 || ssets[s].size() == n) continue;
      //Compute the rank of the subset and remember the union of the bases that intersect the set maximally
      int srank = 0;
      Set<int> maximal_bases_union;
      for(int b = 0; b < bases.size(); b++) {
	int ib = (bases[b] * ssets[s]).size();
	if(ib == srank) maximal_bases_union += bases[b];
	if(ib > srank) {
	  srank = ib;
	  maximal_bases_union = bases[b];
	}
      }//END compute rank
      //Now check if all elements in the complement increase the rank, i.e. are contained in one of the
      //bases that intersect the set maximally
      Set<int> complement = sequence(0,n) - ssets[s];
      if( (complement * maximal_bases_union).size() == complement.size()) {
	flats[srank-1] |= ssets[s];
      }
    }//END compute flats
    
    //Now we need to compute maximal chains
    chains = Vector<Vector<int> >() ;
    //Start with rank-1-flats
    for(int i = 0; i < flats[0].dim(); i++) {
      Vector<int> oneflat; oneflat |= i;
      chains |= oneflat;
    }
    
    //Now subsequently go through rank k - flats and see if you can attach them to a lower chain
    for(int k = 1; k < flats.dim(); k++) {
      Vector<Vector<int> > newchains;
      for(int f = 0; f < flats[k].dim(); f++) {
	//Go through all chains
	for(int c = 0; c < chains.dim(); c++) {
	  Set<int> lastflat = (flats[k-1])[ (chains[c])[chains[c].dim()-1] ];
	  if( (lastflat * (flats[k])[f]).size() == lastflat.size()) {
	    Vector<int> attached_chain =  chains[c];
	    attached_chain |= f;
	    newchains |= attached_chain;
	  }
	}//END iterate rank k-1 chains
      }//END iterate rank k flats
      chains = newchains;
    }//END iterate ranks
  }//END chains_of_flats
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object bergman_by_flats(perl::Object matroid) {
    int rank = matroid.give("RANK");
    int n = matroid.give("N_ELEMENTS");
    
    //Compute flats and chains of flats
    Vector<Vector<Set<int> > > flats(rank-1); 
    Vector<Vector<int> > chains;
    chains_of_flats(matroid, flats,chains);
    
    //Now create the matroid fan
    Matrix<Rational> rays(0,n-1);
    //Go through all flats and create appropriate rays
    for(int r = 0; r < flats.dim(); r++) {
      for(int f = 0; f < flats[r].dim(); f++) {
	//The last element, n-1 is modded out as sum of the rest
	Set<int> flat = (flats[r])[f];
	Vector<Rational> matray(n-1);
	if( flat.contains(n-1) ) {
	  matray.slice(sequence(0,n-1) - flat) = ones_vector<Rational>(n - flat.size());
	}
	else {
	  matray.slice(flat) = (- ones_vector<Rational>(flat.size()));
	}
	rays /= matray;
      }
    }//END create rays
    
    //Now create cones from the list of maximal chains
    Vector<Set<int> > cones;
    for(int c = 0; c < chains.dim(); c++) {
      Set<int> chainset;
      int offset = 0;
      for(int f = 0; f < chains[c].dim(); f++) {
	chainset += ( (chains[c])[f] + offset);
	offset += flats[f].dim();
      }
      cones |= chainset;
    }
    
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << ones_vector<Integer>(cones.dim());
      
    return result;
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
//   //Documentation see perl wrapper
//   perl::Object matroid_intersection_by_flats(perl::Object X, perl::Object Y, perl::Object matroid) {
//     //Extract values
//     int Xdim = X.give("CMPLX_DIM");
//     int Ydim = Y.give("CMPLX_DIM");
//     int n = matroid.give("N_ELEMENTS");
//     Array<Set<int> > bases = matroid.give("BASES");
//     
//     //dbgtrace << "Checking codimension" << endl;
//     
//     //If the codimensions of the varieties add up to something larger then CMPLX_AMBIENT_DIM, return the 0-cycle 
//     if(Xdim + Ydim - (n-1) < 0) {
//       return CallPolymakeFunction("zero_cycle");
//     }
//     
//     //Better: Compute flats of matroid, take sums and compute corr. bergman fan from there together
//     //with rational functions
//     
//     
// //     //Compute the sum matroid + matroid
// //     Vector<Set<int> > bases_pairs;
// //     Map<int,int> shift_map;
// //       for(int i = 0; i < n; i++) {
// // 	shift_map[i] = i+n;
// //       }
// //     for(int b = 0; b< bases.size(); b++) {
// //       for(int c = 0; c < bases.size(); c++) {
// // 	Set<int> shifted_set = attach_operation(bases[c], pm::operations::associative_access<Map<int,int>, int>(&shift_map));
// // 	bases_pairs |= (bases[b] * shifted_set);
// //       }
// //     }
// //     perl::Object sum_matroid("matroid::Matroid");
// //       sum_matroid.take("N_ELEMENTS") << 2*n;
// //       sum_matroid.take("BASES") << bases_pairs;
// //       
// //     //Compute the corresponding bergman fan
// //     perl::Object bergman_fan = bergman_by_flats(sum_matroid);
// //     
// //     //Extract the rays and construct the corresponding diagonal functions
//     
//     
//     
//     
//   }

  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry / Matroid fans"
		    "# This function computes the bergman fan of a matroid in the subdivision associated"
		    "# to its lattice of flats. In particular it computes all flats and all maximal"
		    "# chains of flats. Hence this function is terribly inefficient and should only be used"
		    "# with comparatively small matroids"
		    "# @param matroid::Matroid m Any matroid"
		    "# @return WeightedComplex The bergman fan B(m) in its flat subdivision",
		    &bergman_by_flats,"bergman_fan_flats(matroid::Matroid)");
  
}}