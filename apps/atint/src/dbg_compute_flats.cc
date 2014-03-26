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
#include "polymake/Set.h"
#include "polymake/linalg.h"
#include "polymake/Map.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
//   //using namespace atintlog::dolog;
//   //using namespace atintlog::dotrace;
//   
//   ///////////////////////////////////////////////////////////////////////////////////////

     Vector<Set<int> > compute_fibre_bases(perl::Object matroid, Set<int> E) {
      
       //Extract bases
       Array<Set<int> > m_bases = matroid.give("BASES");
       int N = matroid.give("N_ELEMENTS");
       
       //Determine rank of E and E-maximal bases
       int erank = 0;
       Set<int> max_basis_indices;
       for(int b = 0; b < m_bases.size(); b++) {
	  int sz = (m_bases[b]*E).size();
	  if(sz == erank) max_basis_indices += b;
	  if(sz > erank) {
	    erank = sz;
	    max_basis_indices = Set<int>();
	    max_basis_indices += b;
	  }
       }
       
       //Create maps lambda_i (we attach the second copy of R at the end)
       Set<int> R = sequence(0,N) - E;
       Map<int,int> second_r_map;
       int current_index = N;
       for(Entire<Set<int> >::iterator r = entire(R); !r.at_end(); r++) {
	  second_r_map[*r] = current_index;
	  current_index++;
       }
       for(Entire<Set<int> >::iterator e = entire(E); !e.at_end(); e++) {
	  second_r_map[*e] = *e;
       }
       
       //Create bases
       Set<Set<int> > result;
       for(int b1 = 0; b1 < m_bases.size(); b1++) {
	  for(Entire<Set<int> >::iterator b2 = entire(max_basis_indices); !b2.at_end(); b2++) {
		Set<int> pairset1 = m_bases[b1] + Set<int>(attach_operation(m_bases[*b2] * R,								   pm::operations::associative_access<Map<int,int> ,int >(&second_r_map)));
		Set<int> pairset2 = R * m_bases[*b2] + Set<int>(
		  attach_operation(m_bases[b1],					   pm::operations::associative_access<Map<int,int>,int >(&second_r_map) ) );
		result += pairset1;
// 		result += pairset2;	    
	  }
       }
       
       return Vector<Set<int> >(result);
       
       
       
     }

//   
//   void computeFlats(Matrix<Rational> m) {
//     int r = rank(m);
//     
//     Vector<Vector<Set<int> > > flats(r);
//     
//     Set<int> all_cols = sequence(0,m.cols());
//     Array<Set<int> > ssets = all_subsets(all_cols);    
//     for(int s = 0; s < ssets.size(); s++) {
//       int x = rank(m.minor(All,ssets[s]));
//       if(x == 0) continue;
//       //Check if s is a flat
//       Set<int> complement = all_cols - ssets[s];
//       bool found_bad_element = false;
//       for(Entire<Set<int> >::iterator c = entire(complement); !c.at_end(); c++) {
// 	if(rank(m.minor(All,ssets[s] + *c)) == x) {
// 	  found_bad_element = true; break;
// 	}
//       }
//       if(!found_bad_element) {
// 	flats[x-1] |= ssets[s];
//       }
//     }
//     
//     pm::cout << "Done: " << endl;
//     
//     //Output
//     for(int i = 0; i < flats.size(); i++) {
//       pm::cout << "Rank " << (i+1) << ": " << flats[i] << endl;
//     }
//     
//   }
/*
  perl::Object bergman_by_flats(perl::Object matroid) {
    //First we need to compute all flats of the matroid
    Array<Set<int> > bases = matroid.give("BASES");
    int rank = matroid.give("RANK");
    int n = matroid.give("N_ELEMENTS");
    
    //Contains at position i the rank i+1-flats
    Vector<Vector<Set<int> > > flats(rank-1); 
    
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
    Vector<Vector<int> > chains;
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
    
  }*/

  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
//   Function4perl(&computeFlats, "computeFlats(Matrix<Rational>)");
//   Function4perl(&bergman_by_flats,"bbf(matroid::Matroid)");
  
}}