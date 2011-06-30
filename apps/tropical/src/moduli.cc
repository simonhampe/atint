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
Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>

This file provides the functionality necessary to compute tropical moduli spaces
*/

#include "polymake/client.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/Vector.h"

namespace polymake { namespace tropical {

  using namespace atint::donotlog;
  //using namespace atint::dolog;
  //using namespace atint::dotrace;
  
//Documentation see perl wrapper
perl::Object tropical_mn(int n) {
  if(n < 4) {
    throw std::runtime_error("Number of leaves should be at least 4 for M_0,n computation");
  }
  
  //First we compute the set of all ordered paired Prüfer sequences on n,...,2n-3
  // (i.e. all sequences of length 2n-4 where each element from n,...,2n-3 occurs twice
  // and where removing every second occurence of an element yields an ascending sequence)
  Vector<Vector<int> > sequences;
  
  //Compute the number of sequences = number of maximal cones
  int noOfMax = 1;
  for(int i = 0; i <= n-4; i++) {
    noOfMax *= (2*(n-i) -5);
  }
  
  //Iterate through all Prüfer sequences
  Vector<int> indices = ones_vector<int>(n-2);
  for(int iteration = 0; iteration < noOfMax; iteration++) {
    //Create the sequence currently represented by indices and append it
    Vector<int> baseSequence = zero_vector<int>(2*n -4);
    for(int i = 0; i < n-1; i++) {
      //Go through the non-zero entries of baseSequence. If it is the first or the indices[i]+1-th, 
      //insert an n+i
      int nonzero_count = -1;
      for(int entry = 0; entry < baseSequence.dim(); entry++) {
	if(baseSequence[entry] == 0) {
	    nonzero_count++;
	    if(nonzero_count == 0) {
	      baseSequence[entry] = n+i;
	    }
	    if(nonzero_count == indices[i]) {
	      baseSequence[entry] = n+i;
	      break;
	    }
	}
      }
    }
    sequences = sequences | baseSequence;
    //Increase the indices vector by "1"
    if(iteration < noOfMax-1) {
      int counterindex = n-3;
      while(indices[counterindex] == 2*(n-counterindex)-5) {
	indices[counterindex] = 1;
	counterindex--;
      }
      indices[counterindex]++;
    }
  }
  
  cout << "Found" << sequences.dim() << " prüfer sequences" << endl;
  
  return perl::Object("fan::PolyhedralFan");
  
}

UserFunction4perl("# @category Tropical geometry"
		  "# Creates the moduli space of abstract rational n-marked curves. Its coordinates are"
		  "# given as the coordinates of the bergman fan of the matroid of the complete graph on "
		  "# n-1 nodes (but not computed as such)"
		  "# @param Int n The number of leaves. Should be at least 4"
		  "# @return fan::PolyhedralFan The tropical moduli space M_0,n",
		  &tropical_mn, "tropical_m0n($)");
  
}}
