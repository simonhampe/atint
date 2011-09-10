/*
 T his *program is free software; you can redistribute it and/or
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
 
 Contains a function to create the basis of the matroid of the complete graph 
 on n elements
 */

#include "polymake/client.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/Set.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"


namespace polymake { namespace atint {
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Computes the set of all spanning trees of the complete graph on n nodes
    @param int n The number of nodes
    @return An array of Set<Int> The spanning trees, each encoded as the set of edges, where the edge indices are determined as follows: For i,j = 0,..,n-1 and i < j we sort the edges (i,j) in ascending lexicographical ordering, i.e, (0,1) < (0,2) < .. < (n-2,n-1).
  */
  perl::ListReturn spanning_complete(int n) {
    //First we create the edge index matrix E(i,j) that contains at element i,j the edge index of edge (i,j)
    int nextindex = 0;
    Matrix<int> E(n,n);
    for(int i = 0; i < n; i++) {
      for(int j = i+1; j < n; j++) {
	E(j,i) = E(i,j) = nextindex;
	nextindex++;
      }
    }
    
    dbgtrace << "Computed edge index matrix" << endl;
    
    //Now we generate the set of all Prüfer sequences, i.e. all sequences of length n-2 with elements from 0 to n-1 (repetitions allowed)
    Vector<Vector<int> > prseq;
    Vector<int> counter = zero_vector<int>(n-2);
    prseq = prseq | Vector<int>(counter);
    int iterations = pow(n,n-2)-1;
    for(int i = 0; i < iterations; i++) {
      int counterindex = n-3;
      while(counter[counterindex] == n-1) {
	counter[counterindex] = 0;
	counterindex--;
      }
      counter[counterindex]++;
      prseq = prseq | Vector<int>(counter);
    }
    
    dbgtrace << "Computed all prüfer sequences " << prseq << endl;
    
    //Now we decode each sequence into a tree and save that tree
    perl::ListReturn result;
    Set<int> completeSet = sequence(0,n);
    for(int seqindex = 0; seqindex < prseq.dim(); seqindex++) {
      Vector<int> seq = prseq[seqindex];
      Set<int> edgeset;
      //Initialize V as 0,..,n-1
      Set<int> V(completeSet);
      for(int i = 0; i < n-2; i++) {
	Set<int> pset(seq);
	//Find the smallest element in V that is not in seq
	int smallest = 0;
	for(Entire<Set<int> >::iterator vit = entire(V); !vit.at_end(); vit++) {
	  if(!(pset.contains(*vit))) {
	    smallest = *vit;
	    break;
	  }
	}
	//Create edge [smallest, seq[0]]
	edgeset += E(smallest,seq[0]);
	V = V - smallest;
	seq = seq.slice(1,seq.dim()-1);
      }
      //Finally create an edge from the last two elements left
      Vector<int> last(V);
      edgeset += E(last[0],last[1]);
      result << edgeset;
    }
    
    dbgtrace << "Computed all spanning trees" << endl;
    
    return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Graph theory"
		    "# Computes the set of all spanning trees of the complete graph on n nodes"
		    "# @param int n The number of nodes"
		    "# @return perl::ListReturn An array of Set<Int> The spanning trees, each encoded as"
		    "# the set of edges, where the edge indices are determined as follows: For i,j ="
		    "# 0,..,n-1 and i < j we sort the edges (i,j) in ascending lexicographical ordering"
		    "# , i.e, (0,1) < (0,2) < .. < (n-2,n-1).",
		    &spanning_complete,"spanning_complete($)");
  
}}
