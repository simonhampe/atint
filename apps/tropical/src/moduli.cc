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
#include "polymake/Matrix.h"
#include "polymake/Set.h"
#include "polymake/Rational.h"

namespace polymake { namespace tropical {

  using namespace atint::donotlog;
  //using namespace atint::dolog;
  //using namespace atint::dotrace;
 
  
  //Documentation see perl wrapper
  perl::Object tropical_mn(int n) {
    if(n < 4) {
      throw std::runtime_error("Number of leaves should be at least 4 for M_0,n computation");
    }
    
    //First we create the edge index matrix E(i,j) that contains at element i,j the edge index of edge (i,j)
    //in the complete graph on n-1 nodes
    int nextindex = 0;
    Matrix<int> E(n-1,n-1);
    for(int i = 0; i < n-2; i++) {
      for(int j = i+1; j < n-1; j++) {
	dbgtrace << "Setting E(" << i << "," << j << ") = " << nextindex << endl;
	E(i,j) = nextindex;
	nextindex++;
      }
    }
    
    // We compute the set of all ordered paired Prüfer sequences on n,...,2n-3
    // (i.e. all sequences of length 2n-4 where each element from n,...,2n-3 occurs twice
    // and where removing every second occurence of an element yields an ascending sequence)
    //From each such Prüfer sequence we then construct a maximal cone
    
    //Will contain the rays of the moduli space in matroid coordinates
    int raydim = (n*(n-1))/2 - n;
    Matrix<Rational> rays(0,raydim);
    //Will contain the rays, but as set partitions (better for checking for doubles)
    //Ray i is the i-th row of rays
    Vector<Set<int> > raysAsPartitions;
    //Will contain the set of maximal cones 
    Vector<Set<int> > cones;
      
    //Compute the number of sequences = number of maximal cones
    int noOfMax = 1;
    for(int i = 0; i <= n-4; i++) {
      noOfMax *= (2*(n-i) -5);
    }
    
    //Things we will need:
    Set<int> allLeafs = sequence(0,n);
    Vector<int> onlyones = ones_vector<int>(raydim);
    
    //Iterate through all Prüfer sequences -------------------------------------------------
    
    Vector<int> indices = ones_vector<int>(n-2);
    for(int iteration = 0; iteration < noOfMax; iteration++) {
      
      //Create the sequence currently represented by indices and append it------------------
      dbgtrace << "Creating sequence" << endl;
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
      
      //We now decode the Prüfer sequence to obtain the corresponding cone---------------------
      dbgtrace << "Creating cone for sequence " << baseSequence << endl;
      Set<int> newcone;
      
      Set<int> V = sequence(0,2*n-2);
      dbgtrace << "Initialized sequence to " << V << endl;
      Vector<Set<int> > adjacent(n-2); //These will be the partitions of the edges
      dbgtrace << "Connecting leaves" << endl;
      //First: Connect the leaves
      for(int i = 0; i < n; i++) {
	dbgtrace << "Attaching leaf " << i << " to node " << baseSequence[0] << endl;
	adjacent[baseSequence[0]-n] = adjacent[baseSequence[0]-n] + i;
	V = V - i;
	baseSequence = baseSequence.slice(1,baseSequence.dim()-1);
      }
      //Now create edges:
      int enumber = n-3;
      dbgtrace << "Creating edges" << endl;
      for(int i = 1; i <= enumber; i++) {
	dbgtrace << "Creating edge number " << i << endl;
	//Construct the leaf partition represented by the curve corresponding to the sequence
	Set<int> rayset;
	if(i == enumber) { //If V only has two elements left, simply connect these
	  Vector<int> last(V);
	  dbgtrace << "Only two left: " << last << " created from " << V << endl;
	  rayset = adjacent[last[0]-n];
	}
	else {
	  Set<int> pset(baseSequence); 
	  int smallest = -1;
	  //Find the smallest element in V that is not in P 
	  for(Entire<Set<int> >::iterator vit = entire(V); !vit.at_end(); vit++) {
	    if(!(pset.contains(*vit))) {
	      smallest = *vit; break;
	    }
	  }
	  rayset = adjacent[smallest-n];
	  //Add the leaves of this partition to the adjacency of the newly connected p_i
	  adjacent[baseSequence[0]-n] = adjacent[baseSequence[0]-n] + adjacent[smallest-n];
	  //Remove v and p_i
	  V = V - smallest;
	  baseSequence = baseSequence.slice(1,baseSequence.dim()-1);
	}
	//The new edge is: v_{adjacent[smallest]}. If it containst the last leaf, take the complement
	dbgtrace << "Edge partition is " << rayset << ". Creating matroid coords" << endl;
	if(rayset.contains(n-1)) {
	  rayset = allLeafs - rayset;
	}
	//Now check, if we already have that ray
	dbgtrace << "Checking if ray already exists" << endl;
	
// 	bool found = false;
// 	for(int s = 0; s < raysAsPartitions.dim(); s++) {
// 	  if(raysAsPartitions[s] == rayset) {
// 	      newcone = newcone + s;
// 	      found = true;
// 	      break;
// 	  }
// 	}
	//If not, create the corresponding matroid coordinates
	if(!found) {
	  dbgtrace << "Ray " << rayset << " does not exist. Creating..." << endl;
	  raysAsPartitions = raysAsPartitions | rayset;
	  newcone = newcone + (raysAsPartitions.dim()-1);
	  Vector<int> raylist(rayset);
	  Vector<Rational> newray(raydim);
	  for(int k = 0; k < raylist.dim()-1; k++) {
	      for(int l = k+1; l < raylist.dim(); l++) {
		int newrayindex = E(raylist[k],raylist[l]);
		//If the newrayindex is one higher than the ray dimension, 
		//this means it is first of all the last pair. Also, we don't
		//add -e_n but e_1 + ... + e_{n-1} (as we mod out lineality)
		if(newrayindex < raydim) {
		    newray[newrayindex] = -1;
		}
		else {
		    newray = newray + onlyones;
		}
	      }
	  }
	  rays = rays / newray;
	}
      }
      cones = cones | newcone;
    
      
      dbgtrace << "Increasing counter" << endl;   
      //Increase the indices vector by "1"---------------------------------------------------    
      if(iteration < noOfMax-1) {
	int counterindex = n-3;
	while(indices[counterindex] == 2*(n-counterindex)-5) {
	  indices[counterindex] = 1;
	  counterindex--;
	}
	indices[counterindex]++;
      }
    }
    
    perl::Object result("fan::PolyhedralFan");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << ones_vector<int>(cones.dim());
    return result;
    
  }

  UserFunction4perl("# @category Tropical geometry"
		    "# Creates the moduli space of abstract rational n-marked curves. Its coordinates are"
		    "# given as the coordinates of the bergman fan of the matroid of the complete graph on "
		    "# n-1 nodes (but not computed as such)"
		    "# @param Int n The number of leaves. Should be at least 4"
		    "# @return fan::PolyhedralFan The tropical moduli space M_0,n",
		    &tropical_mn, "tropical_m0n($)");
  
}}
