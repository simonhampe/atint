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
#include "polymake/Integer.h"

namespace polymake { namespace tropical {

  using namespace atint::donotlog;
  //using namespace atint::dolog;
  //using namespace atint::dotrace;
 
  //Documentation see perl wrapper
  Integer count_mn_cones(int n) {
    if(n == 3) {
      return Integer(1);
    }
    Integer result(1);
    Integer nint(n);
    for(Integer i(0); i <= n-4; i++) {
	result = result * (2*(nint-i) -5);
    }
    return result;
  }
  
  //Documentaion see perl wrapper
  Integer count_mn_rays(int n) {
    if(n == 3) {
      return Integer(0);
    }
    Integer result(0);
    Integer nint(n);
    for(long i = 1; i <= n-3; i++) {
      result = result + Integer::binom(nint-1,i);
    }
    return result;
  }
  
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
    int raycount = count_mn_rays(n);
    dbgtrace << "Expecting " << raycount << " rays" << endl;
    Matrix<Rational> rays(raycount,raydim);
    
    //Will contain the rays, but as set partitions (better for checking for doubles)
    //Ray i is the i-th row of rays
    //Vector<Set<int> > raysAsPartitions;
    
    //Will contain value 'true' for each ray that has been computed
    Vector<bool> raysComputed(count_mn_rays(n));
    //Will contain the set of maximal cones 
    Vector<Set<int> > cones;
      
    //Compute the number of sequences = number of maximal cones
    int noOfMax = count_mn_cones(n);   
    
    //Things we will need:
    Set<int> allLeafs = sequence(0,n); //The complete sequence of leaves (for taking complements)
    Vector<int> onlyones = ones_vector<int>(raydim); //A ones vector(for projecting the lineality space)
    Vector<Integer> rayIndices(n-2); //Entry k contains the sum from i = 1 to k of binomial(n-1,i)
      rayIndices[0] = 0;
      for(int i = 1; i < rayIndices.dim(); i++) {
	rayIndices[i] = rayIndices[i-1] + Integer::binom(n-1,i);
      }
    
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

	//Now we compute the index of the ray -----------------------------------------
	// Consider ray as a vector of length n filled with a's and b's where entry i is a iff i is in I
	int k = n - rayset.size();
	int bsleft = k-1; int l = 1;
	int rIndex (rayIndices[k-2]);
	while(bsleft > 1) {
	  if(rayset.contains(n-l-1)) {
	    rIndex += (int)Integer::binom(n-l-1,bsleft-1);
	  }
	  else {
	    bsleft--;
	  }
	  l++;
	}
	int m = 0;
	while(rayset.contains(m)) { m++;}
	//at last we add the difference of the indices of the second b' and the first b (-1)
	rIndex += (n-1-l)-m; 
	newcone = newcone + rIndex;
	
	
	
	//If not, create the corresponding matroid coordinates
	if(!raysComputed[rIndex]) {
	//if(!found) {
	  dbglog << "Ray index of " << rayset << " is " << rIndex << endl;
	  dbgtrace << "Ray " << rayset << " does not exist. Creating..." << endl;
	  //raysAsPartitions = raysAsPartitions | rayset;
	  //newcone = newcone + (raysAsPartitions.dim()-1);
	  raysComputed[rIndex] = true;
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
	  //rays = rays / newray;
	  rays.row(rIndex) = newray;
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
		    "# Computes the number of maximal cones of the tropical moduli space M_0,n"
		    "# @param Int n The number of leaves. Should be >= 3"
		    "# @return Integer The number of maximal cones",
		    &count_mn_cones,"count_mn_cones($)");
		    
  UserFunction4perl("# @category Tropical geometry"
		    "# Computes the number of rays of the tropical moduli space M_0,n"
		    "# @param int n The number of leaves. Should be >= 3"
		    "# @return Integer The number of rays",
		    &count_mn_rays,"count_mn_rays($)");

  UserFunction4perl("# @category Tropical geometry"
		    "# Creates the moduli space of abstract rational n-marked curves. Its coordinates are"
		    "# given as the coordinates of the bergman fan of the matroid of the complete graph on "
		    "# n-1 nodes (but not computed as such)"
		    "# @param Int n The number of leaves. Should be at least 4"
		    "# @return fan::PolyhedralFan The tropical moduli space M_0,n",
		    &tropical_mn, "tropical_m0n($)");
  
}}
