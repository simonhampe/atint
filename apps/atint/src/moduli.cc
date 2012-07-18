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
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/Integer.h"
#include "polymake/PowerSet.h"
#include "polymake/Array.h"
#include "polymake/atint/moduli.h"

namespace polymake { namespace atint {

  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
 
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Helper function to compute (small) binomial coefficients without the help of Integer to avoid casts
   */
  int binomial_int(int n, int k) {
    int num = n;
    for(int j = 1; j <= k-1; j++) { num *= n - j;}
    int den = k;
    for(int j = k-1; j > 1; j--) { den *= j;}
    return num/den;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  int count_mn_cones_int(int n) {
    if(n == 3) {
      return 1;
    }
    int result = 1;
    int nint = n;
    for(int i = 0; i <= n-4; i++) {
	result = result * (2*(nint-i) -5);
    }
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Does exactly the same as count_mn_rays, but returns an int. Only works for n<=12, since larger values produce too large integers
   */
  int count_mn_rays_int(int n) {
    if(n == 3) {
      return 0;
    }
    int result = 0;
    int nint = n;
    for(long i = 1; i <= n-3; i++) {
      result = result + binomial_int(nint-1,i);
    }
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  Matrix<int> pair_index_map(int n) {
    Matrix<int> E(n,n);
    int nextindex = 0;
    for(int i = 0; i < n-1; i++) {
      for(int j = i+1; j < n; j++) {
	E(i,j) = E(j,i) = nextindex;
	nextindex++;
      }
    }
    return E;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  Vector<Set<int> > decodePrueferSequence(const Vector<int> &pseq, int n) {
    //Construct vertex set
    if(n < 0) { n = pseq[0];} //The first element is always the number of leaves
    //Compute number of bounded edges
    int no_of_edges = pseq.dim() - n +1;
    Set<int> V = sequence(0,n + no_of_edges +1);
    Vector<Set<int> > adjacencies(no_of_edges +1); //Which leaves lie "behind" which interior vertex?
    Vector<Set<int> > result;
    Set<int> allLeafs = sequence(0,n);
    
    //dbgtrace << "Connecting leaves" << endl;
    int firstindex = 0; //We pretend that pseq starts at this index
    //Connect leaves
    for(int i = 0; i < n; i++) {
      adjacencies[pseq[firstindex]-n] += i;
      V = V - i;
      firstindex++;
    }//END add leaves
    
    //dbgtrace << "Connecting edges" << endl;
    //dbgtrace << "V: " << V << endl;
    //dbgtrace << "Adjacencies: " << adjacencies << endl;
    
    //Now create edges
    for(int i = 1; i <= no_of_edges; i++) {
      Set<int> rayset;
      //If there are only two vertices left, connect them
      if(i == no_of_edges) {
	Vector<int> lasttwo(V);
	rayset = adjacencies[lasttwo[0]-n];
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
	rayset = Av;
	adjacencies[pseq[firstindex]-n] += Av;
	V = V - smallest;
	firstindex++;
	
	
      }
      
      //If rayset contains the last leaf, take the complement
      if(rayset.contains(n-1)) {
	rayset = allLeafs - rayset;
      }
      
      result |= rayset;
    }//END create edges
    
    return result;
  }//END decodePrueferSequence
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
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
	//dbgtrace << "Setting E(" << i << "," << j << ") = " << nextindex << endl;
	E(i,j) = nextindex;
	E(j,i) = nextindex;
	nextindex++;
      }
    }
    
    // We compute the set of all ordered paired Prüfer sequences on n,...,2n-3
    // (i.e. all sequences of length 2n-4 where each element from n,...,2n-3 occurs twice
    // and where removing every second occurence of an element yields an ascending sequence)
    //From each such Prüfer sequence we then construct a maximal cone
    
    //Will contain the rays of the moduli space in matroid coordinates
    int raydim = (n*(n-1))/2 - n;
    int raycount = count_mn_rays_int(n);
    //dbgtrace << "Expecting " << raycount << " rays" << endl;
    Matrix<Rational> rays(raycount,raydim);
    
    //Will contain the rays, but as set partitions (better for checking for doubles)
    //Ray i is the i-th row of rays
    //Vector<Set<int> > raysAsPartitions;
    
    //Will contain value 'true' for each ray that has been computed
    Vector<bool> raysComputed(count_mn_rays_int(n));
    //Will contain the set of maximal cones 
    Vector<Set<int> > cones;
      
    //Compute the number of sequences = number of maximal cones
    int noOfMax = count_mn_cones_int(n);
    
    //Things we will need:
    Set<int> allLeafs = sequence(0,n); //The complete sequence of leaves (for taking complements)
    Vector<int> onlyones = ones_vector<int>(raydim); //A ones vector(for projecting the lineality space)
    Vector<int> rayIndices(n-2); //Entry k contains the sum from i = 1 to k of binomial(n-1,i)
      rayIndices[0] = 0;
      for(int i = 1; i < rayIndices.dim(); i++) {
	rayIndices[i] = rayIndices[i-1] + binomial_int(n-1,i);
      }
    
    //Iterate through all Prüfer sequences -------------------------------------------------
    
    Vector<int> indices = ones_vector<int>(n-2);
    Vector<int> baseSequence(2*n-4);
    Vector<Set<int> > adjacent(n-2); //These will be the partitions of the edges
    Vector<Rational> newray(raydim); //Container for new rays
    for(int iteration = 0; iteration < noOfMax; iteration++) {
      
      //Create the sequence currently represented by indices and append it------------------
      //dbgtrace << "Creating sequence" << endl;
//       baseSequence = zero_vector<int>(2*n -4);
      baseSequence.resize(2*n-4);
      baseSequence.fill(0);
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
      //dbgtrace << "Creating cone for sequence " << baseSequence << endl;
      Set<int> newcone;
      
      Set<int> V = sequence(0,2*n-2);
      //dbgtrace << "Initialized sequence to " << V << endl;
      adjacent.fill(Set<int>()); 
      //dbgtrace << "Connecting leaves" << endl;
      //First: Connect the leaves
      for(int i = 0; i < n; i++) {
	//dbgtrace << "Attaching leaf " << i << " to node " << baseSequence[0] << endl;
	adjacent[baseSequence[0]-n] = adjacent[baseSequence[0]-n] + i;
	V = V - i;
	baseSequence = baseSequence.slice(1,baseSequence.dim()-1);
      }
      //Now create edges:
      int enumber = n-3;
      //dbgtrace << "Creating edges" << endl;
      for(int i = 1; i <= enumber; i++) {
	//dbgtrace << "Creating edge number " << i << endl;
	//Construct the leaf partition represented by the curve corresponding to the sequence
	Set<int> rayset;
	if(i == enumber) { //If V only has two elements left, simply connect these
// 	  Vector<int> last(V);
	  //dbgtrace << "Only two left: " <<  V << endl;
	  rayset = adjacent[*(V.begin()) - n];//adjacent[last[0]-n];
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
	//dbgtrace << "Edge partition is " << rayset << ". Creating matroid coords" << endl;
	if(rayset.contains(n-1)) {
	  rayset = allLeafs - rayset;
	}
	//Now check, if we already have that ray
	//dbgtrace << "Checking if ray already exists" << endl;
	
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
	    rIndex += binomial_int(n-l-1,bsleft-1);
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
	  //dbgtrace << "Ray index of " << rayset << " is " << rIndex << endl;
	  //dbgtrace << "Ray " << rayset << " does not exist. Creating..." << endl;
	  //raysAsPartitions = raysAsPartitions | rayset;
	  //newcone = newcone + (raysAsPartitions.dim()-1);
	  raysComputed[rIndex] = true;
// 	  Vector<int> raylist(rayset);
	  newray.fill(0);
// 	  for(int k = 0; k < raylist.dim()-1; k++) {
// 	      for(int l = k+1; l < raylist.dim(); l++) {
	  for(pm::Subsets_of_k_iterator<const pm::Set<int>& > raypair = entire(all_subsets_of_k(rayset,2)); !raypair.at_end(); raypair++) {
	      int newrayindex = E((*raypair).front(),(*raypair).back());
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
// 	      }
// 	  }
	  //rays = rays / newray;
	  rays.row(rIndex) = newray;
	}
      }//END iterate edges
      cones |= newcone;
    
      
      //dbgtrace << "Increasing counter" << endl;   
      //Increase the indices vector by "1"---------------------------------------------------    
      if(iteration < noOfMax-1) {
	int counterindex = n-3;
	while(indices[counterindex] == 2*(n-counterindex)-5) {
	  indices[counterindex] = 1;
	  counterindex--;
	}
	indices[counterindex]++;
      }
    }//END iterate cones
    
    std::ostringstream dsc;
      dsc << "Moduli space M_0," << n;
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << ones_vector<int>(cones.dim());
      result.set_description() << dsc.str();
      result.take("IS_UNIMODULAR") << true;
    return result;
    
  }

  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
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
		    "# The lineality space (1,..,1) is modded out by setting the last coordinate"
		    "# equal to minus the sum of the others"
		    "# The isomorphism to the space of curve metrics is obtained by choosing"
		    "# the last leaf as special leaf"
		    "# @param Int n The number of leaves. Should be at least 4"
		    "# @return WeightedComplex The tropical moduli space M_0,n",
		    &tropical_mn, "tropical_m0n($)");
  
  Function4perl(&decodePrueferSequence,"dcp(Vector<Int>)");
//   UserFunction4perl("",&adjacentRays,"adjacentRays(RationalCurve)");
  
}}
