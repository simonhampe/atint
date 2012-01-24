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
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/moduli.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes all paired ordered Prüfer sequences of order n, i.e. all sequences of length 2n-4 where each element from n,...,2n-3 occurs twice and where removing every second occurence of an element yields an ascending sequence)
   @param int n Should be greater equal 3
   @return Vector<Vector<int> > A list of index sequences
   */
  Vector<Vector<int> > computePrueferSequences(int n) {
    if(n < 3) throw std::runtime_error("Cannot compute M_n cones for n < 3");
    
    //Compute number of sequences
    int noOfCones = count_mn_cones(n);
    
    Vector<Vector<int> > result;
    
    //The idea of the algorithm is the following: Indices describes the Prüfer sequence in the 
    //following way:
    //If the entry at position i is j, then the second occurence of i is placed at the (j+1)-st
    // *possible* (i.e. not yet filled) position. Note that the first occurence of each i is alredy
    // determined by the placement of all entries lower than i.
    Vector<int> indices = ones_vector<int>(n-2);
    for(int iteration = 0; iteration < noOfCones; iteration++) {
      Vector<int> baseSequence = zero_vector<int>(2*n -4);
      for(int i = 0; i < n-1; i++) {
	//Go through the zero entries of baseSequence. If it is the first or the indices[i]+1-th, 
	//insert an n+i
	int zero_count = -1;
	for(int entry = 0; entry < baseSequence.dim(); entry++) {
	  if(baseSequence[entry] == 0) {
	      zero_count++;
	      if(zero_count == 0) {
		baseSequence[entry] = n+i;
	      }
	      if(zero_count == indices[i]) {
		baseSequence[entry] = n+i;
		break;
	      }
	  }
	}
      }
      result |= baseSequence;
      //Increase the indices vector by "1"---------------------------------------------------    
      if(iteration < noOfCones-1) {
	int counterindex = n-3;
	while(indices[counterindex] == 2*(n-counterindex)-5) {
	  indices[counterindex] = 1;
	  counterindex--;
	}
	indices[counterindex]++;
      }
    }
    
    return result;
  } //END computePrueferSequences
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Takes a ordered paired Prüfer sequence (as for example computed by computePrueferSequences and decodes it into a list of edge partitions of the combinatorial type
   @param Vector<int> pseq An ordered paired Prüfer sequence
   @return Vector<Set<int> > A list of partitions of the corresponding combinatorial type
   */
  Vector<Set<int> > decodePrueferSequence(const Vector<int> &pseq) {
    //Construct vertex set
    int n = pseq[0]; //The first element is always the number of leaves
    Set<int> V = sequence(0,2*n-2);
    Vector<Set<int> > adjacencies(n-2); //Which leaves lie "behind" which interior vertex?
    Vector<Set<int> > result;
    
    dbgtrace << "Connecting leaves" << endl;
    int firstindex = 0; //We pretend that pseq starts at this index
    //Connect leaves
    for(int i = 0; i < n; i++) {
      adjacencies[pseq[firstindex]-n] += i;
      V = V - i;
      firstindex++;
    }//END add leaves
    
    dbgtrace << "Connecting edges" << endl;
    dbgtrace << "V: " << V << endl;
    dbgtrace << "Adjacencies: " << adjacencies << endl;
    
    //Now create edges
    for(int i = 1; i <= n-3; i++) {
      //If there are only two vertices left, connect them
      if(i == n-3) {
	Vector<int> lasttwo(V);
	result |= adjacencies[lasttwo[0]-n];
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
	result |= Av;
	adjacencies[pseq[firstindex]-n] += Av;
	V = V - smallest;
	firstindex++;
      }
    }//END create edges
    
    return result;
  }//END decodePrueferSequence
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object local_mn(std::vector<perl::Object> curves) {
    
    //Compute the set of all vertex valences > 3
    Set<int> valences;
    for(unsigned int c = 0; c < curves.size(); c++) {
	Vector<int> degrees = curves[c].give("NODE_DEGREES");
	valences += (Set<int>(degrees) - 3);
    }
    dbgtrace << "Valences: " << valences << endl;
    
    dbgtrace << "Computing Prüfer sequences" << endl;
    //Now we compute combinatorially all M_0,v for all possible valences v
    Map<int, Vector<Vector<int> > > combinatorial_mns;
    for(Entire<Set<int> >::iterator v = entire(valences); !v.at_end();v++) {
      combinatorial_mns[*v] = computePrueferSequences(*v);
    }
    dbgtrace << "Done." << endl;
    
    //Now iterate through all curves and compute their adjacent maximal cones
    
    
    
    
    perl::Object result("WeightedComplex");
      return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry / Moduli spaces / Local geometry"
		    "# Computes the moduli space M_0,n locally around a given list of combinatorial"
		    "# types. More precisely: It computes the weighted complex consisting of all"
		    "# maximal cones containing any of the given combinatorial types and localizes "
		    "# at these types "
		    "# @param Int n Should be >= 3"
		    "# @return WeightedComplex The local complex",
		    &local_mn,"local_m0n(;@)");
  
  Function4perl(&decodePrueferSequence,"dcp(Vector<Int>)");
}}