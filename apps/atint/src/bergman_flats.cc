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
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/Map.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/basicoperations.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/diagonal_functions.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  
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
  perl::ListReturn returnFlats(perl::Object matroid) {
    Vector<Vector<Set<int> > > flats;
    Vector<Vector<int> > chains;
    
    chains_of_flats(matroid,flats,chains);
    
    perl::ListReturn result;
    for(int i = 0; i < flats.dim(); i++) {
      result << IncidenceMatrix<>(flats[i]);
    }
    for(int j = 0; j < chains.dim(); j++) {
      result << chains[j];
    }
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief This function takes a list of flats and chains of flats of a matroid and creates the corresponding bergman fan
   @param int n The number of elements of the matroid's ground set
   @param Vector<Vector<Set<int>>> flats A list of flats, as created by chains_of_flats
   @param Vector<Vector<int>> chains A list of chains of flats, as created by chains_of_flats
   @param bool mod_out_lineality Optional. Whether the lineality space should be modded out by setting the last coordinate to be the sum of the remaining ones. True by default.
   @return perl::Object A WeightedComplex object representing the bergman fan of the given flat lattice. Note that the list of rays will be in the same order as the list of corresponding flats
   */
  perl::Object bergman_fan_from_chains(int n, Vector<Vector<Set<int> > > flats, 
							Vector<Vector<int> > chains, bool mod_out_lineality = true) {
    Matrix<Rational> rays(0,mod_out_lineality? n-1 : n);
    //Go through all flats and create appropriate rays
    for(int r = 0; r < flats.dim(); r++) {
      for(int f = 0; f < flats[r].dim(); f++) {
	//The last element, n-1 is modded out as sum of the rest
	Set<int> flat = (flats[r])[f];
	Vector<Rational> matray(mod_out_lineality? n-1 : n);
	if( flat.contains(n-1) && mod_out_lineality) {
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
    
    Matrix<Rational> lineality(0,n);
      lineality /= ones_vector<Rational>(n);
    
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << cones;
      result.take("TROPICAL_WEIGHTS") << ones_vector<Integer>(cones.dim());
      result.take("IS_UNIMODULAR") << true;
      if(!mod_out_lineality) result.take("LINEALITY_SPACE") << lineality;
      
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object bergman_by_flats(perl::Object matroid, bool mod_out_lineality = true) {
    int rank = matroid.give("RANK");
    int n = matroid.give("N_ELEMENTS");
    
    //Compute flats and chains of flats
    Vector<Vector<Set<int> > > flats(rank-1); 
    Vector<Vector<int> > chains;
    chains_of_flats(matroid, flats,chains);
    
    return bergman_fan_from_chains(n, flats,chains, mod_out_lineality);
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object matroid_intersection_by_flats(perl::Object X, perl::Object Y, perl::Object matroid) {
    //Extract values
    int Xdim = X.give("CMPLX_DIM");
    int Xambi = X.give("CMPLX_AMBIENT_DIM");
    int Ydim = Y.give("CMPLX_DIM");
    int n = matroid.give("N_ELEMENTS");
    Array<Set<int> > bases = matroid.give("BASES");
    int rank = matroid.give("RANK");
    
    //Whether we need to lift X and Y
    bool need_to_lift = Xambi < n;
    if(need_to_lift) {
      X = lift_variety(X);
      Y = lift_variety(Y);
    }
        
    //dbgtrace << "Checking codimension" << endl;
    
    //If the codimensions of the varieties add up to something larger then CMPLX_AMBIENT_DIM, return the 0-cycle 
    if(Xdim + Ydim - (rank-1) < 0) {
      return CallPolymakeFunction("zero_cycle");
    }
    
    //dbgtrace << "Computing flats and chains " << endl;
    
    //Compute the flats of the matroid fan
    Vector<Vector<Set<int> > > flats;
    Vector<Vector<int> > chains;
    chains_of_flats(matroid,flats,chains);
    
    //dbgtrace << "Computing pairs of flats " << endl;
    
    //We add the empty set and the complete set for pairing
    Vector<Set<int> > empty_flat; empty_flat |= Set<int>();
    Vector<Vector<Set<int> > > empty_flat_vector; empty_flat_vector |= empty_flat;
    flats = empty_flat_vector | flats;
    Vector<Set<int> > all_flat; all_flat |= sequence(0,n);
    flats |= all_flat;
    
    //We now compute all pairs of flats for the matroid sum (matroid + matroid)
    //and the values of the diagonal functions phi_i
    Vector<Vector<Set<int> > > sum_flats(2*rank - 1);
    Vector<Matrix<Rational> > phi_by_flats(2*rank -1);
      for(int r = 0; r < phi_by_flats.dim(); r++) phi_by_flats[r] = Matrix<Rational>(0,rank);
    
    Map<int,int> shift_map;
    for(int i = 0; i < n; i++) { shift_map[i] = i+n;}
    for(int firstrank = 0; firstrank < flats.dim(); firstrank++) {
      for(int firstflat = 0; firstflat < flats[firstrank].dim(); firstflat++) {
	Set<int> flat1 = (flats[firstrank])[firstflat];
	
	for(int secondrank = 0; secondrank < flats.dim(); secondrank++) {
	  for(int secondflat = 0; secondflat < flats[secondrank].dim(); secondflat++) {
	      Set<int> flat2 = (flats[secondrank])[secondflat];
	      
	      //We ignore the pairing of two empty or two complete sets
	      if(flat1.size() == flat2.size() && (flat1.size() == 0 || flat1.size() == n)) continue;
	      
	      //Compute the rank of the union
	      Set<int> flatsum = flat1 + flat2;
	      int srank = 0;
	      for(int b = 0; b < bases.size(); b++) {
		int ib = (bases[b] * flatsum).size();
		if(ib > srank) srank = ib;
	      }//END compute rank
	      
	      //Insert the sum of flats
	      Set<int> shifted_flat = attach_operation(
		flat2, pm::operations::associative_access<Map<int,int>, int> (&shift_map) );
	      
	      sum_flats[firstrank + secondrank -1] |= (flat1 + shifted_flat);
	      
	      //Compute its function vector, i.e. the value of each phi_i on this flat pair
	      int phi_sum = firstrank + secondrank - srank;
	      Vector<Rational> phi_flat(rank);
	      if(phi_sum > 0) {
		phi_flat.slice(sequence(0, phi_sum)) = (-ones_vector<Rational>(phi_sum));
	      }
	      phi_by_flats[firstrank + secondrank -1] /= phi_flat;	      
	      
	      //dbgtrace << flat1 << "," << flat2 << ", " << phi_flat << endl;
	      
	  }//END iterate second flat
	}//END iterate rank of second flat
	
	
      }//END iterate first flat
    }//END iterate rank of first flat
    
    //dbgtrace << "Computing chains of flat pairs " << endl;
    
    //Now compute chains of flat pairs
    Vector<Vector<int> > sum_chains;
    //Start with rank-1-flats
    for(int i = 0; i < sum_flats[0].dim(); i++) {
      Vector<int> oneflat; oneflat |= i;
      sum_chains |= oneflat;
    }
    
    //Now subsequently go through rank k - flats and see if you can attach them to a lower chain
    for(int k = 1; k < sum_flats.dim(); k++) {
      Vector<Vector<int> > newchains;
      for(int f = 0; f < sum_flats[k].dim(); f++) {
	//Go through all chains
	for(int c = 0; c < sum_chains.dim(); c++) {
	  Set<int> lastflat = (sum_flats[k-1])[ (sum_chains[c])[sum_chains[c].dim()-1] ];
	  if( (lastflat * (sum_flats[k])[f]).size() == lastflat.size()) {
	    Vector<int> attached_chain =  sum_chains[c];
	    attached_chain |= f;
	    newchains |= attached_chain;
	  }
	}//END iterate rank k-1 chains
      }//END iterate rank k flats
      sum_chains = newchains;
    }//END iterate ranks
   
    //Now we create the matroid fan B(M + M)
    perl::Object sum_fan = bergman_fan_from_chains(2*n, sum_flats, sum_chains,false);
    
    //dbgtrace << "Computing divisors of diagonal functions " << endl;
    
    //Compute the product of X and Y
    std::vector<perl::Object> XandY;
	XandY.push_back(X); XandY.push_back(Y);
    perl::Object Z = compute_product_complex_lattice(XandY);
    
    //Concatentate function values
    Matrix<Rational> phi_matrix(0,rank);
    for(int r = 0; r < phi_by_flats.dim(); r++) {
      phi_matrix /= phi_by_flats[r];
    }
    
    //Now subsequently compute the divisor of phi_i
    for(int i = 1; i <= rank; i++) {
      //dbgtrace << "Applying function " << i << endl;
      perl::Object phi_i("RationalFunction");
	phi_i.take("DOMAIN") << sum_fan;
	phi_i.take("RAY_VALUES") << phi_matrix.col(i-1);
	phi_i.take("LIN_VALUES") << ones_vector<Rational>(1);
	
      Z = divisor_rational(Z,phi_i);
    }
    
    //Finally we project
    Matrix<Rational> zrays = Z.give("RAYS");
    Matrix<Rational> zlineality = Z.give("LINEALITY_SPACE");
    IncidenceMatrix<> zcones = Z.give("MAXIMAL_CONES");
    Vector<Integer> zweights = Z.give("TROPICAL_WEIGHTS");
    bool uses_homog = Z.give("USES_HOMOGENEOUS_C");
    
    
    //dbgtrace << "Rays: " << zrays << endl;
    //dbgtrace << "Cones: " << zcones << endl;
    //dbgtrace << "Lineality: " << zlineality << endl;
    //dbgtrace << "Weights: " << zweights << endl;
    
    zrays = zrays.minor(All, sequence(0, n - (uses_homog? 0 : 1)));
    zlineality = zlineality.minor(All, sequence(0, n - (uses_homog? 0 : 1)));
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << zrays;
      result.take("LINEALITY_SPACE") << zlineality;
      result.take("MAXIMAL_CONES") << zcones;
      result.take("TROPICAL_WEIGHTS") << zweights;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
    
    //If we lifted, project back
    return need_to_lift? project_variety(result) : result;
  }//END matroid_intersection_by_flats

  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Matroid fans"
		    "# This function computes the bergman fan of a matroid in the subdivision associated"
		    "# to its lattice of flats. In particular it computes all flats and all maximal"
		    "# chains of flats. Hence this function is terribly inefficient and should only be used"
		    "# with comparatively small matroids"
		    "# @param matroid::Matroid m Any matroid"
		    "# @param Bool mod_out_lineality Optional. Whether the lineality space should be "
		    "# modded out by setting the last coordinate to be the sum of the remaining ones. True "
		    "# by default."
		    "# @return WeightedComplex The bergman fan B(m) in its flat subdivision",
		    &bergman_by_flats,"bergman_fan_flats(matroid::Matroid; $=1)");
  
  UserFunction4perl("# @category Intersection products"
		    "# Computes an intersection product of two cycles in a bergman fan"
		    "# @param WeightedComplex X A tropical cycle living either in B(M) or B(M)/L, where"
		    "# L = <(1,...,1)>"
		    "# @param WeightedComplex Y A tropical cycle, also living in B(M) or B(M)/L (both "
		    "# cycles should of course live in the same fan)"
		    "# @param matroid::Matroid M A matroid. X and Y live in B(M) or B(M)/L. The function "
		    "# will detect automatically which of these cases applies"
		    "# @return WeightedComplex The intersection product of X and Y in B(M) (possibly mod L)", 
		    &matroid_intersection_by_flats, "matroid_intersection_by_flats(WeightedComplex, WeightedComplex, matroid::Matroid)");
  
  UserFunction4perl("# @category Matroids"
	            "# This function computes all flats of a matroid and all maximal chains of flats. "
		    "# The function does so via a brute-force method and it is not recommended to try this "
		    "# on larger matroids"
		    "# @param matroid::Matroid M"
		    "# @return IncidenceMatrix An array containing first rank(M)-1 IncidenceMatrix objects, "
		    "# I_1,...,I_k, where I_j contains the rank-j-flats. Then the array contains Vector<Int>,"
		    "# c_0,...,c_l, where c_i[j] contains the row index of the j-th flat in I_{j+1}",
		    &returnFlats, "compute_matroid_flats(matroid::Matroid)");
  
  
  
}}