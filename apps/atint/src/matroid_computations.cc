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
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 This file provides convenience methods for creating locally restricted
 tropical varieties from given varieties
 */


#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  perl::Object matroidFromFan(perl::Object fan) {
    
    //Dehomogenize, if necessary.
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    if(uses_homog) {
	fan = fan.CallPolymakeMethod("dehomogenize");
    }
    
    //Extract values
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    IncidenceMatrix<> cones = fan.give("MAXIMAL_CONES");
    int rank = fan.give("CMPLX_DIM");
    int no_of_el = fan.give("CMPLX_AMBIENT_DIM");
    
    //Compute all potential bases
    Set<Set<int> > bases(pm::all_subsets_of_k(sequence(0,no_of_el),rank));
    
    //Now go through all potential flat vectors of cardinality >= rank. The corr. set is a flat, iff
    // the fan contains that vector. If it is a flat, remove all bases contained in it.
    for(int card = rank; card < no_of_el; card++) {
	Array<Set<int> > potential_flats = pm::all_subsets_of_k(sequence(0,no_of_el),card);
	for(int flat = 0; flat < potential_flats.size(); flat++) {
	    Vector<Rational> inc_vec(no_of_el);
	      inc_vec.slice(potential_flats[flat]) = -ones_vector<Rational>(card);
	    bool is_in_fan = false;
	    for(int mc = 0; mc < cones.rows(); mc++) {
	      if(is_ray_in_cone(rays.minor(cones.row(mc),All),linspace,inc_vec)) {
		  is_in_fan = true; break;
	      }
	    }//END iterate maximal cones
	    if(is_in_fan) {
	      //Compute all subsets of cardinality the rank and remove them from the bases
	      bases -= Set<Set<int> >(pm::all_subsets_of_k(potential_flats[flat],rank));	      
	    }
	}//END iterate flats of given cardinality
	
    }//END iterate flat cardinalities
    
    Vector<Set<int> > final_bases(bases);
    
    perl::Object result("matroid::Matroid");
      result.take("N_ELEMENTS") << no_of_el;
      result.take("BASES") << final_bases;
      
    return result;
    
  }//END matroidFromFan
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object matroid_sum(perl::Object m1, perl::Object m2) {
    //Extract values
    int n1 = m1.give("N_ELEMENTS");
    int n2 = m2.give("N_ELEMENTS");
    int N = n1 + n2;
    
    Array<Set<int> > base1 = m1.give("BASES");
    Array<Set<int> > base2 = m2.give("BASES");
    
    Map<int,int> base_shift;
      for(int i = 0; i < n2; i++) { base_shift[i] = n1+i;}

    Vector<Set<int> > result_bases;
    for(int b1 = 0; b1 < base1.size(); b1++) {
      for(int b2 = 0; b2 < base2.size(); b2++) {
	result_bases |= Set<int>(base1[b1] + 
			Set<int>(attach_operation(base2[b2],
			pm::operations::associative_access<Map<int,int>,int>(&base_shift))));
      }
    }
    
    perl::Object result("matroid::Matroid");
      result.take("N_ELEMENTS") << N;
      result.take("BASES") << result_bases;
      
    return result;
  }//END matroid_sum
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::ListReturn find_matroid_amalgam(perl::Object m1, perl::Object m2, perl::Object n, Set<int> m1_set,Set<int> m2_set) {
    
    //Extract values
    int N1 = m1.give("N_ELEMENTS");
    int N2 = m2.give("N_ELEMENTS");
    int base_rank = n.give("RANK");
    
    //Compute flats of the base fan
    perl::ListResult n_flats = ListCallPolymakeFunction("compute_matroid_flats",n);
    
    pm::cout << "Computing product fan" << endl;
    
    //Compute the bergman fan of the sum of the two extensions
    perl::Object sum_matroid = matroid_sum(m1,m2);
    perl::Object sum_fan = CallPolymakeFunction("bergman_fan_flats",sum_matroid,0);
    
    Matrix<Rational> sum_rays = sum_fan.give("RAYS");
    
    //The product complex has lineality space <(1,...,1)>, which has always value 1
    Matrix<Rational> lin_values(base_rank,0);
      lin_values |= ones_vector<Rational>(base_rank);
    //We compute ray values by projecting onto N + N and computing the corresponding function value
    Matrix<Rational> ray_values(base_rank,sum_rays.rows());
    
    pm::cout << "Computing ray values" << endl;
    
    for(int r = 0; r < sum_rays.rows(); r++) {
      //Separate ray and project each part
      Vector<Rational> part1 = sum_rays.row(r).slice(sequence(0,N1)).slice(m1_set);
      Vector<Rational> part2 = sum_rays.row(r).slice(sequence(N1,N2)).slice(m2_set);
      //Compute the corresponding flats
      Set<int> flat1;
      Set<int> flat2;
      for(int c1 = 0; c1 < part1.dim(); c1++) {if(part1[c1] != 0) flat1 += c1;}
      for(int c2 = 0; c2 < part2.dim(); c2++) {if(part2[c2] != 0) flat2 += c2;}
      Set<int> flat_union = flat1 + flat2;
      //Compute rank(flat1),rank(flat2),rank(flat1 + flat2)
      int rank1 = 0;
      int rank2 = 0;
      //The union is not nec. a flat, so we need to find the smallest flat containing it
      int rankunion = base_rank; 
      for(int f = 0; f < base_rank-1; f++) {
	IncidenceMatrix<> fflats = n_flats[f];
	for(int compflat = 0; compflat < fflats.rows(); compflat++) {
	  if( (fflats.row(compflat) * flat1).size() == flat1.size()) rank1 = f+1;
	  if( (fflats.row(compflat) * flat2).size() == flat2.size()) rank2 = f+1;
	  if( (fflats.row(compflat) * flat_union).size() == flat_union.size()) {
	    if(f+1 < rankunion) rankunion = f+1;
	  }
	}//END iterate rank f+1 flats
      }//END iterate flat ranks
      if(rank1 == 0) rank1 = flat1.size() == 0? 0 : base_rank;
      if(rank2 == 0) rank2 = flat2.size() == 0? 0 : base_rank;
      if(flat_union.size() == 0) rankunion = 0;
      
      //Compute rank-delta and ray values
      int rank_delta = rank1 + rank2 - rankunion;
      //Insert ray values
      for(int i = 1; i <= ray_values.rows(); i++) {
	if(i <= rank_delta) ray_values(i-1,r) = -1;
      }
    }//END iterate sum fan rays
    
    //dbgtrace << "Ray values " << ray_values << endl;
    
    //Compute divisor of diagonal functions
    pm::cout << "Computing divisor" << endl;    
    Matrix<Rational> value_matrix = ray_values | lin_values;
    perl::Object fibre = divisorByValueMatrix(sum_fan, value_matrix);
    fibre = fibre.CallPolymakeMethod("dehomogenize");
    
    //Compute matroid
    pm::cout << "Computing matroid" << endl;
    perl::Object matroid = matroidFromFan(fibre);
    
    
    perl::ListReturn result;
      result << fibre;
      result << matroid;
    return result;
    
  }//END find_matroid_amalgam
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Matroids"
		    "# Takes a tropical fan (that should have degree 1) and computes the matroid that has this fan"
		    "# as bergman fan"
		    "# @param WeightedComplex fan A tropical fan that should be the Bergman fan of a matroid (NOT"
		    "# modulo lineality space). Be aware that this implementation uses (-1/0)-incidence vectors"
		    "# for flats, not (1/0)."
		    "# @return matroid::Matroid The corresponding matroid (given in terms of bases)",
		    &matroidFromFan,"matroid_from_fan(WeightedComplex)");
  
  UserFunction4perl("# @category Matroids"
		    "# Computes the sum of two matroids"
		    "# @param matroid::Matroid M1"
		    "# @param matroid::Matroid M2"
		    "# @return matroid::Matroid The sum, i.e. defined on the disjoint union of the ground"
		    "# sets and the bases are disjoint unions of bases",
		    &matroid_sum, "matroid_sum(matroid::Matroid, matroid::Matroid)");
  
  UserFunction4perl("# @category Matroids"
		    "# Takes two extensions M_1 and M_2 of a matroid N and tries to find a potential "
		    "# amalgam by computing the intersection-theoretic fibre product of the projections"
		    "# of the Bergman fans, then computing the corresponding matroid (assuming it is smooth)"
		    "# @param matroid::Matroid M1 An extension of N"
		    "# @param matroid::Matroid M2 An extension of N"
		    "# @param matroid::Matroid N A matroid"
		    "# @param Set<Int> M1_Set Tells which subset of the ground set of M_1 we have to "
		    "# restrict to obtain N"
		    "# @param Set<Int> M2_Set Tells which subset of the ground set of M_2 we have to "
		    "# restrict to obtain N"
		    "# @return An array, containing first the fibre product as WeightedComplex, then the"
		    "# corresponding matroid, as computed by [[matroid_from_fan]]",
		    &find_matroid_amalgam,"find_matroid_amalgam(matroid::Matroid, matroid::Matroid, matroid::Matroid, Set<Int>, Set<Int>)");
  
}}