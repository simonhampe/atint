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
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/cdd_helper_functions.h"

namespace polymake { namespace atint { 
  
    using namespace atintlog::donotlog;
//   using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;

//   //Checks circuit axioms
//   bool check_circuit_axioms(Vector<Set<int> > circuits) {
//     
//     //First we check pairwise non-containment and that no set is empty
//     for(int i = 0; i < circuits.dim(); i++) {
//       if(circuits[i].size() == 0) {
// 	  dbgtrace << "Empty set cannot be a circuit" << endl;
// 	  return false;
//       }
//       for(int j = i+1; j < circuits.dim(); j++) {
// 	Set<int> inter = circuits[i] * circuits[j];
// 	int isize = inter.size();
// 	if(isize == circuits[i].size() || isize == circuits[j].size()) {
// 	    dbgtrace << "Circuits cannot contain one another: " << circuits[i] << ", " << circuits[j] << endl;
// 	    return false;
// 	}
// 	//While we're at it, we go through the elements of the intersection and check circuit elimination
// 	for(Entire<Set<int> >::iterator x = entire(inter); !x.at_end(); x++) {
// 	  Set<int> must_contain = circuits[i] + circuits[j] - *x;
// 	  bool found_one = false;
// 	  for(int k = 0; k < circuits.dim(); k++) {
// 	    if(k != i && k != j) {
// 	      if(circuits[k].size() == (circuits[k] * must_contain).size()) {
// 		   found_one = true;
// 		   break;
// 	      }
// 	    }
// 	  }//END try to find a contained set
// 	  if(found_one == false) {
// 	      dbgtrace << "Circuit elimination fails for " << circuits[i] << ", " << circuits[j] << " and " << *x << endl;
// 	      return false;
// 	  }
// 	}//END check circuit elimination
// 	
//       }//END iterate j
//     }//END iterate i
//     
//     return true;
//   }//END check_circuit_axioms
  
//   /*
//    * In addition to a list of matroids that we will with all possibilities and the list of circuits we want, this takes
//    * the following arguments:
//    * - possible_sets: This list does not change through iterations. It contains all the subsets of 1...n that do not contain and are not contained in a set of list_of_circuits
//    * - associated_sets: Also does not change. Contains for each set of possible_sets the indices of 
//    * sets that are contained in it or contain it.
//    * - used_sets: The sets (as indices in possible_sets) which we currently add to our matroid.
//    * - allowed sets: The sets we are still allowed to add.
//    */
//   void recursiveCircuitListFiller(Vector<perl::Object> &listToFill, int n, const Vector<Set<int> > &list_of_circuits,
// 				  const Vector<Set<int> > &possible_sets, const Map<int,Set<int> > &associated_sets,
// 				  Set<int> used_sets, Set<int> allowed_sets) {
//     
//     //First of all, check if we already have a matroid
//     Vector<Set<int> > to_check(list_of_circuits);
//     for(Entire<Set<int> >::iterator addc = entire(used_sets); !addc.at_end(); addc++) {
// 	to_check |= possible_sets[*addc];
//     }
//     if(check_circuit_axioms(to_check)) {
// 	dbgtrace << "Checked " << to_check << " and succeeded" << endl;
// 	perl::Object m("matroid::Matroid");
// 	  m.take("N_ELEMENTS") << n;
// 	  m.take("CIRCUITS") << to_check;
// 	listToFill |= m;
// 	pm::cout << "Solutions found: " << listToFill.dim() << ", newest solution: " << to_check << endl;
//     }
//     
//     //Now we iterate through the allowed sets, add one element, then do a recursive call
//     Set<int> already_tested;
//     for(Entire<Set<int> >::iterator as = entire(allowed_sets); !as.at_end(); as++) {
//       
//       recursiveCircuitListFiller(listToFill,n, list_of_circuits, possible_sets, associated_sets,
// 				 used_sets + *as, allowed_sets - already_tested - associated_sets[*as] - *as);
//       
//       //We don't want this element allowed in future iterations
//       already_tested += *as;
//     }//END iterate allowed sets
//     
//     
//   } //END recursiveCircuitListFiller
  
//   /*
//    * Computes all matroids whose list of circuits contains a given list. Circuits should be given as 
//    * subsets of 0,..,n-1
//    */
//   perl::ListReturn computePossibleMatroids(int n, Vector<Set<int> > list_of_circuits) {
//     
//     //Compute a list of all subsets of n
//     Vector<Set<int> > subsets = Vector<Set<int> >(pm::all_subsets(sequence(0,n)));
//     //Now throw out all the ones that are contained in or contain an elment of list_of_circuits
//     for(int lc = 0; lc < list_of_circuits.dim(); lc++) {
//       Set<int> badsets;
//       for(int ss = 0; ss < subsets.dim(); ss++) {
// 	Set<int> inter = list_of_circuits[lc] * subsets[ss];
// 	if(inter.size() == list_of_circuits[lc].size() || inter.size() == subsets[ss].size() ) {
// 	    badsets += ss;
// 	}
//       }
//       subsets = subsets.slice(~badsets);
//     }//END iterate subsets
//     
//     //Now compute for each subset the list of (allowed) subsets that are contained in it or contain it
//     Map<int, Set<int> > associated_sets;
//     for(int s = 0; s < subsets.dim(); s++) {
// 	associated_sets[s] = Set<int>();
//     }
//     for(int i = 0; i < subsets.dim(); i++) {
//       for(int j = i+1; j < subsets.dim(); j++) {
// 	Set<int> inter = subsets[i] * subsets[j];
// 	if(inter.size() == subsets[i].size() || inter.size() == subsets[j].size()) {
// 	  associated_sets[i] += j;
// 	  associated_sets[j] += i;
// 	}
//       }//END iterate possible assoc subsets
//     }//END iterate subsets for finding associated sets
//     
//     //Now call the filler Function
//     Vector<perl::Object> listToFill;
//     
//     recursiveCircuitListFiller(listToFill, n, list_of_circuits, subsets, associated_sets, Set<int>(), sequence(0,subsets.dim()));
//     
//     perl::ListReturn result;
//       for(int r = 0; r < listToFill.dim(); r++) {
// 	  result << listToFill[r];
//       }
//     return result;
//     
//   }//END computePossibleMatroids
 
   /**
    * Takes a list of sets and checks the basis exchange axiom.
    * If force_loopfree = true, it checks that the corresponding matroid is also loopfree
    * The two references at the end will contain the minimal indices i and j and the minimal element x in B_i - B_j
    * such that the basis exchange axiom fails for these.
    * If the function fails for any other reason, it will set minimal_element to -1
    */
   bool check_basis_axioms(const Vector<Set<int> > &list_of_bases, int n, bool force_loopfree,
			  std::pair<int,int> &minimal_pair, int &minimal_element) {
     if(list_of_bases.dim() == 0) {
// 	dbgtrace << "There must be at least one bases" << endl;
	minimal_element = -1;
	return false;
     }
     if(force_loopfree) {
       Set<int> total = accumulate(list_of_bases,operations::add());
       if(total.size() != n) {
	 minimal_element = -1;
	  return false;
       }
     }//END check
     //Iterate over all pairs
     for(int i = 0; i < list_of_bases.dim(); i++) {
	for(int j = 0; j < list_of_bases.dim(); j++) {
	  if(i != j) {
	    //Take all elements in B_i - B_j and look for replacements in B_j - B_i
	    Set<int> iwithoutj = list_of_bases[i] - list_of_bases[j];
	    Set<int> jwithouti = list_of_bases[j] - list_of_bases[i];
	    for(Entire<Set<int> >::iterator x = entire(iwithoutj); !x.at_end(); x++) {
	      bool replacement_found = false;
	      for(Entire<Set<int> >::iterator y = entire(jwithouti); !y.at_end(); y++) {
		Set<int> newbasis = list_of_bases[i] - *x + *y;
  // 	      dbgtrace << "Trying " << newbasis << endl;
		for(int otherbases = 0; otherbases < list_of_bases.dim(); otherbases++) {
		  if(list_of_bases[otherbases] == newbasis) {
		    replacement_found = true; break;
		  }
		}//END compare to existing bases
		if(replacement_found) break;
	      }//END iterate possible replacements
	      
	      if(!replacement_found) {
//   	      dbgtrace << "No exchange element found for " << *x << " from " << list_of_bases[i] << " in " << list_of_bases[j] << endl;
		minimal_pair.first = i;
		minimal_pair.second = j;
		minimal_element = *x;
		return false;
	      }
		
	    }//END iterate iwithoutj
	  }//END if i != j
	}//END iterate j
     }//END iterate i
     
     
     return true;     
   }//END check_basis_axioms
   
   bool cbatest(int n, Vector<Set<int> > bases) {
     std::pair<int,int> p;
     int q;
      return check_basis_axioms(bases,n,false,p,q);
   }
   
   bool is_already_in_there(const Vector<Array<Set<int> > > &matroids, const Vector<Set<int> >&bases) {
      Set<Set<int> > bases_set(bases);
      dbgtrace << "Comparing " << matroids << " to " << bases << endl;
      for(int i = 0; i < matroids.dim(); i++) {
	  Set<Set<int> > mbaseset(matroids[i]);
	  if(mbaseset.size() == bases_set.size()) {
	    if( (bases_set * mbaseset).size() == bases_set.size()) {
		dbgtrace << "Is equal to number " << i << endl;
		return true;
	    }
	  }
      }
      return false;
   }
   
//   /*
//     * Takes a recursively filled list of matroids, a list of all possible bases,
//     * the list of bases we currently use and the list of bases we might still use.
//   */
//   void recursiveBasesListFiller(Vector<perl::Object> &listToFill, const Vector<Set<int> > &list_of_bases,
// 				Set<int> currently_used, Vector<int> allowed_sets, int n, int force_loopfree ) {
//     //First of all test, if the currently used list gives a valid matroid, if so: add it
// //     if( listToFill.dim() == 8) return;
//     Vector<Set<int> > current_matroid = list_of_bases.slice(currently_used);
//     
//     if(check_basis_axioms(current_matroid, n, force_loopfree)) {
// 	perl::Object matroid("matroid::Matroid");
// 	matroid.take("N_ELEMENTS") << n;
// 	matroid.take("BASES") << current_matroid;
// 	listToFill |= matroid;
// 	pm::cout << "Found matroid number " << listToFill.dim() <<" with bases " << current_matroid << endl;
// 	
//     }
//     
//     //Now we add the next basis on our list and recursively check
//     for(int as = 0; as < allowed_sets.dim(); as++) {
// 	recursiveBasesListFiller(listToFill, list_of_bases,
// 				 currently_used + allowed_sets[as],
// 				 allowed_sets.slice(~sequence(0,as+1)),n,force_loopfree);
//     }//END add next basis
//     
//   } //END recursiveBasesListFiller

   inline Array<Set<int> > arrayconv(Vector<Set<int> > v) {
      Array<Set<int> > result(v.dim());
      for(int i = 0; i < v.dim(); i++) {
	   result[i] = v[i];
      }
      return result;
   }

   void recursiveBasesListFiller(Vector<Array<Set<int> > > &listToFill, const Vector<Set<int> > &list_of_bases,
				 Set<int> currently_used,Set<int> allowed_sets, int n, bool force_loopfree) {
      //First of all test, if the currently used list gives a valid matroid, if so: add it
  //     if( listToFill.dim() == 8) return;
      Vector<Set<int> > current_matroid = list_of_bases.slice(currently_used);
      dbgtrace << "Iteration " << current_matroid << endl;
      std::pair<int,int> minimal_pair;
      int minimal_element;
      Set<int> elements_to_add; 
      if(check_basis_axioms(current_matroid, n, force_loopfree,minimal_pair, minimal_element)) {
	dbgtrace << "Is a matroid!" << endl;
	  if(!is_already_in_there(listToFill, current_matroid)) {
// 	    perl::Object matroid("matroid::Matroid");
// 	    matroid.take("N_ELEMENTS") << n;
	    Array<Set<int> > basisconv = arrayconv(current_matroid);
// 	    matroid.take("BASES") << basisconv;	    
	    listToFill = listToFill | basisconv;
// 	    pm::cout << "Found matroid number " << listToFill.dim() <<" with bases " << current_matroid << endl;
	  }
	  elements_to_add = allowed_sets;
	  
	  
      }
      else {
	//See if the test failed for another reason then the BEC
	if(minimal_element == -1) {
	  //If we have no basis yet, just add anything
	  if(currently_used.size() == 0) {
	      dbgtrace << "No bases yet" << endl;
	      elements_to_add = allowed_sets;
	  }
	  else {
	    //Otherwise find a loop and all sets containing it
	    Set<int> total_set = accumulate(current_matroid,operations::add());
	    int loop = *( (sequence(0,n) -  total_set).begin());
	    dbgtrace << "Have loop " << loop << endl;
	    for(Entire<Set<int> >::iterator as = entire(allowed_sets); !as.at_end(); as++) {
	      if(list_of_bases[*as].contains(loop)) elements_to_add += *as;
	    }
	  }
	}
	else {
	  dbgtrace << "Fixing basis exchange problem " << minimal_pair.first << ", " << minimal_pair.second << ", " << minimal_element << endl;
	  //If the basis exchange axiom is violated, find all allowed sets that fix this for the
	  // specific pair we got back.
  // 	dbgtrace << "Looking for fixing elements for " << current_matroid[minimal_pair.first] << ", " << 
  // 		current_matroid[minimal_pair.second] << ", " << minimal_element << endl;
	  for(Entire<Set<int> >::const_iterator fs = entire(allowed_sets); !fs.at_end(); fs++) {
	    if(!list_of_bases[*fs].contains(minimal_element)) {
	      Set<int> diff_with_bi = list_of_bases[*fs] - current_matroid[minimal_pair.first];
	      if(diff_with_bi.size() == 1) {
		  if( current_matroid[minimal_pair.second].contains( *(diff_with_bi.begin()))) {
  // 		    dbgtrace << list_of_bases[*fs] << " matches " << endl;
		      elements_to_add += *fs;
		  }
	      }
	    }
	  }//END find fixing sets
	}
      }
      
      dbgtrace << "Elements to consider for next iteration: " << list_of_bases.slice(elements_to_add) << endl;
      
      Set<int> not_to_add;
      for(Entire<Set<int> >::iterator ea = entire(elements_to_add); !ea.at_end(); ea++) {
	not_to_add += *ea;
	recursiveBasesListFiller(listToFill, list_of_bases,
				 currently_used + *ea, allowed_sets - not_to_add,n, force_loopfree);
      }//END add elements recursively
      
      
     
     
   }//END recursiveBasesListFiller
   
   /**
    * Computes a list of all possible matroids of rank r on n elements whose list of dependent sets
    * contains a specified list and whose bases contain a specified list.
    */
   perl::ListReturn computePossibleMatroids(int n, int r, bool force_loopfree, Vector<Set<int> > list_of_bases ,Vector<Set<int> > list_of_dependents) {
     
     //Compute a list of all possible bases - take all r-element subsets of (0..n) and remove
     // the ones that contain a set from the list. Also find the ones we definitely want to be bases.
     Vector<Set<int> > bases = Vector<Set<int> >(pm::all_subsets_of_k(sequence(0,n),r));
     Set<int> nonbases;
     for(int b = 0; b < bases.dim(); b++) {
       for(int d = 0; d < list_of_dependents.dim(); d++) {
	 if ((bases[b] * list_of_dependents[d]).size() == list_of_dependents[d].size()) {
	   nonbases += b; break;
	 }
       }//END iterate dependents
     }//END iterate r-subsets
     bases = bases.slice(~nonbases);
     Set<int> definite_bases;
     for(int b = 0; b < bases.dim(); b++) {
       for(int db = 0; db < list_of_bases.dim(); db++) {
	  if( (bases[b] * list_of_bases[db]).size() == list_of_bases[db].size()) {
	      definite_bases += b; break;
	  }
       }
     }
     
     
     
     //Now fill a list of matroids
     dbgtrace << "List of bases " << bases << endl;
     Vector<Array<Set<int> > > resultList;
     recursiveBasesListFiller(resultList, bases,
			      definite_bases,sequence(0,bases.dim()) - definite_bases, n,force_loopfree);
     
     perl::ListReturn result;
      for(int r = 0; r < resultList.dim(); r++) {
	  perl::Object m("matroid::Matroid");
	  m.take("N_ELEMENTS") << n;
	  m.take("BASES") << resultList[r];
	  result << m;
      }
     
     return result;
     
     
   }
   
   /**
    * Compute the set-theoretic intersection of the tropical circuits 
    */
   perl::ListReturn intersect_circuits(Matrix<Rational> circuits, int infty) {
     perl::Object linearspace = CallPolymakeFunction("atint::linear_nspace",circuits.cols());
     linearspace = linearspace.CallPolymakeMethod("homogenize");
     
     
     //Create the full space in homog. coordinates
     int n = circuits.cols();
     Matrix<Rational> rays = linearspace.give("RAYS");
     IncidenceMatrix<> cones = linearspace.give("MAXIMAL_CONES");
     Matrix<Rational> linspace = linearspace.give("LINEALITY_SPACE");
     
     perl::ListReturn ret;
     for(int c = 0; c < circuits.rows(); c++) {
//        pm::cout << "Intersecting circuit " << c+1 << " of " << circuits.rows() << endl;
       //Create equation
       Set<int> support;
       for(int i = 0; i < n; i++) {
	 if(circuits(c,i) != infty) support += i;
       }
       Matrix<Rational> fmatrix = unit_matrix<Rational>(n);
       perl::Object function("atint::MinMaxFunction");
	  function.take("LINEAR_COEFFICIENTS") << fmatrix.minor(support,All);
	  function.take("CONSTANT_COEFFICIENTS") << circuits.row(c).slice(support);
	  function.take("USES_MIN") << true;
       perl::Object hyperplane = CallPolymakeFunction("divisor", linearspace, function);
       ret << hyperplane;
       Matrix<Rational> hrays = hyperplane.give("RAYS");
       IncidenceMatrix<> hcones = hyperplane.give("MAXIMAL_CONES");
       Matrix<Rational> hlin = hyperplane.give("LINEALITY_SPACE");
       
//        fan_intersection_result inter = cdd_fan_intersection(	
// 	      rays, linspace, cones,
// 	      hrays, hlin, hcones,
// 	      true);
//        
//        rays = inter.rays;
//        cones = inter.cones;
//        linspace = inter.lineality_space;       
       
     }//END iterate circuits
   
   return ret;
   
   //Clear up redundant maximal cones
//    Set<int> redundant;
//    for(int i = 0; i < cones.rows(); i++) {
//      if(!redundant.contains(i)) {
// 	for(int j = i+1; j < cones.rows(); j++) {
// 	    Set<int> inter = cones.row(i) * cones.row(j);
// 	    if(inter.size() == cones.row(i).size()) redundant += i;
// 	    if(inter.size() == cones.row(j).size()) redundant += j;
// 	}
//      }     
//    }//END clear up cones
//    cones = cones.minor(~redundant,All);
//    
//    //Create result
//    perl::Object result("WeightedComplex");
//     result.take("RAYS") << rays;
//     result.take("MAXIMAL_CONES") << cones;
//     result.take("LINEALITY_SPACE") << linspace;
//     result.take("USES_HOMOGENEOUS_C") << true;
//    
//    return result;
   
   }//END intersect_circuits
   
  
  
  // --------------- PERL WRAPPERS ------------------------
  
//   Function4perl(&check_circuit_axioms,"cca(Vector<Set<Int> >)");
  Function4perl(&cbatest,"cba($,Vector<Set<Int> >)");
  Function4perl(&computePossibleMatroids,"cpm($,$,$,Vector<Set<Int> >,Vector<Set<Int> >)");
  Function4perl(&intersect_circuits,"isc(Matrix<Rational>,$)");
  
  
}}