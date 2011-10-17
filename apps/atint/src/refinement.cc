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
 
 This file contains all functionality for computing refinements of non-compatible cones
 */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
    
  using polymake::polytope::cdd_interface::solver;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  perl::Object complexify(Matrix<Rational> rays, Vector<Set<int> > max_cones, Vector<Integer> weights, bool uses_homog) {
    solver<Rational> sv;
    //First we compute all the H-representations of the cones
    dbgtrace << "Computing H-rep of cones" << endl;
    Vector<Matrix<Rational> > inequalities(max_cones.dim());
    Vector<Matrix<Rational> > equalities(max_cones.dim());
    for(int c = 0; c < max_cones.dim(); c++) {
      //To make computations come out right, we always have to add a zero in front
      std::pair<Matrix<Rational>, Matrix<Rational> > fanEqs =
		    sv.enumerate_facets(zero_vector<Rational>() | rays.minor(max_cones[c],All), Matrix<Rational>(0,rays.cols()+1), !uses_homog,false);
      inequalities[c] = fanEqs.first;
      equalities[c] = fanEqs.second;
    }
    //Compute the dimension
    int dimension = rays.cols() - equalities[0].rows() - (uses_homog? 1 : 0);
    
    //The following matrix marks pairs of cones (i,j) that don't need to be made compatible since they already are
    //More precisely: a value of true at position (i,j), i < j means that cone i and cone j are compatible. A value of false or
    // any value on or below the diagonal carries no information
    Matrix<bool> markedPairs(max_cones.dim(), max_cones.dim());
    
    //The following value indicates that new cones have been created through refinement 
    bool created;
    //The following values indicate that during an instance of the for(i,j)-loop below, 
    //the cone i and/or j has been removed and replaced by other cones. Note that here
    // assigning a different weight does count as "new"
//     bool replacedI = false;
//     bool replacedJ = false;
    
    //We go through all cone pairs i < j and check if we need to refine anything. As soon as we actually made some changes, we start all over
    do {
      created = false;
      for(int i = 0; i < max_cones.dim() -1 && !created; i++) {
	for(int j = i+1; j < max_cones.dim() && !created;j++) {
	    //replacedI = false; replacedJ = false;
	    dbgtrace << "Refining " << i << " and "<< j << endl;
	    //We only intersect non-marked pairs
	    if(!markedPairs[i][j]) {
	      //Compute an irredundant H-rep of the intersection of cones i and j
	      Matrix<Rational> isIneq = inequalities[i]/inequalities[j];
	      Matrix<Rational> isEq = equalities[i] / equalities[j];
	      Matrix<Rational> isMatrix = isIneq / isEq;
	      
	      std::pair<Bitset,Bitset> isection = sv.canonicalize(isIneq,isEq,1);
	      
	      int isDimension = isMatrix.cols() - isection.second.size() - (uses_homog? 2 : 1);
	      
	      //If the dimension is < 0 (this can only occur in the homog. case), then the two
	      //cells don't intersect at all and we don't need to refine
	      if(isDimension < 0) {
		continue;
	      }
	      
	      //This set contains the indices in isMatrix of equations coming from i
	      Set<int> IIndices = sequence(0,inequalities[i].rows()) +
				  sequence(isIneq.rows(),equalities[i].rows());
	      //We now put together the equations for the intersection, such that the equations from i
	      //come first
	      //Note: If the intersection is full-dimensional, we don't need to refine along the equalities
	      Matrix<Rational> hsEquations = isMatrix.minor(Set<int>(isection.first) * IIndices,All);
	      if(isDimension < dimension) hsEquations /= isMatrix.minor(Set<int>(isection.second) * IIndices,All);
	      int equationsFromI = hsEquations.rows();
	      hsEquations /= isMatrix.minor(Set<int>(isection.first) - IIndices,All);
	      if(isDimension < dimension) hsEquations /= isMatrix.minor(Set<int>(isection.second) - IIndices,All);
	      
	      //These variables remember the refinements of i and j (as indices in max_cones)
	      Set<int> newconesI, newconesJ;	
	       
	      for(int coneindex = i; coneindex <= j; coneindex += (j-i)) {
		//First we have to determine the relevant sign choices. We only want to change
		//the signs of equations NOT coming from the cone we currently refine
		//A sign choice's semantics is that x is in signChoices[s], iff the equation at index
		//s in relevantEquations has sign +1, otherwise it has sign -1
		Matrix<Rational> relevantEquations;
		if(coneindex == j) { 
		  relevantEquations = hsEquations.minor(sequence(0,equationsFromI),All);
		}
		else {
		  relevantEquations = hsEquations.minor(~sequence(0,equationsFromI),All);
		}
		Array<Set<int> > signChoices = 
		    pm::AllSubsets<Set<int> >(sequence(0,relevantEquations.rows()));
		
		dbgtrace << "Relevant equations: " << relevantEquations << endl;
		dbgtrace << "Cone index: " << coneindex << " of " << i << ", " << j << endl;
		for(int s = 0; s < signChoices.size(); s++) {
		  dbgtrace << "Sign choice: " << signChoices[s] << endl;
		  //If the intersection is full-dimensional and the signs are all +1's, then we
		  //only compute this intersection for coneindex == i
		  if(isDimension == dimension && signChoices[s].size() == relevantEquations.rows() && coneindex == j) {
		    dbgtrace << "Full-dimension intersection not for " << coneindex << endl;
		    continue;
		  }
		  
		  //Compute intersection of coneindex with current sign choice
		  
		  Matrix<Rational> refIneqs = relevantEquations.minor(signChoices[s],All);
		  refIneqs /= (-relevantEquations.minor(~signChoices[s],All));	
		  
		  Matrix<Rational> refineIneqs = inequalities[coneindex] / refIneqs;
		  try {
		    Matrix<Rational> ref = sv.enumerate_vertices(uses_homog? 
					  refineIneqs.minor(All,~scalar2set(0)) :
					  refineIneqs,
					  uses_homog? (equalities[coneindex].minor(All,~scalar2set(0))) :
					  equalities[coneindex],!uses_homog,true).first;
		    if(!uses_homog) {
		      ref = ref.minor(All,~scalar2set(0));
		    }
		    dbgtrace << "Intersection rays: " << ref << endl;
		    //If the refinement is full-dimensional, we get a new cone 
		    if(rank(ref) - (uses_homog? 1 :0 ) == dimension) {
		      dbgtrace << "is refinement" << endl;
		      //First we canonicalize the directional rays
		      for(int rw = 0; rw < ref.rows(); rw++) {
			if(!uses_homog || ref(rw,0) == 0) {
			    for(int cl = 0; cl < ref.cols();cl++) {
			      if(ref(rw,cl) != 0) {
				ref.row(rw) /= abs(ref(rw,cl));
				break;
			      }
			    }
			}
		      }
		      
		      //Add as new cone
		      //Go through all rays and check if they already exist
		      //Assign appropriate indices
		      Set<int> coneRays;
		      int noOfRays = rays.rows();
		      for(int nr = 0; nr < ref.rows(); nr++) {
			  for(int r = 0; r < noOfRays; r++) {
			    if(rays.row(r) == ref.row(nr)) {
			      coneRays += r;
			      break;
			    }
			    if(r == noOfRays-1) {
			      rays /= ref.row(nr);
			      int newrayindex = rays.rows()-1;
			      coneRays += newrayindex;
			    }
			  }
		      }
		      dbgtrace << "Cone rays are " << coneRays << endl;
		      
		      //We have to check for doubles: If i fulfills one of the relevant equations with
		      //equality, we can change its sign and still get the same refinement cone
		      bool isDouble = false;
		      for(Entire<Set<int> >::iterator existing = entire(coneindex == i? newconesI : newconesJ); !existing.at_end(); existing++) {
			if(max_cones[*existing] == coneRays) {
			  isDouble = true; break;
			}
		      }
		      if(isDouble) continue;
		      
		      //Add the cone
		      max_cones |= coneRays;
		      int newconeindex = max_cones.dim()-1;
		      if(coneindex == i) newconesI += newconeindex;
		      else newconesJ += newconeindex;
		      dbgtrace << "Cone gets index " << newconeindex << endl;
		      
		      //Weight is the weight of coneindex (or the sum of both, if s is everything)
		      weights |= (isDimension < dimension || 
				  signChoices[s].size() < relevantEquations.rows()? 
				  weights[coneindex] : weights[i] + weights[j]);
		      
		      dbgtrace << "Cone weight is " << weights[weights.dim()-1] << endl;
		      
		      inequalities |= (inequalities[coneindex] / refIneqs);
		      equalities |= equalities[coneindex];
		      markedPairs |= Vector<bool>(markedPairs.rows());
		      markedPairs /= Vector<bool>(markedPairs.cols());		
		    } //END if refDimension = dimension
		  }
		  catch(...) {
		    dbgtrace << "Empty refinement intersection - can be ignored" << endl;
		  }
		}//END iterate sign choices
	      }//END iterate i and j

	      //Now if we did create anything new, we remove i and/or j
	      //(depending on whether any proper refinement of each cone was created)
	      //We only set 'created' to true, if at least one of i and j has 
	      //been replaced by at least two cones.
	      //If both have been replaced by a single cone (i.e. by themselves),
	      //we throw out the new cones instead of the old ones
	      //(Except, if the intersection is full-dimensional. In this case i and j agree
	      // and they are replaced by the single new cone with weight the sum of their weights)
	      if(newconesI.size() + newconesJ.size() > 0) {
		dbgtrace << "Cleaning up indices and markings" << endl;
		
		dbgtrace << "Created cones " << newconesI << " in " << i << " and " << newconesJ << " in " << j << endl;
		
		created = true;//newconesI.size() >= 2 || newconesJ.size() >= 2;
		//Will contain the indices of cones to be marked compatible
		Vector<int> markIndices;
		//Will contain the indices of cones to be removed
		Set<int> removedIndices;
		
		
		//Find out indices to be marked and removed:
		
		//If the cones intersect in full dimension, then we keep the new cones in any case
		//and discard the old ones
		if(isDimension == dimension) {
		  removedIndices += i; removedIndices += j;
		  markIndices |= Vector<int>(newconesI);
		  markIndices |= Vector<int>(newconesJ);
		  //replacedI = replacedJ = true;
		}
		else {
		  //Otherwise we discard the old cones iff they have been refined by at least 2 cones
		  if(newconesI.size() >= 2) {
		    markIndices |= Vector<int>(newconesI);
		    removedIndices += i;
		    //replacedI = true;
		  }
		  else {
		    markIndices |= i;
		    removedIndices += newconesI;
		  }		  
		  
		  if(newconesJ.size() >= 2) {
		    markIndices |= Vector<int>(newconesJ);
		    removedIndices += j;
		    //replacedJ = true;
		  }
		  else {
		    markIndices |= j;
		    removedIndices += newconesJ;
		  }
		}
		
// 		for(int coneindex = i; coneindex <= j; coneindex += (j-i)) {
// 		  Set<int> st = coneindex == i? newconesI : newconesJ;
// 		  for(Entire<Set<int> >::iterator ni = entire(st); !ni.at_end(); ni++) {
// 		    if(max_cones[*ni].contains(3) && max_cones[*ni].contains(8) && !removedIndices.contains(*ni)) {
// 		      dbglog << "Keeping 3,8 cone " << *ni << endl;
// 		      dbglog << "Created in intersecting " << i << "," << j << endl;
// 		      dbglog << "Cone of " << coneindex << endl;
// 		      dbglog << "cones in i: " << newconesI << endl;
// 		      dbglog << "cones in j: " << newconesJ << endl;
// 		      dbglog << "Created: " << created << endl;
// 		      for(int k = 0; k < *ni; k++) {
// 			if(max_cones[k] == max_cones[*ni]) {
// 			    dbglog << "Same as cone " << k << endl;
// 			    dbglog << "Compatible: " << markedPairs(k,*ni) << endl;
// 			}
// 		      }
// 		    }
// 		  }
// 		}
	
		dbgtrace << "Marking indices " << markIndices << " and removing " << removedIndices << endl;
		
		for(int v = 0; v < markIndices.dim()-1; v++) {
		  for(int w = v+1; w < markIndices.dim(); w++) {
		    markedPairs(markIndices[v],markIndices[w]) = true;
		  }
		}
		
		dbgtrace << "Done marking" << endl;
		
		//Now remove them
		max_cones = max_cones.slice(~removedIndices);
		weights = weights.slice(~removedIndices);
		inequalities = inequalities.slice(~removedIndices);
		equalities = equalities.slice(~removedIndices);
		markedPairs = markedPairs.minor(~removedIndices,~removedIndices);
		dbgtrace << "Done." << endl;
		dbgtrace << "Maximal cones: " << endl;
		dbgtrace << max_cones.dim() << endl;
	      }
	      
	      
	    }//END refine i and j
	}//END iterate j
      }//END iterate i
//       for(int k = 0; k < max_cones.dim();k++) {
// 	 for(int l = k+1; l < max_cones.dim(); l++) {
// 	    if(max_cones[k] == max_cones[l] && max_cones[k].contains(3) && max_cones[k].contains(8)) {
// 	      dbglog << "Cones " << k << " and " << l << " agree. Compatible: " << markedPairs(k,l) << endl;
// 	      dbglog << "Restarting: " << created << endl;
// 	    }
// 	 }
//       }
    } while(created);
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << max_cones;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
    return result;
  }
  
  
//   perl::Object complexify(Matrix<Rational> rays, Vector<Set<int> > max_cones, Vector<Integer> weights, bool uses_homog) {
//       solver<Rational> sv;
//       //First we compute all the H-representations of the cones
//       dbgtrace << "Computing H-rep of cones" << endl;
//       Vector<Matrix<Rational> > inequalities(max_cones.dim());
//       Vector<Matrix<Rational> > equalities(max_cones.dim());
//       for(int c = 0; c < max_cones.dim(); c++) {
// 	std::pair<Matrix<Rational>, Matrix<Rational> > fanEqs =
// 		      solver<Rational>().enumerate_facets(rays.minor(max_cones[c],All), Matrix<Rational>(0,rays.cols()), true);
// 	inequalities[c] = fanEqs.first;
// 	equalities[c] = fanEqs.second;	
//       }
//       //Compute the dimension
//       int dimension = rays.cols() - equalities[0].rows();
//       dbglog << "Dimension: " << dimension << endl;
//       //The following matrix marks pairs of cones (i,j) that don't need to be made compatible since they already are
//       //More precisely: a value of true at position (i,j), i < j means that cone i and cone j are compatible. A value of false or
//       // any value on or below the diagonal carries no information
//       Matrix<bool> markedPairs(max_cones.dim(), max_cones.dim());
//       
//       //The following value indicates that new cones have been created through refinement 
//       bool created;
//       dbgtrace << "Starting intersection... " << endl;
//             
//       do {
// 	created = false;
// 	dbgtrace << "Starting refinement from scratch " << endl;
// 	dbgtrace << "Maximal cones: " << max_cones.dim() << endl;
// 	//We go through all cone pairs i < j and check if we need to refine anything. As soon as we actually made some changes, we start all over
// 	for(int i = 0; i < max_cones.dim()-1 && !created; i++) {
// 	    for(int j = i+1; j < max_cones.dim() && !created; j++) {
// 		dbgtrace << "Are " << i << " and " << j << " marked: " << endl;
// 		dbgtrace << markedPairs[i][j] << endl;
// 		//We only intersect cones that may be incompatible
// 		if(!markedPairs[i][j]) {
// 		      dbgtrace << "Refining cones " << i << " and "<< j << endl;
// 		      //Compute an irredundant H-representation of the intersection:
// 		      Matrix<Rational> isIneq = inequalities[i]/inequalities[j];
// 		      Matrix<Rational> isEq = equalities[i] / equalities[j];
// 		      Matrix<Rational> isMatrix = isIneq / isEq;
// 		      
// 		      std::pair<Bitset,Bitset> isection = sv.canonicalize(isIneq,isEq,1);
// 		      int isDimension = isMatrix.cols() - isection.second.size() ;
// 		      dbgtrace << "intersection dimension: " << isDimension << endl;
// 		      
// 		      //This set contains the indices in isMatrix of equations coming from i
// 		      Set<int> IIndices = sequence(0,inequalities[i].rows()) +
// 					  sequence(isIneq.rows(),equalities[i].rows());
// 		      //We now put together the equations for the intersection, such that the equations from i
// 		      //come first
// 		      //Note: If the intersection is full-dimensional, we don't need to refine along the equalities
// 		      Matrix<Rational> hsEquations = isMatrix.minor(Set<int>(isection.first) * IIndices,All);
// 		      if(isDimension < dimension) hsEquations /= isMatrix.minor(Set<int>(isection.second) * IIndices,All);
// 		      int equationsFromI = hsEquations.rows();
// 		      hsEquations /= isMatrix.minor(Set<int>(isection.first) - IIndices,All);
// 		      if(isDimension < dimension) hsEquations /= isMatrix.minor(Set<int>(isection.second) - IIndices,All);
// 		      
// 		      dbgtrace << "Computed intersection " << endl;
// 		      dbgtrace << "Has equations: " << hsEquations << endl;
// 		      dbgtrace << "Coming from i: " << equationsFromI << endl;
// 		      
// 		      //These variables remember the refinements of i and j (as indices in max_cones)
// 		      Set<int> newconesI, newconesJ;	
// 		      //These remember the row indices of the newly added rays
// 		      //Set<int> newraysI, newraysJ;
// 		      
// 		      for(int coneindex = i; coneindex != j; coneindex = j) {
// 			  //First we have to determine the relevant sign choices. We only want to change
// 			  //the signs of equations NOT coming from the cone we currently refine
// 			  //A sign choice's semantics is that x is in signChoices[s], iff the equation at index
// 			  //s in relevantEquations has sign +1, otherwise it has sign -1
// 			  Matrix<Rational> relevantEquations;
// 			  if(coneindex == j) { 
// 			    relevantEquations = hsEquations.minor(sequence(0,equationsFromI),All);
// 			  }
// 			  else {
// 			    relevantEquations = hsEquations.minor(~sequence(0,equationsFromI),All);
// 			  }
// 			  Array<Set<int> > signChoices = 
// 			      pm::AllSubsets<Set<int> >(sequence(0,relevantEquations.rows()));
// 			  
// 			  for(int s = 0; s < signChoices.size(); s++) {
// 			    //If the intersection is full-dimensional and the signs are all +1's, then we
// 			    //only compute this intersection for coneindex == i
// 			    if(isDimension == dimension && signChoices.size() == relevantEquations.rows() && coneindex == j) continue;
// 			    
// 			    //Compute intersection of coneindex with current sign choice
// 			    
// 			    Matrix<Rational> refIneqs = relevantEquations.minor(signChoices[s],All);
// 			    refIneqs /= (-relevantEquations.minor(~signChoices[s],All));
// 			    try {
// 			      Matrix<Rational> ref = sv.enumerate_vertices(inequalities[coneindex] / refIneqs,equalities[coneindex],true).first;
// 			      //Check the dimension of the intersection
// 			      if(rank(ref) == dimension) {
// 				  dbgtrace << "Found refinement of " << coneindex << endl;
// 				  dbgtrace << "Rays are " << ref << endl;
// 				  //Add as new cone
// 				  //Go through all rays and check if they already exist
// 				  Set<int> coneRays;
// 				  int noOfRays = rays.rows();
// 				  for(int nr = 0; nr < ref.rows(); nr++) {
// 				      for(int r = 0; r < noOfRays; r++) {
// 					if(rays.row(r) == ref.row(nr)) {
// 					  coneRays += r;
// 					  break;
// 					}
// 					if(r == noOfRays-1) {
// 					  rays /= ref.row(nr);
// 					  int newrayindex = rays.rows()-1;
// 					  coneRays += newrayindex;
// 					  //(coneindex == i? newraysI : newraysJ) += newrayindex;
// 					}
// 				      }
// 				  }
// 				  dbgtrace << "Added rays" << endl;
// 				  //Add the cone and remember its index
// 				  max_cones |= coneRays;
// 				  int newconeindex = max_cones.dim()-1;
// 				  if(coneindex == i) newconesI += newconeindex;
// 				  else newconesJ += newconeindex;
// 				  //Weight is the weight of coneindex (or the sum of both, if s is everything)
// 				  weights |= (isDimension < dimension || 
// 					      signChoices[s].size() < relevantEquations.rows()? 
// 					      weights[coneindex] : weights[i] + weights[j]);
// 				  
// 				  inequalities |= (inequalities[coneindex] / refIneqs);
// 				  equalities |= equalities[coneindex];
// 				  markedPairs |= zero_vector<bool>();
// 				  markedPairs /= zero_vector<bool>();
// 				  dbgtrace << "Inserted new cone with weight " << weights[weights.dim()-1] << endl;
// 			      }
// 			    }
// 			    catch(...) { //Is called whenever the intersection is empty
// 			      continue;
// 			    }
// 			    
// 			  }
// 		      }
// 		      
// 		      
// 		      //Now if we did create anything new, we remove i and/or j
// 		      //(depending on whether any proper refinement of each cone was created)
// 		      //We only set 'created' to true, if at least one of i and j has 
// 		      //been replaced by at least two cones.
// 		      //If both have been replaced by a single cone (i.e. by themselves),
// 		      //we throw out the new cones instead of the old ones
// 		      //(Except, if the intersection is full-dimensional. In this case i and j agree
// 		      // and they are replaced by the single new cone with weight the sum of their weights)
// 		      if(newconesI.size() + newconesJ.size() > 0) {
// 			dbgtrace << "Cleaning up indices and markings" << endl;
// 			created = newconesI.size() >= 2 || newconesJ.size() >= 2;
// 			//Will contain the indices of cones to be marked compatible
// 			Vector<int> markIndices;
// 			//Will contain the indices of cones to be removed
// 			Set<int> removedIndices;
// 			
// 			
// 			//Find out indices to be marked and removed:
// 			if(newconesI.size() >= 2) {
// 			  markIndices |= Vector<int>(newconesI);
// 			  removedIndices += i;
// 			}
// 			else {
// 			  markIndices |= i;
// 			  if(isDimension == dimension && newconesJ.size() < 2) removedIndices += i;
// 			  else removedIndices += newconesI;
// 			}
// 			if(newconesJ.size() >= 2) {
// 			  markIndices |= Vector<int>(newconesJ);
// 			  removedIndices += j;
// 			}
// 			else {
// 			  markIndices |= j;
// 			  if(isDimension == dimension && newconesI.size() < 2) removedIndices += j;
// 			  else removedIndices += newconesJ;
// 			}
// 			
// 			dbgtrace << "Marking indices " << markIndices << " and removing " << removedIndices << endl;
// 			
// 			for(int v = 0; v < markIndices.dim()-1; v++) {
// 			  for(int w = 0; w < markIndices.dim() -1; w++) {
// 			    markedPairs(markIndices[v],markIndices[w]) = true;
// 			  }
// 			}
// 			//Now remove them
// 			max_cones = max_cones.slice(~removedIndices);
// 			weights = weights.slice(~removedIndices);
// 			inequalities = inequalities.slice(~removedIndices);
// 			equalities = equalities.slice(~removedIndices);
// 			markedPairs = markedPairs.minor(~removedIndices,~removedIndices);
// 			dbgtrace << "Done." << endl;
// 			dbgtrace << "Maximal cones: " << endl;
// 			dbgtrace << max_cones.dim() << endl;
// 		      }
// 		}
// 	    }
// 	}
//       } while(created);
//       
//       dbgtrace << "Done with refinement." << endl;
//       
//       perl::Object result("WeightedComplex");
// 	result.take("RAYS") << rays;
// 	result.take("MAXIMAL_CONES") << max_cones;
// 	result.take("TROPICAL_WEIGHTS") << weights;
// 	result.take("USES_HOMOGENEOUS_C") << uses_homog;
// 	
//       return result;
//       
//   }
  

  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&complexify,"complexify(Matrix<Rational>,IncidenceMatrix,Vector<Integer>, $)");
  //Function4perl(&testcone, "testcone(polytope::Polytope)");
  
}}