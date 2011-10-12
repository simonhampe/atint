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
  
  //using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
    
  using polymake::polytope::cdd_interface::solver;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  perl::Object complexify(Matrix<Rational> rays, Vector<Set<int> > max_cones, Vector<Integer> weights, bool uses_homog) {
      solver<Rational> sv;
      //First we compute all the H-representations of the cones
      dbgtrace << "Computing H-rep of cones" << endl;
      Vector<Matrix<Rational> > inequalities(max_cones.dim());
      Vector<Matrix<Rational> > equalities(max_cones.dim());
      for(int c = 0; c < max_cones.dim(); c++) {
	std::pair<Matrix<Rational>, Matrix<Rational> > fanEqs =
		      solver<Rational>().enumerate_facets(zero_vector<Rational>() | rays.minor(max_cones[c],All), Matrix<Rational>(0,rays.cols()+1), true);
	inequalities[c] = fanEqs.first;
	equalities[c] = fanEqs.second;	
      }
      //Compute the dimension
      int dimension = rank(rays.minor(max_cones[0],All));
      //The following matrix marks pairs of cones (i,j) that don't need to be made compatible since they already are
      //More precisely: a value of true at position (i,j), i < j means that cone i and cone j are compatible. A value of false or
      // any value on or below the diagonal carries no information
      Matrix<bool> markedPairs(max_cones.dim(), max_cones.dim());
      //The following value indicates that new cones have been created through refinement 
      bool created = false;
      dbgtrace << "Starting intersection... " << endl;
      do {
	created = false;
	//We go through all cone pairs i < j and check if we need to refine anything. As soon as we actually made some changes, we start all over
	for(int i = 0; i < max_cones.dim()-1 && !created; i++) {
	    for(int j = i+1; j < max_cones.dim() && !created; j++) {
		//We only intersect cones that may be incompatible
		if(!markedPairs[i][j]) {
		      dbgtrace << "Refining cones " << i << " and "<< j << endl;
		      //Compute an irredundant H-representation of the intersection:
		      Matrix<Rational> isIneq = inequalities[i]/inequalities[j];
		      Matrix<Rational> isEq = equalities[i] / equalities[j];
		      Matrix<Rational> isMatrix = isIneq / isEq;
		      std::pair<Bitset,Bitset> isection = sv.canonicalize(isIneq,isEq,1);
		      Matrix<Rational> hsEquations = isMatrix.minor(isection.first,All) / isMatrix.minor(isection.second,All);
		      dbgtrace << "Computed intersection " << endl;
		      dbgtrace << "Has equations: " << hsEquations << endl;
		      //Go through all sign choices for the equations and compute the corr. refinement of i and j
		      
		      //This iterates all sign choices: x is in a set of this list, iff the corr
		      // sign is +1
		      Array<Set<int> > signChoices = pm::AllSubsets<Set<int> >(sequence(0,hsEquations.rows()));
		      dbgtrace << "Going through all sign choices " <<  signChoices << endl;
		      for(int s = 0; s < signChoices.size(); s++) {
			//If s is not the whole set (which is simply the intersectio of i and j)
			//then we compute the refinements of both i and j, otherwise one is sufficient
			for(int coneindex = i; coneindex != j && signChoices[s].size() < hsEquations.rows(); coneindex = j) {
			    Matrix<Rational> refIneqs = hsEquations.minor(signChoices[s],All);
			    refIneqs /= (-hsEquations.minor(~signChoices[s],All));
			    Matrix<Rational> ref = sv.enumerate_vertices(inequalities[coneindex] / refIneqs,equalities[coneindex],true).first;
			    //Check the dimension of the intersection
			    if(rank(ref) == dimension) {
				dbgtrace << "Found refinement" << endl;
				created = true;
				//Add as new cone
				//Go through all rays and check if they already exist
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
					coneRays += (rays.rows()-1);
				      }
				    }
				}
				max_cones |= coneRays;
				//Weight is the weight of coneindex (or the sum of both, if s is everything)
				weights |= (signChoices[s].size() < hsEquations.rows()? weights[coneindex] : weights[i] + weights[j]);
				
				inequalities |= refIneqs;
				equalities |= equalities[coneindex];
			    }			    
			}
		      } //End refine i and j
		      
		      //Now if we did create anything new, we remove i and j
		      if(created) {
			dbgtrace << "Removing cones " << i << " and " << j << endl;
			Set<int> removedIndices; removedIndices += i; removedIndices += j;
			max_cones = max_cones.slice(~removedIndices);
			weights = weights.slice(~removedIndices);
		      }
		}
	    }
	}
      } while(created);
      
      perl::Object result("WeightedComplex");
	result.take("RAYS") << rays;
	result.take("MAXIMAL_CONES") << max_cones;
	result.take("TROPICAL_WEIGHTS") << weights;
	result.take("USES_HOMOGENEOUS_C") << uses_homog;
	
      return result;
      
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
//   //Documentation see header 
//   perl::Object complexify(Matrix<Rational> rays, IncidenceMatrix<> max_cones, Vector<Integer> weights, bool uses_homog) {
//       Vector<Set<int> > cones(max_cones.rows());
// 	for(int i = 0; i < max_cones.rows(); i++) { cones[i] = max_cones.row(i);}
//       dbgtrace << "Precomputing (in)equalities" << endl;
//       //We precompute the defining (in)equalities of all cones
//       Vector<Matrix<Rational> > inequalities(cones.dim());
//       Vector<Matrix<Rational> > equalities(cones.dim());
//       for(int c = 0; c < cones.dim(); c++) {
// 	std::pair<Matrix<Rational>, Matrix<Rational> > fanEqs =
// 		      solver<Rational>().enumerate_facets(zero_vector<Rational>() | rays.minor(cones[c],All), Matrix<Rational>(0,rays.cols()+1), true);
// 	inequalities[c] = fanEqs.first;
// 	equalities[c] = fanEqs.second;	
//       }
//       
//       dbgtrace << "Intersecting cones" << endl;
//     
//       //We go through all pairs of cones and refine them along each other. As long as we create
//       //a new cone in the process, we restart this process from the beginning
//       bool created = true;
//       while(created) {
// 	created = false;
// 	for(int i = 0; i < cones.dim() && !created; i++) {
// 	    for(int j = 0; j < cones.dim() && !created; j++) {
// 	      if(i == j) continue;
// 	      //We refine cone i along the equations of cone j
// 	      Matrix<Rational> intersections = inequalities[j] / equalities[j];
// 	      for(Entire< Rows< Matrix<Rational> > >::iterator r = entire(rows(intersections)); !r.at_end() && !created; r++) {
// 		try {
// 		  dbgtrace << "Intersecting cone " << i << " along equation " << *r << " of cone " << j << endl;
// 		  std::pair<Matrix<Rational>, Matrix<Rational> > onehalf = 
// 		    solver<Rational>().enumerate_vertices(inequalities[i] / *r, equalities[i],true);
// 		  std::pair<Matrix<Rational>, Matrix<Rational> > otherhalf =
// 		    solver<Rational>().enumerate_vertices(inequalities[i] / (-*r), equalities[i],true);
// 		  //If this is a relevant intersection then both cones have full dimension
// 		  // (<=> they have the same dimension)
// 		  if(rank(onehalf.first) == rank(otherhalf.first)) {
// 		    dbgtrace << "Has correct dimension, creating new halfs" << endl;
// 		    //Normalize ray matrices
// 		    Matrix<Rational> firstMatrix = onehalf.first.minor(All,~scalar2set(0));
// 		    Matrix<Rational> secondMatrix = otherhalf.first.minor(All,~scalar2set(0));
// 		    if(uses_homog) {
// 		      for(int k = 0; k < firstMatrix.rows(); k ++) {
// 			if(firstMatrix(k,0) != 0) firstMatrix.row(k) /= firstMatrix(k,0);
// 		      }
// 		      for(int k = 0; k < secondMatrix.rows(); k++) {
// 			if(secondMatrix(k,0) != 0) secondMatrix.row(k) /= secondMatrix(k,0);
// 		      }
// 		    }
// 		    dbgtrace << "First half\n" << firstMatrix<< endl;
// 		    dbgtrace << "Second half\n" << secondMatrix << endl;
// 		    created = true;
// 		    //Then we compute the ray indices of the new cones
// 		    Set<int> onehalfIndices, otherhalfIndices;
// 		    for(int x = 1; x <= 2; x++) {
// 		      Matrix<Rational> checkRays = x==1? firstMatrix: secondMatrix;
// 		      for(int r = 0; r < checkRays.rows(); r++) {
// 			int foundIndex = -1;
// 			for(int ray = 0; ray < rays.rows(); ray++) {
// 			    if(rays.row(ray) == checkRays.row(r)) {
// 			      foundIndex = ray; break;
// 			    }
// 			}
// 			if(foundIndex < 0) {
// 			    rays /= checkRays.row(r);
// 			    foundIndex = rays.rows()-1;
// 			}
// 			(x==1? onehalfIndices : otherhalfIndices) += foundIndex;
// 		      }
// 		    }
// 		    cones |= onehalfIndices; cones |= otherhalfIndices;
// 		    dbgtrace << "Computed new indices" << endl;
// 		    dbgtrace << "Index sets are: " << onehalfIndices << ", " << otherhalfIndices << endl;
// 		    //We also insert the inequalities / equalities of the new cones
// 		    // and its new weight (which is just the weight of cone i)
// 		    inequalities |= (inequalities[i] / *r);
// 		    inequalities |= (inequalities[i] / (-*r));
// 		    equalities |= (equalities[i]);
// 		    equalities |= (equalities[i]);
// 		    weights |= weights[i];
// 		    weights |= weights[i];
// 		    dbgtrace << "Added new cones " << endl;
// 		    //Finally we remove cone i
// 		    cones = cones.slice(~scalar2set(i));
// 		    inequalities = inequalities.slice(~scalar2set(i));
// 		    equalities = equalities.slice(~scalar2set(i));
// 		    weights = weights.slice(~scalar2set(i));
// 		    dbgtrace << "Replaced cone " << endl;
// 		  }
// 		}
// 		//If an error occurs, then one of the halfs is empty and we dont want it anyway
// 		catch(...) { 
// 		  //--
// 		}
// 	      } //END iterate equations
// 	    } //END iterate j
// 	} //END iterate i
// 	
//       } //END while cones are created
//       
//       dbgtrace << "Checking for doubles" << endl;
//       
//       //At this point we more or less have a weighted complex, except that some cones might 
//       //be identical.
//       for(int i = 0; i < cones.dim()-1; i++) {
// 	for(int j = i+1; j < cones.dim(); j++) {
// 	  //Check if both cone sets are equal
// 	  if((cones[i] * cones[j]).size() == cones[i].size()) {
// 	    //Add j's weight to i and remove j
// 	    weights[i] += weights[j];
// 	    cones = cones.slice(~scalar2set(j));
// 	    weights = weights.slice(~scalar2set(j));
// 	    j--;
// 	  } 
// 	}
//       }
//       
//       //Now return the resulting complex
//       perl::Object result("WeightedComplex");
// 	result.take("RAYS") << rays; 
// 	result.take("MAXIMAL_CONES") << cones;
// 	result.take("TROPICAL_WEIGHTS") << weights;
//       return result;
//   }

  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&complexify,"complexify(Matrix<Rational>,IncidenceMatrix,Vector<Integer>)");
  //Function4perl(&testcone, "testcone(polytope::Polytope)");
  
}}