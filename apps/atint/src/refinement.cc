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
  
  //Documentation see header 
  perl::Object complexify(Matrix<Rational> rays, IncidenceMatrix<> max_cones, Vector<Integer> weights, bool uses_homog) {
      Vector<Set<int> > cones(max_cones.rows());
	for(int i = 0; i < max_cones.rows(); i++) { cones[i] = max_cones.row(i);}
      dbgtrace << "Precomputing (in)equalities" << endl;
      //We precompute the defining (in)equalities of all cones
      Vector<Matrix<Rational> > inequalities(cones.dim());
      Vector<Matrix<Rational> > equalities(cones.dim());
      for(int c = 0; c < cones.dim(); c++) {
	std::pair<Matrix<Rational>, Matrix<Rational> > fanEqs =
		      solver<Rational>().enumerate_facets(zero_vector<Rational>() | rays.minor(cones[c],All), Matrix<Rational>(0,rays.cols()+1), true);
	inequalities[c] = fanEqs.first;
	equalities[c] = fanEqs.second;	
      }
      
      dbgtrace << "Intersecting cones" << endl;
    
      //We go through all pairs of cones and refine them along each other. As long as we create
      //a new cone in the process, we restart this process from the beginning
      bool created = true;
      while(created) {
	created = false;
	for(int i = 0; i < cones.dim() && !created; i++) {
	    for(int j = 0; j < cones.dim() && !created; j++) {
	      if(i == j) continue;
	      //We refine cone i along the equations of cone j
	      Matrix<Rational> intersections = inequalities[j] / equalities[j];
	      for(Entire< Rows< Matrix<Rational> > >::iterator r = entire(rows(intersections)); !r.at_end() && !created; r++) {
		try {
		  dbgtrace << "Intersecting cone " << i << " along equation " << *r << " of cone " << j << endl;
		  std::pair<Matrix<Rational>, Matrix<Rational> > onehalf = 
		    solver<Rational>().enumerate_vertices(inequalities[i] / *r, equalities[i],true);
		  std::pair<Matrix<Rational>, Matrix<Rational> > otherhalf =
		    solver<Rational>().enumerate_vertices(inequalities[i] / (-*r), equalities[i],true);
		  //If this is a relevant intersection then both cones have full dimension
		  // (<=> they have the same dimension)
		  if(rank(onehalf.first) == rank(otherhalf.first)) {
		    dbgtrace << "Has correct dimension, creating new halfs" << endl;
		    //Normalize ray matrices
		    Matrix<Rational> firstMatrix = onehalf.first.minor(All,~scalar2set(0));
		    Matrix<Rational> secondMatrix = otherhalf.first.minor(All,~scalar2set(0));
		    if(uses_homog) {
		      for(int k = 0; k < firstMatrix.rows(); k ++) {
			if(firstMatrix(k,0) != 0) firstMatrix.row(k) /= firstMatrix(k,0);
		      }
		      for(int k = 0; k < secondMatrix.rows(); k++) {
			if(secondMatrix(k,0) != 0) secondMatrix.row(k) /= secondMatrix(k,0);
		      }
		    }
		    dbgtrace << "First half\n" << firstMatrix<< endl;
		    dbgtrace << "Second half\n" << secondMatrix << endl;
		    created = true;
		    //Then we compute the ray indices of the new cones
		    Set<int> onehalfIndices, otherhalfIndices;
		    for(int x = 1; x <= 2; x++) {
		      Matrix<Rational> checkRays = x==1? firstMatrix: secondMatrix;
		      for(int r = 0; r < checkRays.rows(); r++) {
			int foundIndex = -1;
			for(int ray = 0; ray < rays.rows(); ray++) {
			    if(rays.row(ray) == checkRays.row(r)) {
			      foundIndex = ray; break;
			    }
			}
			if(foundIndex < 0) {
			    rays /= checkRays.row(r);
			    foundIndex = rays.rows()-1;
			}
			(x==1? onehalfIndices : otherhalfIndices) += foundIndex;
		      }
		    }
		    cones |= onehalfIndices; cones |= otherhalfIndices;
		    dbgtrace << "Computed new indices" << endl;
		    dbgtrace << "Index sets are: " << onehalfIndices << ", " << otherhalfIndices << endl;
		    //We also insert the inequalities / equalities of the new cones
		    // and its new weight (which is just the weight of cone i)
		    inequalities |= (inequalities[i] / *r);
		    inequalities |= (inequalities[i] / (-*r));
		    equalities |= (equalities[i]);
		    equalities |= (equalities[i]);
		    weights |= weights[i];
		    weights |= weights[i];
		    dbgtrace << "Added new cones " << endl;
		    //Finally we remove cone i
		    cones = cones.slice(~scalar2set(i));
		    inequalities = inequalities.slice(~scalar2set(i));
		    equalities = equalities.slice(~scalar2set(i));
		    weights = weights.slice(~scalar2set(i));
		    dbgtrace << "Replaced cone " << endl;
		  }
		}
		//If an error occurs, then one of the halfs is empty and we dont want it anyway
		catch(...) { 
		  //--
		}
	      } //END iterate equations
	    } //END iterate j
	} //END iterate i
	
      } //END while cones are created
      
      dbgtrace << "Checking for doubles" << endl;
      
      //At this point we more or less have a weighted complex, except that some cones might 
      //be identical.
      for(int i = 0; i < cones.dim()-1; i++) {
	for(int j = i+1; j < cones.dim(); j++) {
	  //Check if both cone sets are equal
	  if((cones[i] * cones[j]).size() == cones[i].size()) {
	    //Add j's weight to i and remove j
	    weights[i] += weights[j];
	    cones = cones.slice(~scalar2set(j));
	    weights = weights.slice(~scalar2set(j));
	    j--;
	  } 
	}
      }
      
      //Now return the resulting complex
      perl::Object result("WeightedComplex");
	result.take("RAYS") << rays; 
	result.take("MAXIMAL_CONES") << cones;
	result.take("TROPICAL_WEIGHTS") << weights;
      return result;
  }

  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&complexify,"complexify(Matrix<Rational>,IncidenceMatrix,Vector<Integer>)");
  
}}