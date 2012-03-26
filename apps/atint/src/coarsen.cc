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
 * This file contains a function to compute the coarsest polyhedral structure of a
 * tropical variety
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  using polymake::polytope::cdd_interface::solver;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object coarsen(perl::Object complex) {
    //Extract values
    Matrix<Rational> rays = complex.give("RAYS");
    Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
    IncidenceMatrix<> maximalCones = complex.give("MAXIMAL_CONES");
    IncidenceMatrix<> codimOneCones = complex.give("CODIM_1_FACES");
    IncidenceMatrix<> codimInMaximal = complex.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximalOverCodim = T(codimInMaximal);
    int noOfCones = maximalCones.rows();
    
    //If the fan has no rays, it has only lineality space and we return the complex itself
    if(rays.rows() == 0)  return complex;
    
    //Compute equivalence classes of maximal cones
    Vector<Set<int> > equivalenceClasses;
    Vector<bool> hasBeenAdded(noOfCones); //contains whether a cone has been added to an equivalence class
//     Map<int, int> conesInClasses; //Maps cone indices to the index of their class in equivalenceClasses
    
    //dbgtrace << "Computing equivalence classes" << endl;
    
    for(int mc = 0; mc < noOfCones; mc++) {
      if(!hasBeenAdded[mc]) {
	    Set<int> newset; newset = newset + mc;
	    //Do a breadth-first search for all other cones in the component
	    std::list<int> queue;
	      queue.push_back(mc);
	      //Semantics: Elements in that queue have been added but their neighbours might not
	    while(queue.size() != 0) {
	      //Take the first element and find its neighbours
	      int node = queue.front(); 
		queue.pop_front();
	      //dbgtrace << "Checking node " << node << endl;
	      Set<int> cdset = maximalOverCodim.row(node);
	      for(Entire<Set<int> >::iterator cd = entire(cdset); !cd.at_end(); cd++) {
		  Set<int> otherMaximals = codimInMaximal.row(*cd) - node;
		  //We are only interested in codim-one-faces that have exactly two adjacent maximal cones
		  if(otherMaximals.size() == 1) {
		    int othermc = *(otherMaximals.begin());
		    if(!hasBeenAdded[othermc]) {
		      //Now add the cone to the equivalence class of mc
		      newset += othermc;
		      queue.push_back(othermc);
		      hasBeenAdded[othermc] = true;
		    }
// 		    conesInClasses[othermc] = equivalenceClasses.dim();
		  }
	      }		
	    }//End iterate queue
	    equivalenceClasses |= newset;
// 	    conesInClasses[mc] = equivalenceClasses.dim();
      }	
    } //END iterate maximal cones
    
    //dbgtrace << "Equivalence classes: " << equivalenceClasses << endl;
    
    //Now compute the new cones as unions of the cones in each equivalence class
    Matrix<Rational> newlin; bool newlin_computed = false;
    Vector<Set<int> > newcones;
    Vector<Integer> newweights;
    solver<Rational> sv;
    Matrix<Rational> complete_matrix = rays / linspace;
    Set<int> used_rays;
    
    for(int cl = 0; cl < equivalenceClasses.dim(); cl++) {
      Matrix<Rational> union_rays(0,rays.cols());
      Set<int> conesInClass = equivalenceClasses[cl];
      Vector<int> union_ray_list;
      for(Entire<Set<int> >::iterator mc = entire(conesInClass); !mc.at_end(); mc++) {
	union_rays /= rays.minor(maximalCones.row(*mc),All);
	union_ray_list |= Vector<int>(maximalCones.row(*mc));
      }
      std::pair<Bitset,Bitset> union_cone = sv.canonicalize(union_rays, linspace,0);
      //dbgtrace << "Class " << cl << ": " << union_cone << endl;
      
      //Compute lineality if it hasn't been computed yet
      if(!newlin_computed) {
	newlin = (union_rays/linspace).minor(union_cone.second,All);
	newlin_computed = true;
      }
      
      //Convert indices of rays in union_rays to indices in rays
      Set<int> ray_set(union_ray_list.slice(union_cone.first));
      
      newcones |= ray_set;
      used_rays += ray_set;
      newweights |= weights[*(conesInClass.begin())];      
    }
    IncidenceMatrix<> newcone_matrix(newcones);
    newcone_matrix = newcone_matrix.minor(All,used_rays);
    complete_matrix = complete_matrix.minor(used_rays,All);
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << complete_matrix;
      result.take("MAXIMAL_CONES") << newcone_matrix;
      result.take("LINEALITY_SPACE") << newlin;
      result.take("TROPICAL_WEIGHTS") << newweights;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
    
    return result;
    
  }//END coarsen
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry / Polyhedral geometry"
		    "# Takes a tropical variety on which a coarsest polyhedral structure exists"
		    "# and computes this structure."
		    "# @param WeightedComplex complex A tropical variety which has a unique "
		    "# coarsest polyhedral structre (e.g. it is locally connected in codimension "
		    "# one."
		    "# @return WeightedComplex The corresponding coarse complex",
		    &coarsen, "coarsen(WeightedComplex)");
  
  
}}