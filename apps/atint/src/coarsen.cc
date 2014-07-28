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
#include "polymake/linalg.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  using polymake::polytope::cdd_interface::solver;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   * @brief Takes a polyhedron in its V-description and tests whether it fulfills a given inequality
   * @param Matrix<Rational> rays The rays/vertices of the polyhedron
   * @param Matrix<Rational> linspace The lineality space generators of the polyhedron
   * @param Vector<Rational> facet The inequality, in standard polymake format. Note that all three should be given in the same format, i.e. in homog. or non-homog. coordinates.
   * @return True, if the polyhedron fulfills the inequality; false otherwise.
   */
  bool coneInHalfspace(const Matrix<Rational> &rays, const Matrix<Rational> &linspace,
			const Vector<Rational> &facet) {
      Matrix<Rational> allgenerators = rays / linspace;
      Vector<Rational> product = allgenerators * facet;
      for(int i = 0; i < product.dim(); i++) {
	if(product[i] < 0) return false;
      }
      return true;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object coarsen(perl::Object complex, bool testFan = false) {
    //Extract values
    Matrix<Rational> rays = complex.give("RAYS");
    Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
    int cmplx_dim = complex.give("CMPLX_DIM");
    int cmplx_ambient_dim = complex.give("CMPLX_AMBIENT_DIM");
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    bool weights_exist = complex.exists("TROPICAL_WEIGHTS");
    Vector<Integer> weights;
    if(weights_exist) {
	complex.give("TROPICAL_WEIGHTS") >> weights;
    }    
    IncidenceMatrix<> maximalCones = complex.give("MAXIMAL_CONES");
    IncidenceMatrix<> codimOneCones = complex.give("CODIM_1_FACES");
    IncidenceMatrix<> codimInMaximal = complex.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximalOverCodim = T(codimInMaximal);
    int noOfCones = maximalCones.rows();
    
    //For fan testing, we need the equations of all non-twovalent codimension one faces
    Vector<Matrix<Rational> > codimOneEquations(codimOneCones.rows());
    Set<int> highervalentCodim;
    for(int cd = 0; cd < codimOneCones.rows(); cd++) { 
	if(codimInMaximal.row(cd).size() > 2) {
	  codimOneEquations[cd] = null_space( rays.minor(codimOneCones.row(cd),All) / linspace);
	  highervalentCodim += cd;
	}//END only for >2-valent
    }//END compute codim one cone equations.
    
    //If the fan has no rays, it has only lineality space and we return the complex itself
    if(rays.rows() == 0)  return complex;
    
    //If the fan has no codim 1 faces, it must be 0- dimensional or a lineality space and we return
    // the complex itself
    if(codimOneCones.rows() == 0) return complex;
    
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
    
    /*Test equivalence classes for convexit:
    An equivalence class S = (s1,...,sr) is convex, iff it fulfills the following criterion: Let T be the set of all outer (i.e. not two-valent) codim one facets and write f_t >= a_t for the corresponding facet inequality. Then for each pair t,t' in T the polyhedron t must fulfill the inequality f_t' >= a_t'.
    */
    if(testFan) {
	//Check each equivalence class
	for(int cl = 0; cl < equivalenceClasses.dim(); cl++) {
	  //dbgtrace << "Checking convexity of class " << cl << endl;
	  //First, we compute the equation of the subdivision class
	  Matrix<Rational> classEquation =
	    null_space(rays.minor(maximalCones.row(*(equivalenceClasses[cl].begin())),All) / linspace);
	  //dbgtrace << "Class equations are " << classEquation << endl;
	  //Compute all outer codim one faces and their facet inequality;
	  Set<int> outerFacets;
	  Map<int,Vector<Rational> > facetInequality;
	  for(Entire<Set<int> >::iterator mc = entire(equivalenceClasses[cl]); !mc.at_end(); mc++) {
	    //dbgtrace << "Checking maximal cone " << (*mc) << endl;
	    Set<int> outerCodimInMc = maximalOverCodim.row(*mc) * highervalentCodim;
	    //dbgtrace << "Outer facets are " << outerCodimInMc << endl;
	    outerFacets += outerCodimInMc;
	    for(Entire<Set<int> >::iterator fct = entire(outerCodimInMc); !fct.at_end(); fct++) {
		//Find the one equation of fct that is not an equation of the class
		//dbgtrace << "Computing facet inequality of cone " << (*mc) << " and facet " << (*fct) << endl;
		//dbgtrace << "Facet equations are " << codimOneEquations[*fct] << endl;
		for(int r = 0; r < codimOneEquations[*fct].rows(); r++) {
		  if(rank(classEquation / codimOneEquations[*fct].row(r)) > cmplx_ambient_dim - cmplx_dim) {
		    //Find the right sign and save the inequality
		    //dbgtrace << "Row " << r << " is the additional equation" << endl;
		    int additionalRay = *( (maximalCones.row(*mc) - codimOneCones.row(*fct)).begin());
		    facetInequality[*fct] = codimOneEquations[*fct].row(r);
		    if(facetInequality[*fct] * rays.row(additionalRay) < 0) {
		      facetInequality[*fct] = - facetInequality[*fct];
		    }		    
// 		    dbgtrace << "Right inequality is " << facetInequality[*fct] << endl;
		    break;
		  }
		}//END iterate equations of fct
	    }//END iterate all outer facets of the cone
	  }//END iterate all maximal cones in the class
	  //dbgtrace << "Checking convexity criterion" << endl;
	  //Now go through all pairs of outer facets
	  for(Entire<Set<int> >::iterator out1 = entire(outerFacets); !out1.at_end(); out1++) {
	      Set<int> remainingOuter = outerFacets - *out1;
	      for(Entire<Set<int> >::iterator out2 = entire(remainingOuter); !out2.at_end(); out2++) {
		  //Check if out2 fulfills the facet inequality of out1
		  if(!coneInHalfspace(rays.minor(codimOneCones.row(*out2),All), linspace, facetInequality[*out1])) {
		    throw std::runtime_error("The equivalence classes are not convex. There is no coarsest structure.");
		    
		  }		  
	      }//END iterate second facet.
	  }//END iterate first facet
	}//END iterate classes
    }//END test convexity of equivalence classes
    
    //Now compute the new cones as unions of the cones in each equivalence class
    Matrix<Rational> newlin;
      int newlindim = 0;
      bool newlin_computed = false;
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
	if(testFan) newlindim = rank(newlin);
	newlin_computed = true;
      }
      //If we need to test fan-ness, we need to check if different classes have different lin. spaces
      if(testFan && newlin_computed) {
	Matrix<Rational> otherlin = (union_rays/linspace).minor(union_cone.second,All);
	if(newlindim != rank(otherlin)) {
	    throw std::runtime_error("Different classes have different lineality spaces. There is no coarsest structure.");
	}
	else {
	    if(rank(newlin / otherlin) > rank(newlin)) {
	      throw std::runtime_error("Different classes have different lineality spaces. There is no coarsest structure.");
	    }
	}
      }//END check equality of lineality spaces.
      
      //Convert indices of rays in union_rays to indices in rays
      Set<int> ray_set(union_ray_list.slice(union_cone.first));
      
      newcones |= ray_set;
      used_rays += ray_set;
      if(weights_exist) newweights |= weights[*(conesInClass.begin())];      
    }
    
    //
    
    //Some rays might become equal (modulo lineality space) when coarsening, so we have to clean up
    Map<int,int> ray_index_conversion;
    int next_index = 0;
    Matrix<Rational> final_rays(0,complete_matrix.cols());
    int linrank = rank(newlin);
    
    for(Entire<Set<int> >::iterator r1 = entire(used_rays); !r1.at_end(); r1++) {
      if(!keys(ray_index_conversion).contains(*r1)) {
	ray_index_conversion[*r1] = next_index;
	final_rays /= complete_matrix.row(*r1);
	
	for(Entire<Set<int> >::iterator r2 = entire(used_rays); !r2.at_end(); r2++) {
	  if(!keys(ray_index_conversion).contains(*r2)) {
	      //Check if both rays are equal mod lineality space
	      if( rank(newlin / (complete_matrix.row(*r1) - complete_matrix.row(*r2))) == linrank) {
		ray_index_conversion[*r2] = next_index;
	      }
	  }//END if second not mapped yet
	}//END iterate second ray index
	next_index++;
      }//END if first not mapped yet
    }//END iterate first ray index
    
    //Now convert the cones
    for(int nc = 0; nc < newcones.dim(); nc++) {
      newcones[nc] = attach_operation(newcones[nc],pm::operations::associative_access<Map<int,int>,int>(&ray_index_conversion));
    }
    
    
    //If necessary, test for fan-ness
    if(testFan) {
      //This function throws an error if the result is not a fan.
      try {
	CallPolymakeFunction("fan::check_fan",final_rays, newcones);
      }
      catch(...) {
	throw std::runtime_error("The equivalence classes do not form a fan. There is no coarsest structure.");
      }
    }
    
    //Produce final result
    perl::Object result("WeightedComplex");
      result.take("RAYS") << final_rays;
      result.take("MAXIMAL_CONES") << newcones;
      result.take("LINEALITY_SPACE") << newlin;
      if(weights_exist) {
	result.take("TROPICAL_WEIGHTS") << newweights;
      }
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
    
    return result;
    
  }//END coarsen
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Basic polyhedral operations"
		    "# Takes a tropical variety on which a coarsest polyhedral structure exists"
		    "# and computes this structure."
		    "# @param WeightedComplex complex A tropical variety which has a unique "
		    "# coarsest polyhedral structre "
		    "# @param Bool testFan Whether the algorithm should test if the result"
		    "# is actually a polyhedral complex. This is false by default and should "
		    "# only be set to true, if the user is not actually certain, whether the "
		    "# complex in question has a coarsest structure. If the result is not a"
		    "# fan, the algorithm throws an exception."
		    "# @return WeightedComplex The corresponding coarse complex. Throws an "
		    "# exception if testFan = True and the result is not a fan",
		    &coarsen, "coarsen(WeightedComplex; $=0)");
  
}}