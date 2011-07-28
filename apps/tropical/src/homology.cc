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

This file contains functionality to compute homological properties of tropical
varieties.
*/

#include "polymake/client.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"

namespace polymake { namespace tropical {

    using namespace atint::donotlog;
    //using namespace atint::dolog;
    //using namespace atint::dotrace;
    
    perl::Object adjacencyComplex(perl::Object fan) {
	 //Extract values from fan
	 IncidenceMatrix<> maximalByCodim = fan.give("CODIM_1_IN_MAXIMAL_CONES");
	 IncidenceMatrix<> codimInMaximal = T(maximalByCodim);
	 IncidenceMatrix<> mcones = fan.give("MAXIMAL_CONES");
	 int noOfCones = mcones.rows();
	 
	 //Create adjacency list
	 Vector<Set<int> > facets;
	 for(int mc = 0; mc < noOfCones; mc++) {
	    //Go through all codim-1-faces of mc
	    Set<int> cdSet = codimInMaximal.row(mc);
	    for(Entire<Set<int> >::iterator cd = entire(cdSet); !cd.at_end(); cd++) {
	      //For all maximal cones adjacent to cd that are greater than mc, add
	      //an adjacency
	      Set<int> adjMax = maximalByCodim.row(*cd);
	      for(Entire<Set<int> >::iterator othermc = entire(adjMax); !othermc.at_end(); othermc++) {
		if(*othermc > mc) {
		    Set<int> s;
		    s += mc; s += *othermc;
		    facets |= s;
		}
	      }
	    }
	 }
	   
	 //Create and return the simplicial complex
	 perl::Object complex("topaz::SimplicialComplex");
	  complex.take("FACETS") << facets;
	 return complex;
	   
    }
    
    perl::Object equivalencyComplex(perl::Object fan) {
      //Extract values
      IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
      IncidenceMatrix<> codimOneCone = fan.give("CODIM_1_FACES");
      IncidenceMatrix<> maximalByCodim = fan.give("CODIM_1_IN_MAXIMAL_CONES");
      IncidenceMatrix<> codimInMaximal = T(maximalByCodim);
      int noOfCones = maximalCones.rows();
      Matrix<Rational> rays = fan.give("RAYS");
      Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
      Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
      //Vector<Integer> weights = fan.give("TROPICAL_WEIGHTS");
      
      dbgtrace << "Computing adjacency complex" << endl;
      
      //First we compute the adjacency complex
      perl::Object complex = adjacencyComplex(fan);
      
      //Now we compute the equivalence classes of maximal cones
      Vector<Set<int> > equivalenceClasses;
      Vector<bool> hasBeenAdded(noOfCones); //contains whether a cone has been added to an equivalence class
      
      dbgtrace << "Computing equivalence classes" << endl;
      
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
		dbgtrace << "Checking node " << node << endl;
		Set<int> cdset = codimInMaximal.row(node);
		for(Entire<Set<int> >::iterator cd = entire(cdset); !cd.at_end(); cd++) {
		    Set<int> otherMaximals = maximalByCodim.row(*cd) - node;
		    //We are only interested in codim-one-faces that have exactly two adjacent maximal cones
		    if(otherMaximals.size() == 1) {
		      int othermc = *(otherMaximals.begin());
		      if(!hasBeenAdded[othermc]) {
			//Now add the cone to the equivalence class of mc
			newset += othermc;
			queue.push_back(othermc);
			hasBeenAdded[othermc] = true;
		      }
		    }
		    
		    
		    //For each neighbour, if it has not been assigned a component yet, check if it should be added
// 		    for(Entire<Set<int> >::iterator othermc = entire(otherMaximals); !othermc.at_end(); othermc++) {
// 		      if(!hasBeenAdded[*othermc]) {
// // 			Matrix<Rational> vtau = rays.minor(codimOneCone.row(*cd),All) / linspace;
// // 			Vector<Rational> sum((latticeNormals[*cd])[node] + (latticeNormals[*cd])[*othermc]);
// // 			if(rank(vtau/sum) == rank(vtau)) {
// 			
// 			    //Now add the cone to the equivalence class of mc
// 			    newset += *othermc;
// 			    queue.push_back(*othermc);
// 			    hasBeenAdded[*othermc] = true;
// // 			}
// 		      }
// 		    }
		}		
	      }//End iterate queue
	if(newset.size() > 2) {equivalenceClasses |= newset;}
	}	
      } //End compute equivalence classes
      
      dbgtrace << "Equivalence classes are " << equivalenceClasses << endl;
      dbgtrace << "Removing edges from equivalence classes" << endl;
      
      //Now we remove all edges contained in the equivalence classes and add the classes as facets
      Vector<Set<int> > facets = complex.give("FACETS");
	dbgtrace << "Facets before are " << facets << endl;
      for(int ec = 0; ec < equivalenceClasses.dim(); ec++) {
	  Set<int> facetsToRemove;
	  for(int fc = 0; fc < facets.dim(); fc++) {
	    //Fancy way of saying facets[fc] contained in equivClasses[ec]
	    if((facets[fc]*equivalenceClasses[ec]).size() == facets[fc].size()) {
	      facetsToRemove += fc;
	    }
	  }
	  facets = facets.slice(~facetsToRemove);
      }
      facets |= equivalenceClasses;
      
      perl::Object result("topaz::SimplicialComplex");
	result.take("FACETS") << facets;
      return result;
    }
    
    UserFunction4perl("# @category Homology"
		      "# Computes the codimension one adjacency graph of a tropical variety and returns it as a "
		      "# simplicial complex. The nodes of the graph are the maximal cones, which are connected, iff "  
		      "# they share a codimension one face"
		      "# @param fan::PolyhedralFan fan A tropical variety"
		      "# @return topaz::SimplicialComplex The adjacency graph as a simplicial complex",
		      &adjacencyComplex, "adjacencyComplex(fan::PolyhedralFan)");
		      
    UserFunction4perl("# @category Homology"
		      "# Computes the canonical complex of a tropical variety. The one-skeleton is simply the "
		      "# adjacency graph of the variety. Higher dimensional faces are determined by equivalence"
		      "#  classes of cones. Two cones a,b are equivalent, if there exists a sequence of "
		      "# maximal cones a=c0,...,cr=b, such that ci,ci+1 intersect in a codimension one face of "
		      "# which they are the only neighbours"
		      "# @param fan::PolyhedralFan fan A tropical variety"
		      "# @return topaz::SimplicialComplex The canonical complex of the variety.",
		      &equivalencyComplex, "equivalenceComplex(fan::PolyhedralFan)");
    
}}