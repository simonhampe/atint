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
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"

namespace polymake { namespace atint {

    using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
    
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
//       Matrix<Rational> rays = fan.give("RAYS");
//       Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
//       Map<int, Map<int, Vector<Integer> > > latticeNormals = fan.give("LATTICE_NORMALS");
      //Vector<Integer> weights = fan.give("TROPICAL_WEIGHTS");
      
      dbgtrace << "Computing adjacency complex" << endl;
      
      //First we compute the adjacency complex
      perl::Object complex = adjacencyComplex(fan);
      
      //Now we compute the equivalence classes of maximal cones
      Vector<Set<int> > equivalenceClasses;
      Vector<bool> hasBeenAdded(noOfCones); //contains whether a cone has been added to an equivalence class
      Map<int, int> conesInClasses; //Maps cone indices to the index of their class in equivalenceClasses
      
      dbgtrace << "Computing equivalence classes" << endl;
      
      for(int mc = 0; mc < noOfCones; mc++) {
	if(!hasBeenAdded[mc]) {
	      Set<int> newset; newset = newset + mc;
	      conesInClasses[mc] = equivalenceClasses.dim();
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
			conesInClasses[othermc] = equivalenceClasses.dim();
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
	//if(newset.size() > 2) {equivalenceClasses |= newset;}
	equivalenceClasses |= newset;
	}	
      } //End compute equivalence classes
      
      dbglog << "Equivalence classes are " << equivalenceClasses << endl;
      dbglog << "Cones in classes reads " << conesInClasses << endl;
      //dbgtrace << "Removing edges from equivalence classes" << endl;
      
      //Now we take each equivalence class, iterate over all cones in this class and add as connected
      // vertices all classes containing a cone connected to one of the cones in this class in the adjacencyComplex
      Vector<Set<int> > facets;
      Vector<Set<int> > adjacency = complex.give("FACETS");
      IncidenceMatrix<> vertexInEdge = T(IncidenceMatrix<>(adjacency));
      dbgtrace << "vertex-in-edge matrix reads " << vertexInEdge << endl;
      for(int cls = 0; cls < equivalenceClasses.dim(); cls++) {
	dbgtrace << "Finding connected classes for class " << cls << endl;
	Set<int> c = equivalenceClasses[cls];
	Set<int> connectedClasses;
	dbgtrace << "Cones in that class are " << c << endl;
	//Iterate over all cones in the class
	for(Entire<Set<int> >::iterator cone = entire(c); !cone.at_end(); cone++) {
	    //Extract all edges containing that cone
	    Set<int> edgeSet = vertexInEdge.row(*cone);
	    dbgtrace << "Cone " << *cone << " contained in edges " << edgeSet << endl;
	    //Now add all neighbour cones
	    for(Entire<Set<int> >::iterator nb = entire(edgeSet); !nb.at_end(); nb++) {
	      connectedClasses += conesInClasses[*((adjacency[*nb] - *cone).begin())];
	    }
	}
	dbgtrace << "Connected classes are " << connectedClasses << endl;
	//Now that we know all connected classes, add the corresponding edges (only if the index is > )
	for(Entire<Set<int> >::iterator cc = entire(connectedClasses); !cc.at_end(); cc++) {
	    if(*cc > cls) {
	      Set<int> edge; edge += cls; edge += *cc;
	      facets |= edge;
	    }
	}
	//If there are no connected classes, create a vertex
	if(connectedClasses.size() == 0) {
	    Set<int> vertex; vertex += cls;
	    facets |= vertex;
	}
      }
      
      dbgtrace << "Done. Facets are " << facets << endl;
      
      //Now we remove all edges contained in the equivalence classes and add the classes as facets
//       Vector<Set<int> > facets = complex.give("FACETS");
// 	dbgtrace << "Facets before are " << facets << endl;
//       for(int ec = 0; ec < equivalenceClasses.dim(); ec++) {
// 	  Set<int> facetsToRemove;
// 	  for(int fc = 0; fc < facets.dim(); fc++) {
// 	    //Fancy way of saying facets[fc] contained in equivClasses[ec]
// 	    if((facets[fc]*equivalenceClasses[ec]).size() == facets[fc].size()) {
// 	      facetsToRemove += fc;
// 	    }
// 	  }
// 	  facets = facets.slice(~facetsToRemove);
//       }
//       facets |= equivalenceClasses;
      
      perl::Object result("topaz::SimplicialComplex");
	result.take("FACETS") << facets;
      return result;
    }
    
    UserFunction4perl("# @category Homology"
		      "# Computes the codimension one adjacency graph of a tropical variety and returns it as a "
		      "# simplicial complex. The nodes of the graph are the maximal cones, which are connected, iff "  
		      "# they share a codimension one face"
		      "# @param WeightedComplex fan A tropical variety"
		      "# @return topaz::SimplicialComplex The adjacency graph as a simplicial complex",
		      &adjacencyComplex, "adjacencyComplex(WeightedComplex)");
		      
    UserFunction4perl("# @category Homology"
		      "# Computes the canonical complex of a tropical variety. The one-skeleton is determined "
		      "# in the following way: Two cones a,b are equivalent, if there exists a sequence of "
		      "# maximal cones a=c0,...,cr=b, such that ci,ci+1 intersect in a codimension one face of "
		      "# which they are the only neighbours"
		      "# The vertices are the equivalence classes and two of them are connected if any cones in "
		      "# them are connected."
		      "# @param WeightedComplex fan A tropical variety"
		      "# @return topaz::SimplicialComplex The canonical complex of the variety.",
		      &equivalencyComplex, "equivalenceComplex(WeightedComplex)");
    
}}
