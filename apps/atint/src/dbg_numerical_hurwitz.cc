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
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
//   using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Vector<int> findShelling(perl::Object cycle) {
    
    //Extract values
    IncidenceMatrix<> cones = cycle.give("MAXIMAL_CONES");
    IncidenceMatrix<> codim = cycle.give("CODIM_1_FACES");
    IncidenceMatrix<> maxInCodim = cycle.give("CODIM_1_IN_MAXIMAL_CONES");
      maxInCodim = T(maxInCodim);
    
    //Initialize backtrack variables
    Vector<int> currentSequence;
    Vector<Set<int> > triedCones(cones.rows());
    
    while(currentSequence.dim() < cones.rows()) {
      dbgtrace << "Current sequence: " << currentSequence << endl;
      dbgtrace << "Tried cones: " << triedCones[currentSequence.dim()] << endl;
      //Find next untried cone
      int tryIndex = 0;
      while( (Set<int>(currentSequence).contains(tryIndex) || triedCones[currentSequence.dim()].contains(tryIndex)) && tryIndex < triedCones.dim()) {
	tryIndex++;
      }
      
      //If we find none, go back a step
      if(tryIndex >= triedCones.dim()) {
	if(currentSequence.dim() > 0) {
	  triedCones[currentSequence.dim()] = Set<int>();
	  currentSequence = currentSequence.slice(sequence(0,currentSequence.dim()-1));
	}
	continue;
      }
      
      dbgtrace << "Trying to add " << tryIndex << endl;
      
      //Otherwise, check if this works
      triedCones[currentSequence.dim()] += tryIndex;
      
      //Compute union of all previous cones
      Set<int> prevUnion;
      for(int i = 0; i < currentSequence.dim(); i++) {
	prevUnion += cones.row(currentSequence[i]);
      }
      prevUnion = prevUnion * cones.row(tryIndex);
      
      dbgtrace << "Union of previous cones: " << prevUnion << endl;
      
      //Then see if it is a union of codim 1 faces of tryIndex - cone
      Set<int> tryCodims = maxInCodim.row(tryIndex);
      Set<int> codimUnion;
      for(Entire<Set<int> >::iterator cd = entire(tryCodims); !cd.at_end(); cd++) {
	if( (codim.row(*cd) * prevUnion).size() == codim.row(*cd).size()) {
	    codimUnion += codim.row(*cd);
	}
      }
      
      dbgtrace << "Union of codim 1 faces: " << codimUnion << endl;
      
      if(codimUnion.size() == prevUnion.size()) {
	currentSequence |= tryIndex;
      }
    }
    
    return currentSequence;
    
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&findShelling,"fshell(WeightedComplex)");
  
}}