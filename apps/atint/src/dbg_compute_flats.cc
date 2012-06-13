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
#include "polymake/PowerSet.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  void computeFlats(Matrix<Rational> m) {
    int r = rank(m);
    
    Vector<Vector<Set<int> > > flats(r);
    
    Set<int> all_cols = sequence(0,m.cols());
    Array<Set<int> > ssets = all_subsets(all_cols);    
    for(int s = 0; s < ssets.size(); s++) {
      int x = rank(m.minor(All,ssets[s]));
      if(x == 0) continue;
      //Check if s is a flat
      Set<int> complement = all_cols - ssets[s];
      bool found_bad_element = false;
      for(Entire<Set<int> >::iterator c = entire(complement); !c.at_end(); c++) {
	if(rank(m.minor(All,ssets[s] + *c)) == x) {
	  found_bad_element = true; break;
	}
      }
      if(!found_bad_element) {
	flats[x-1] |= ssets[s];
      }
    }
    
    pm::cout << "Done: " << endl;
    
    //Output
    for(int i = 0; i < flats.size(); i++) {
      pm::cout << "Rank " << (i+1) << ": " << flats[i] << endl;
    }
    
  }

  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  Function4perl(&computeFlats, "computeFlats(Matrix<Rational>)");
  
}}