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
   Copyright (C) 2011 - 2015, Simon Hampe <simon.hampe@googlemail.com>

   Contains functions to compute the affine transform of a cycle 
   */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/PowerSet.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/thomog.h"
#include "polymake/matroid/check_axioms.h"
#include "polymake/Ring.h"
#include "polymake/Polynomial.h"

namespace polymake { namespace tropical {

   Set<int> minimal_nonempty(const int total, const Vector<Set<int> > &supports, const Set<int> &indices) {
      Set<int> result(sequence(0,total));
      for(Entire<Set<int> >::const_iterator iit = entire(indices); !iit.at_end(); iit++) {
         if( supports[*iit].size() != 0 && supports[*iit].size() < result.size()) {
            result = supports[*iit];
         }
      }
      return result;
   }

   perl::ListReturn chain_basis(int n, int r) {

      Set<Array<Set<int> > > collection;

      //Start by finding all chains starting with the emtpy set 
      perl::Object unnfan = CallPolymakeFunction("ufan_max",n);
      perl::Object skel = CallPolymakeFunction("skeleton_complex",unnfan,n-r-1,1);

      Matrix<Rational> vertices = unnfan.give("PROJECTIVE_VERTICES");
      Set<int> far_vertices = unnfan.give("FAR_VERTICES");
      int origin = *((sequence(0, vertices.rows()) - far_vertices).begin());
      vertices = vertices.minor(All,~scalar2set(0));
      Vector<Set<int> > supports (vertices.rows());
      Set<int> too_large;
      for(int r = 0; r < vertices.rows(); r++) {
         supports[r] = support(vertices.row(r));
         if(supports[r].size() > n-2) too_large += r;
      }
      IncidenceMatrix<> chains = skel.give("MAXIMAL_POLYTOPES");

      for(Entire<Rows<IncidenceMatrix<> > >::iterator cn = entire(rows(chains)); !cn.at_end(); cn++) {
         if( (*cn * too_large).size() == 0) {
            perl::Object cn_matroid = CallPolymakeFunction("chain_matroid", n, supports.slice(*cn));
            Array<Set<int> > cn_bases = cn_matroid.give("BASES");
            collection += cn_bases;
            //Also compute coloop-containing ones by replacing the first (empty set)
            Set<int> minne = minimal_nonempty(n, supports, *cn);
            Array<Set<int> > replacers = all_subsets(minne);
            for(Entire<Array<Set<int> > >::iterator rp = entire(replacers); !rp.at_end(); rp++) {
               if(!((*rp).size() == 0 || (*rp).size() == minne.size())) {
                  perl::Object cl_matroid = CallPolymakeFunction("chain_matroid", n, supports.slice(*cn - origin) | Set<int>(*rp)); 
                  Array<Set<int> > cl_bases = cl_matroid.give("BASES");
                  collection += cl_bases;
               }
            }
         }
      }

      perl::ListReturn result;
      for(Entire<Set<Array<Set<int> > > >::iterator cl = entire(collection); !cl.at_end(); cl++) {
         perl::Object mat("matroid::Matroid");
         mat.take("N_ELEMENTS") << n;
         mat.take("BASES") << *cl;
         result << mat;
      }

      return result;

   }

 /*  Vector<Set<int> > find_chain( int n, int r, perl::Object matroid) {
      Set<Array<Set<int> > > collection;

      Array<Set<int> > mbases = matroid.give("BASES");

      //Start by finding all chains starting with the emtpy set 
      perl::Object unnfan = CallPolymakeFunction("ufan_max",n);
      perl::Object skel = CallPolymakeFunction("skeleton_complex",unnfan,n-r-1,1);

      Matrix<Rational> vertices = unnfan.give("PROJECTIVE_VERTICES");
      vertices = vertices.minor(All,~scalar2set(0));
      Vector<Set<int> > supports (vertices.rows());
      Set<int> too_large;
      for(int r = 0; r < vertices.rows(); r++) {
         supports[r] = support(vertices.row(r));
         if(supports[r].size() > n-2) too_large += r;
      }
      IncidenceMatrix<> chains = skel.give("MAXIMAL_POLYTOPES");

      for(Entire<Rows<IncidenceMatrix<> > >::iterator cn = entire(rows(chains)); !cn.at_end(); cn++) {
         if( (*cn * too_large).size() == 0) {
            perl::Object cn_matroid = CallPolymakeFunction("chain_matroid", n, supports.slice(*cn));
            Array<Set<int> > cn_bases = cn_matroid.give("BASES");
            if(cn_bases == mbases) return supports.slice(*cn);
         }
      }
      return Vector<Set<int> >();
   }
*/
   Integer eulerian(int n, int k) {
      if(n == k && n == 0) return 1;
      if(n == 0) return 0;
      if(k == 0) return 1;
      return Integer(n-k)*eulerian(n-1,k-1) + Integer(k+1)*eulerian(n-1,k);
   }

  /* Integer nested_number(int n, int r) {
      if(n == r && n == 0) return 1;
      if(n == 0 || r == 0) return 0;
      if(n < r) return 0;
      if(r == 1) return 1;
      Integer result(1);
      for(int k = 1; k <= r-1; k++) {
         for(int s = k+1; s <= k+n-r; s++) {
            result += Integer::binom(n,s) * nested_number(n-s,r-k);
         }
      }
      return result;
   }*/

   UserFunction4perl("", &chain_basis, "chain_basis($,$)");
//   UserFunction4perl("", &find_chain, "find_chain($,$,matroid::Matroid)");
   UserFunction4perl("", &eulerian, "eulerian($,$)");
 //  UserFunction4perl("", &nested_number,"nested_number($,$)");

}}
