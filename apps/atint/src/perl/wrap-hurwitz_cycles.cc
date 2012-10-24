/* Copyright (c) 1997-2010
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
   http://www.polymake.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*/

#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"

namespace polymake { namespace atint {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( perl::Object (int, pm::Vector<int>, pm::Vector<pm::Rational>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0, arg1.get< perl::TryCanned< const Vector< int > > >(), arg2.get< perl::TryCanned< const Vector< Rational > > >() );
   }
   FunctionWrapperInstance4perl( perl::Object (int, pm::Vector<int>, pm::Vector<pm::Rational>) );

   OperatorInstance4perl(assign, pm::incidence_line<pm::AVL::tree<pm::sparse2d::traits<pm::sparse2d::traits_base<pm::nothing, true, false, (pm::sparse2d::restriction_kind)2>, false, (pm::sparse2d::restriction_kind)2> > >, perl::Canned< const pm::Series<int, true> >);
   FunctionWrapper4perl( perl::Object (pm::Vector<int>) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Vector< int > > >() );
   }
   FunctionWrapperInstance4perl( perl::Object (pm::Vector<int>) );

   FunctionWrapper4perl( pm::Integer (pm::Vector<int>) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Vector< int > > >() );
   }
   FunctionWrapperInstance4perl( pm::Integer (pm::Vector<int>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
