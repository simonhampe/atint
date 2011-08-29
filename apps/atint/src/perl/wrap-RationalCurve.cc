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

namespace polymake { namespace atint {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( int (int) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0 );
   }
   FunctionWrapperInstance4perl( int (int) );

   FunctionWrapper4perl( void () ) {
      IndirectWrapperReturnVoid(  );
   }
   FunctionWrapperInstance4perl( void () );

   FunctionWrapper4perl( pm::Vector<pm::Rational> (pm::IncidenceMatrix<pm::NonSymmetric>, pm::Vector<pm::Rational>, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg1.get< perl::TryCanned< const Vector< Rational > > >(), arg2 );
   }
   FunctionWrapperInstance4perl( pm::Vector<pm::Rational> (pm::IncidenceMatrix<pm::NonSymmetric>, pm::Vector<pm::Rational>, int) );

   FunctionWrapper4perl( pm::Vector<pm::Rational> (perl::Object, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1 );
   }
   FunctionWrapperInstance4perl( pm::Vector<pm::Rational> (perl::Object, int) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
