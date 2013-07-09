/* Copyright (c) 1997-2010
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*/

namespace polymake { namespace atint { namespace {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( perl::Object (perl::Object, perl::Object, pm::Vector<int>, int, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]), arg4(stack[4]);
      IndirectWrapperReturn( arg0, arg1, arg2.get< perl::TryCanned< const Vector< int > > >(), arg3, arg4 );
   }
   FunctionWrapperInstance4perl( perl::Object (perl::Object, perl::Object, pm::Vector<int>, int, int) );

   FunctionWrapper4perl( pm::Vector<pm::Integer> (perl::Object, perl::Object) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1 );
   }
   FunctionWrapperInstance4perl( pm::Vector<pm::Integer> (perl::Object, perl::Object) );

   FunctionWrapper4perl( pm::Matrix<pm::Rational> (perl::Object, perl::Object, pm::Vector<int>, int, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]), arg4(stack[4]);
      IndirectWrapperReturn( arg0, arg1, arg2.get< perl::TryCanned< const Vector< int > > >(), arg3, arg4 );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Rational> (perl::Object, perl::Object, pm::Vector<int>, int, int) );

   FunctionWrapper4perl( pm::Matrix<pm::Rational> (perl::Object, perl::Object) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1 );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Rational> (perl::Object, perl::Object) );

   FunctionWrapper4perl( pm::Map<int, int, pm::operations::cmp> (perl::Object, perl::Object) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1 );
   }
   FunctionWrapperInstance4perl( pm::Map<int, int, pm::operations::cmp> (perl::Object, perl::Object) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
