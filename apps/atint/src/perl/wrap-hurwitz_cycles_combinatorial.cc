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
   FunctionWrapper4perl( pm::perl::ListReturn (int, pm::Vector<int>, pm::Vector<pm::Rational>, std::vector<perl::Object, std::allocator<perl::Object> >) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Vector< int > > >(), arg2.get< perl::TryCanned< const Vector< Rational > > >(), arg3 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (int, pm::Vector<int>, pm::Vector<pm::Rational>, std::vector<perl::Object, std::allocator<perl::Object> >) );

   FunctionWrapper4perl( pm::perl::ListReturn (int, pm::Vector<int>, pm::Vector<pm::Rational>, perl::Object) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Vector< int > > >(), arg2.get< perl::TryCanned< const Vector< Rational > > >(), arg3 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (int, pm::Vector<int>, pm::Vector<pm::Rational>, perl::Object) );

   FunctionWrapper4perl( pm::perl::ListReturn (int, pm::Vector<int>, perl::Object) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Vector< int > > >(), arg2 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (int, pm::Vector<int>, perl::Object) );

   FunctionWrapper4perl( pm::perl::ListReturn (int, pm::Vector<int>, pm::Vector<pm::Rational>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Vector< int > > >(), arg2.get< perl::TryCanned< const Vector< Rational > > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (int, pm::Vector<int>, pm::Vector<pm::Rational>) );

   FunctionWrapper4perl( perl::Object (perl::Object, pm::Vector<int>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1.get< perl::TryCanned< const Vector< int > > >() );
   }
   FunctionWrapperInstance4perl( perl::Object (perl::Object, pm::Vector<int>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
