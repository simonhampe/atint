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

namespace polymake { namespace fan {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( void (perl::Object, pm::Array<pm::Matrix<pm::Rational>, void>, pm::Array<pm::Matrix<pm::Rational>, void>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Array< Matrix< Rational > > > >(), arg2.get< perl::TryCanned< const Array< Matrix< Rational > > > >() );
   }
   FunctionWrapperInstance4perl( void (perl::Object, pm::Array<pm::Matrix<pm::Rational>, void>, pm::Array<pm::Matrix<pm::Rational>, void>) );

   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object const&, pm::Rational const&, bool) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Rational > >(), arg2 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object const&, pm::Rational const&, bool) );

   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object const&, pm::Rational const&, bool, bool) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Rational > >(), arg2, arg3 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object const&, pm::Rational const&, bool, bool) );

   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object const&, pm::Rational const&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Rational > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object const&, pm::Rational const&) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
