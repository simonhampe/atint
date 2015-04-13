/* Copyright (c) 1997-2015
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Berlin, Germany)
   http://www.polymake.org

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

namespace polymake { namespace tropical { namespace {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( pm::Matrix<pm::Rational> (pm::Matrix<pm::Rational>, pm::Rational, bool) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix< Rational > > >(), arg1.get< perl::TryCanned< const Rational > >(), arg2 );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Rational> (pm::Matrix<pm::Rational>, pm::Rational, bool) );

   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object, pm::Vector<pm::Integer>, pm::Matrix<pm::Rational>, pm::Array<std::string, void>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturnVoid( arg0, arg1.get< perl::TryCanned< const Vector< Integer > > >(), arg2.get< perl::TryCanned< const Matrix< Rational > > >(), arg3.get< perl::TryCanned< const Array< std::string > > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object, pm::Vector<pm::Integer>, pm::Matrix<pm::Rational>, pm::Array<std::string, void>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
