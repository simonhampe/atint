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
   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object, bool, bool, pm::Rational, bool, pm::Matrix<pm::Rational>, pm::Array<std::string, void>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]), arg4(stack[4]), arg5(stack[5]), arg6(stack[6]);
      IndirectWrapperReturnVoid( arg0, arg1, arg2, arg3.get< perl::TryCanned< const Rational > >(), arg4, arg5.get< perl::TryCanned< const Matrix< Rational > > >(), arg6.get< perl::TryCanned< const Array< std::string > > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object, bool, bool, pm::Rational, bool, pm::Matrix<pm::Rational>, pm::Array<std::string, void>) );

   FunctionWrapper4perl( pm::Matrix<pm::Rational> (pm::Matrix<pm::Rational>, pm::Rational) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix< Rational > > >(), arg1.get< perl::TryCanned< const Rational > >() );
   }
   FunctionWrapperInstance4perl( pm::Matrix<pm::Rational> (pm::Matrix<pm::Rational>, pm::Rational) );

   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object, bool, bool, pm::Rational, pm::Matrix<pm::Rational>, pm::Array<std::string, void>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]), arg4(stack[4]), arg5(stack[5]);
      IndirectWrapperReturnVoid( arg0, arg1, arg2, arg3.get< perl::TryCanned< const Rational > >(), arg4.get< perl::TryCanned< const Matrix< Rational > > >(), arg5.get< perl::TryCanned< const Array< std::string > > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object, bool, bool, pm::Rational, pm::Matrix<pm::Rational>, pm::Array<std::string, void>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
