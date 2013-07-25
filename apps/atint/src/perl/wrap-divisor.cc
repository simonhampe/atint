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
   FunctionWrapper4perl( perl::Object (perl::Object, perl::Object, bool) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0, arg1, arg2 );
   }
   FunctionWrapperInstance4perl( perl::Object (perl::Object, perl::Object, bool) );

   FunctionWrapper4perl( perl::Object (perl::Object, perl::Object, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0, arg1, arg2 );
   }
   FunctionWrapperInstance4perl( perl::Object (perl::Object, perl::Object, int) );

   FunctionWrapper4perl( perl::Object (perl::Object, pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1.get< perl::TryCanned< const Matrix< Rational > > >() );
   }
   FunctionWrapperInstance4perl( perl::Object (perl::Object, pm::Matrix<pm::Rational>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
