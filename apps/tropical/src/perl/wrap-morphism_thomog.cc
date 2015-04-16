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
   FunctionWrapper4perl( std::pair<pm::Matrix<pm::Rational>, pm::Vector<pm::Rational> > (pm::Matrix<pm::Rational> const&, pm::Vector<pm::Rational> const&, int, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix< Rational > > >(), arg1.get< perl::TryCanned< const Vector< Rational > > >(), arg2, arg3 );
   }
   FunctionWrapperInstance4perl( std::pair<pm::Matrix<pm::Rational>, pm::Vector<pm::Rational> > (pm::Matrix<pm::Rational> const&, pm::Vector<pm::Rational> const&, int, int) );

   FunctionWrapper4perl( std::pair<pm::Matrix<pm::Rational>, pm::Vector<pm::Rational> > (pm::Matrix<pm::Rational> const&, pm::Vector<pm::Rational> const&, int, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix< Rational > > >(), arg1.get< perl::TryCanned< const Vector< Rational > > >(), arg2, arg3 );
   }
   FunctionWrapperInstance4perl( std::pair<pm::Matrix<pm::Rational>, pm::Vector<pm::Rational> > (pm::Matrix<pm::Rational> const&, pm::Vector<pm::Rational> const&, int, int) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
