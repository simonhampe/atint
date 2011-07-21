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

#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"

namespace polymake { namespace fan {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object, bool, pm::Rational, pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturnVoid( arg0, arg1, arg2.get< perl::TryCanned< const Rational > >(), arg3.get< perl::TryCanned< const Matrix< Rational > > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object, bool, pm::Rational, pm::Matrix<pm::Rational>) );

   OperatorInstance4perl(assign, pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational>&>, pm::Series<int, true>, void>, perl::Canned< const pm::incidence_line<pm::AVL::tree<pm::sparse2d::traits<pm::sparse2d::traits_base<pm::nothing, true, false, (pm::sparse2d::restriction_kind)0>, false, (pm::sparse2d::restriction_kind)0> > const&> >);
   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object, bool, pm::Rational, pm::Matrix<pm::Rational>, bool) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]), arg4(stack[4]);
      IndirectWrapperReturnVoid( arg0, arg1, arg2.get< perl::TryCanned< const Rational > >(), arg3.get< perl::TryCanned< const Matrix< Rational > > >(), arg4 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object, bool, pm::Rational, pm::Matrix<pm::Rational>, bool) );

   FunctionWrapper4perl( pm::perl::ListReturn (perl::Object, bool, bool, pm::Rational, pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]), arg4(stack[4]);
      IndirectWrapperReturnVoid( arg0, arg1, arg2, arg3.get< perl::TryCanned< const Rational > >(), arg4.get< perl::TryCanned< const Matrix< Rational > > >() );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (perl::Object, bool, bool, pm::Rational, pm::Matrix<pm::Rational>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
