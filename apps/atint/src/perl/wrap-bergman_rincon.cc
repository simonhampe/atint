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
   FunctionWrapper4perl( pm::Vector<pm::Set<int, pm::operations::cmp> > (pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix<Rational> > >() );
   }
   FunctionWrapperInstance4perl( pm::Vector<pm::Set<int, pm::operations::cmp> > (pm::Matrix<pm::Rational>) );

   FunctionWrapper4perl( pm::IncidenceMatrix<pm::NonSymmetric> (pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix<Rational> > >() );
   }
   FunctionWrapperInstance4perl( pm::IncidenceMatrix<pm::NonSymmetric> (pm::Matrix<pm::Rational>) );

   FunctionWrapper4perl( pm::Set<int, pm::operations::cmp> (pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix<Rational> > >() );
   }
   FunctionWrapperInstance4perl( pm::Set<int, pm::operations::cmp> (pm::Matrix<pm::Rational>) );

   FunctionWrapper4perl( pm::Set<int, pm::operations::cmp> (pm::IncidenceMatrix<pm::NonSymmetric>, int, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg1, arg2 );
   }
   FunctionWrapperInstance4perl( pm::Set<int, pm::operations::cmp> (pm::IncidenceMatrix<pm::NonSymmetric>, int, int) );

   FunctionWrapper4perl( void (pm::IncidenceMatrix<pm::NonSymmetric>, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturnVoid( arg0.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg1 );
   }
   FunctionWrapperInstance4perl( void (pm::IncidenceMatrix<pm::NonSymmetric>, int) );

   FunctionWrapper4perl( void (pm::IncidenceMatrix<pm::NonSymmetric>, int, pm::Matrix<pm::Rational>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturnVoid( arg0.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg1, arg2.get< perl::TryCanned< const Matrix<Rational> > >() );
   }
   FunctionWrapperInstance4perl( void (pm::IncidenceMatrix<pm::NonSymmetric>, int, pm::Matrix<pm::Rational>) );

   FunctionWrapper4perl( pm::Set<int, pm::operations::cmp> (pm::IncidenceMatrix<pm::NonSymmetric> const&, int, int, pm::Matrix<pm::Rational> const&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg1, arg2, arg3.get< perl::TryCanned< const Matrix<Rational> > >() );
   }
   FunctionWrapperInstance4perl( pm::Set<int, pm::operations::cmp> (pm::IncidenceMatrix<pm::NonSymmetric> const&, int, int, pm::Matrix<pm::Rational> const&) );

   FunctionWrapper4perl( perl::Object (int, pm::IncidenceMatrix<pm::NonSymmetric> const&, bool, pm::Matrix<pm::Rational> const&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturn( arg0, arg1.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg2, arg3.get< perl::TryCanned< const Matrix<Rational> > >() );
   }
   FunctionWrapperInstance4perl( perl::Object (int, pm::IncidenceMatrix<pm::NonSymmetric> const&, bool, pm::Matrix<pm::Rational> const&) );

   FunctionWrapper4perl( perl::Object (pm::Matrix<pm::Rational>, bool, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Matrix<Rational> > >(), arg1, arg2 );
   }
   FunctionWrapperInstance4perl( perl::Object (pm::Matrix<pm::Rational>, bool, int) );

   FunctionWrapper4perl( perl::Object (int, pm::IncidenceMatrix<pm::NonSymmetric>, bool, int) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      IndirectWrapperReturn( arg0, arg1.get< perl::TryCanned< const IncidenceMatrix< NonSymmetric > > >(), arg2, arg3 );
   }
   FunctionWrapperInstance4perl( perl::Object (int, pm::IncidenceMatrix<pm::NonSymmetric>, bool, int) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
