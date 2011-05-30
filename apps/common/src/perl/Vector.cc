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

///==== this line controls the automatic file splitting: max.instances=40

#include "polymake/client.h"
#include "polymake/Vector.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/linalg.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   Class4perl("Polymake::common::Vector__Vector__Integer", Vector< Vector< Integer > >);
   OperatorInstance4perl(BinaryAssign_add, perl::Canned< Wary< Vector< Integer > > >, perl::Canned< const Vector< Rational > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< Vector< Integer > > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< pm::SameElementVector<pm::Rational> > >, perl::Canned< const Vector< Rational > >);
   OperatorInstance4perl(Binary__or, int, perl::Canned< const Vector< Integer > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< pm::VectorChain<pm::SingleElementVector<pm::Integer>, pm::Vector<pm::Integer> const&> > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const Matrix< Rational > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational>&>, pm::Series<int, true>, void> > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const Vector< Rational > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational>&>, pm::Series<int, true>, void> > >, int);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const pm::RowChain<pm::RowChain<pm::Matrix<pm::Rational> const&, pm::SingleRow<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&> > const&, pm::SingleRow<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&> > >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< pm::SameElementVector<pm::Integer> > >, perl::Canned< const Vector< Rational > >);
   Class4perl("Polymake::common::Vector__Set__Int", Vector< Set< int > >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational const&> > >, int);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, int);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > >, int);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational const&> > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< Vector< Rational > > >, int);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, false>, void>, pm::Vector<pm::Rational> const&>, void> >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const Matrix< Rational > >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SameElementVector<pm::Rational> const&> >);
   OperatorInstance4perl(Binary__eq, perl::Canned< const Wary< pm::IndexedSlice<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, false>, void>, pm::incidence_line<pm::AVL::tree<pm::sparse2d::traits<pm::sparse2d::traits_base<pm::nothing, true, false, (pm::sparse2d::restriction_kind)0>, false, (pm::sparse2d::restriction_kind)0> > const&> const&, void> > >, perl::Canned< const pm::SameElementVector<pm::Rational> >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const Rational >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::IndexedSlice<pm::Vector<pm::Rational>&, pm::Array<int, void> const&, void> >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::IndexedSlice<pm::Vector<pm::Rational>&, pm::Array<int, void> const&, void> >, int);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::IndexedSlice<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&, pm::Array<int, void> const&, void> >, int);
   OperatorInstance4perl(Binary_sub, perl::Canned< const Wary< pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> > >, perl::Canned< const pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void>, pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&, pm::Vector<pm::Rational> const&> >, void> >);
   FunctionInstance4perl(new_X, Vector< double >, perl::Canned< const pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<double>&>, pm::Series<int, true>, void> >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<double> const&>, pm::Series<int, true>, void>, pm::Vector<double> const&>, void> > >, perl::Canned< const pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<double> const&>, pm::Series<int, true>, void>, pm::Vector<double> const&>, void> >);
   FunctionInstance4perl(new_X, Vector< double >, perl::Canned< const Vector< double > >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void>, pm::Vector<pm::Rational> const&>, void> > >, perl::Canned< const pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void>, pm::Vector<pm::Rational> const&>, void> >);
   OperatorInstance4perl(Binary_mul, int, perl::Canned< const Wary< pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void>, pm::Vector<pm::Rational> const&>, void> > >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< pm::ContainerUnion<pm::cons<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void>, pm::Vector<pm::Rational> const&>, void> > >, perl::Canned< const Vector< Rational > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
