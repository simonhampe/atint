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
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
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
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
