/* Copyright (c) 1997-2014
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

///==== this line controls the automatic file splitting: max.instances=40

#include "polymake/IncidenceMatrix.h"
#include "polymake/Integer.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Set.h"
#include "polymake/SparseMatrix.h"
#include "polymake/Vector.h"
#include "polymake/client.h"
#include "polymake/linalg.h"

namespace polymake { namespace common { namespace {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   Class4perl("Polymake::common::Vector__Set__Int", Vector< Set< int > >);
   Class4perl("Polymake::common::Vector__Bool", Vector< bool >);
   OperatorInstance4perl(Binary__ora, perl::Canned< const Vector< Rational > >, perl::Canned< const pm::SameElementVector<pm::Rational const&> >);
   OperatorInstance4perl(Binary__ora, perl::Canned< const Vector< Rational > >, perl::Canned< const pm::ColChain<pm::SingleCol<pm::Vector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> >);
   Class4perl("Polymake::common::Vector__String", Vector< std::string >);
   OperatorInstance4perl(Binary__ora, perl::Canned< const Vector< Rational > >, perl::Canned< const SparseMatrix< Rational, NonSymmetric > >);
   OperatorInstance4perl(Binary__ora, perl::Canned< const Vector< Rational > >, perl::Canned< const pm::ColChain<pm::SingleCol<pm::Vector<pm::Rational> const&>, pm::SparseMatrix<pm::Rational, pm::NonSymmetric> const&> >);
   FunctionInstance4perl(new_X, Vector< Integer >, perl::Canned< const pm::SameElementVector<pm::Integer const&> >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, false>, void> >);
   OperatorInstance4perl(Binary__ora, perl::Canned< const pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, false>, void> >, perl::Canned< const Matrix< Rational > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
