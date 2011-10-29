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
#include "polymake/linalg.h"
#include "polymake/Rational.h"
#include "polymake/SparseMatrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/Set.h"
#include "polymake/SparseVector.h"
#include "polymake/Vector.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   OperatorInstance4perl(Binary__or, perl::Canned< const Wary< pm::RowChain<pm::SingleRow<pm::SameElementVector<pm::Rational> const&>, pm::DiagMatrix<pm::SameElementVector<pm::Rational>, true> const&> > >, perl::Canned< const pm::SameElementVector<pm::Rational> >);
   FunctionInstance4perl(new_X, SparseMatrix< Rational, NonSymmetric >, perl::Canned< const pm::ColChain<pm::RowChain<pm::SingleRow<pm::SameElementVector<pm::Rational> const&>, pm::DiagMatrix<pm::SameElementVector<pm::Rational>, true> const&> const&, pm::SingleCol<pm::SameElementVector<pm::Rational> const&> > >);
   FunctionInstance4perl(new_X, SparseMatrix< Rational, NonSymmetric >, perl::Canned< const pm::MatrixMinor<pm::SparseMatrix<pm::Rational, pm::NonSymmetric> const&, pm::all_selector const&, pm::Complement<pm::SingleElementSet<int const&>, int, pm::operations::cmp> const&> >);
   OperatorInstance4perl(Unary_neg, perl::Canned< const Wary< pm::ColChain<pm::SingleCol<pm::SparseVector<pm::Rational, pm::conv<pm::Rational, bool> > const&>, pm::SparseMatrix<pm::Rational, pm::NonSymmetric> const&> > >);
   FunctionInstance4perl(new_X, SparseMatrix< Rational, NonSymmetric >, perl::Canned< const pm::ColChain<pm::SingleCol<pm::Vector<pm::Rational> const&>, pm::SparseMatrix<pm::Rational, pm::NonSymmetric> const&> >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Wary< pm::DiagMatrix<pm::SameElementVector<pm::Rational>, true> > >, perl::Canned< const Matrix<Rational> >);
   FunctionInstance4perl(new_X, SparseMatrix< int, NonSymmetric >, perl::Canned< const SparseMatrix< int, NonSymmetric > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
