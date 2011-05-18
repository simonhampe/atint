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
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/linalg.h"
#include "polymake/Integer.h"
#include "polymake/Set.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::MatrixMinor<pm::ColChain<pm::ColChain<pm::ColChain<pm::ColChain<pm::Matrix<pm::Rational> const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::all_selector const&, pm::Array<int, void> const&> >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::ColChain<pm::SingleCol<pm::Vector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< Matrix< Rational > > >, perl::Canned< const pm::ColChain<pm::SingleCol<pm::SameElementVector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::RowChain<pm::Matrix<pm::Rational> const&, pm::ColChain<pm::SingleCol<pm::SameElementVector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> const&> >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< Matrix< Rational > > >, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SameElementVector<pm::Rational> const&> >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::RowChain<pm::Matrix<pm::Rational> const&, pm::SingleRow<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SameElementVector<pm::Rational> const&> const&> > >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::RowChain<pm::RowChain<pm::RowChain<pm::Matrix<pm::Rational> const&, pm::SingleRow<pm::Vector<pm::Rational> const&> > const&, pm::SingleRow<pm::Vector<pm::Rational> const&> > const&, pm::SingleRow<pm::Vector<pm::Rational> const&> > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< Matrix< Rational > > >, perl::Canned< const pm::IndexedSlice<pm::Vector<pm::Rational>&, pm::Array<int, void> const&, void> >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< Matrix< Rational > > >, perl::Canned< const pm::VectorChain<pm::IndexedSlice<pm::Vector<pm::Rational>&, pm::Array<int, void> const&, void> const&, pm::SingleElementVector<pm::Rational> > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< Matrix< Rational > > >, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > >);
   OperatorInstance4perl(Unary_neg, perl::Canned< const Wary< Matrix< Rational > > >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::MatrixMinor<pm::ColChain<pm::ColChain<pm::ColChain<pm::Matrix<pm::Rational> const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::all_selector const&, pm::Array<int, void> const&> >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::MatrixMinor<pm::ColChain<pm::ColChain<pm::Matrix<pm::Rational> const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::all_selector const&, pm::Array<int, void> const&> >);
   OperatorInstance4perl(Unary_neg, perl::Canned< const Wary< pm::ColChain<pm::SingleCol<pm::Vector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> > >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< Matrix< Rational > > >, perl::Canned< const pm::VectorChain<pm::IndexedSlice<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&, pm::Array<int, void> const&, void> const&, pm::SingleElementVector<pm::Rational> > >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::RowChain<pm::RowChain<pm::RowChain<pm::Matrix<pm::Integer> const&, pm::SingleRow<pm::Vector<pm::Integer> const&> > const&, pm::SingleRow<pm::Vector<pm::Integer> const&> > const&, pm::SingleRow<pm::Vector<pm::Integer> const&> > >);
   FunctionInstance4perl(new_X, Matrix< Rational >, perl::Canned< const pm::RowChain<pm::Matrix<pm::Rational> const&, pm::MatrixMinor<pm::Matrix<pm::Rational> const&, pm::Set<int, pm::operations::cmp> const&, pm::all_selector const&> const&> >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
