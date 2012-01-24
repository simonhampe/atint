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
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/linalg.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
#include "polymake/Integer.h"
#include "polymake/SparseMatrix.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0>
   FunctionInterface4perl( new, T0 ) {
      WrapperReturnNew(T0, () );
   };

   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   Class4perl("Polymake::common::Vector__Set__Int", Vector< Set< int > >);
   Class4perl("Polymake::common::Vector__String", Vector< std::string >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const pm::SameElementVector<pm::Rational> >);
   OperatorInstance4perl(Binary_mul, perl::Canned< const Wary< pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> > >, perl::Canned< const Rational >);
   OperatorInstance4perl(assign, pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational>&>, pm::Series<int, true>, void>, perl::Canned< const pm::incidence_line<pm::AVL::tree<pm::sparse2d::traits<pm::sparse2d::traits_base<pm::nothing, true, false, (pm::sparse2d::restriction_kind)0>, false, (pm::sparse2d::restriction_kind)0> > const&> >);
   OperatorInstance4perl(convert, Vector< Set< int > >, perl::Canned< const IncidenceMatrix< NonSymmetric > >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, false>, void> >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const Array< Integer > >);
   FunctionInstance4perl(new_X, Vector< Integer >, perl::Canned< const Vector< int > >);
   OperatorInstance4perl(convert, Vector< Set< int > >, perl::Canned< const Array< Set< int > > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::SameElementVector<pm::Rational> > >, perl::Canned< const pm::DiagMatrix<pm::SameElementVector<pm::Rational>, true> >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const SparseMatrix< Rational, NonSymmetric > >);
   OperatorInstance4perl(Binary__ne, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const pm::SameElementVector<pm::Rational> >);
   FunctionInstance4perl(new, Vector< Set< int > >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const Rational >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational const&> > >, perl::Canned< const Rational >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational const&> > const&, pm::SingleElementVector<pm::Rational const&> > >);
   OperatorInstance4perl(Binary_add, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const Vector< Rational > >);
   OperatorInstance4perl(Binary_mul, int, perl::Canned< const Wary< Vector< Rational > > >);
   OperatorInstance4perl(Binary_sub, perl::Canned< const Wary< Vector< Rational > > >, perl::Canned< const Vector< Rational > >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::Vector<pm::Rational> const&> >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const pm::ColChain<pm::SingleCol<pm::Vector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> >);
   FunctionInstance4perl(new_X, Vector< Rational >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, int);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > >, int);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   OperatorInstance4perl(Binary__or, perl::Canned< const pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::VectorChain<pm::Vector<pm::Rational> const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > const&, pm::SingleElementVector<pm::Rational> > >, perl::Canned< const Integer >);
   FunctionInstance4perl(new_X, Vector< Integer >, perl::Canned< const Vector< Integer > >);
   FunctionInstance4perl(new_X, Vector< int >, perl::Canned< const Vector< int > >);
   OperatorInstance4perl(Binary__or, perl::Canned< const Vector< Rational > >, perl::Canned< const Matrix<Rational> >);
   Class4perl("Polymake::common::Vector__Vector__Int", Vector< Vector< int > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
