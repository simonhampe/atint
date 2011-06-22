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
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0>
   FunctionInterface4perl( null_space_X, T0 ) {
      perl::Value arg0(stack[0]);
      WrapperReturn( null_space(arg0.get<T0>()) );
   };

   FunctionInstance4perl(null_space_X, perl::Canned< const pm::RowChain<pm::SingleRow<pm::Vector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> >);
   FunctionInstance4perl(null_space_X, perl::Canned< const pm::Transposed<pm::RowChain<pm::SingleRow<pm::Vector<pm::Rational> const&>, pm::Matrix<pm::Rational> const&> > >);
   FunctionInstance4perl(null_space_X, perl::Canned< const pm::Transposed<pm::RowChain<pm::SingleRow<pm::Vector<pm::Rational> const&>, pm::RowChain<pm::RowChain<pm::Matrix<pm::Rational> const&, pm::SingleRow<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&> > const&, pm::SingleRow<pm::IndexedSlice<pm::masquerade<pm::ConcatRows, pm::Matrix_base<pm::Rational> const&>, pm::Series<int, true>, void> const&> > const&> > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
