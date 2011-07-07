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
#include "polymake/IncidenceMatrix.h"
#include "polymake/SparseMatrix.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0>
   FunctionInterface4perl( col_x_f5, T0 ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      WrapperReturnLvalue( T0, arg0.get<T0>().col(arg1) );
   };

   FunctionInstance4perl(col_x_f5, perl::Canned< const Wary< pm::ColChain<pm::ColChain<pm::ColChain<pm::ColChain<pm::Matrix<pm::Rational> const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > > >);
   FunctionInstance4perl(col_x_f5, perl::Canned< const Wary< pm::MatrixMinor<pm::Matrix<pm::Rational> const&, pm::incidence_line<pm::AVL::tree<pm::sparse2d::traits<pm::sparse2d::traits_base<pm::nothing, true, false, (pm::sparse2d::restriction_kind)0>, false, (pm::sparse2d::restriction_kind)0> > const&> const&, pm::all_selector const&> > >);
   FunctionInstance4perl(col_x_f5, perl::Canned< const Wary< pm::ColChain<pm::ColChain<pm::ColChain<pm::Matrix<pm::Rational> const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > > >);
   FunctionInstance4perl(col_x_f5, perl::Canned< const Wary< pm::ColChain<pm::ColChain<pm::Matrix<pm::Rational> const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > const&, pm::SingleCol<pm::Vector<pm::Rational> const&> > > >);
   FunctionInstance4perl(col_x_f5, perl::Canned< const Wary< SparseMatrix< Rational, NonSymmetric > > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
