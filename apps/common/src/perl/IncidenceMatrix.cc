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
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   FunctionInstance4perl(new_X, IncidenceMatrix< NonSymmetric >, perl::Canned< const Vector< Set< int > > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > > >, perl::Canned< const Set< int > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::RowChain<pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > > >, perl::Canned< const Set< int > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > > >, perl::Canned< const Set< int > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > > >, perl::Canned< const Set< int > >);
   FunctionInstance4perl(new_X, IncidenceMatrix< NonSymmetric >, perl::Canned< const pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > >);
   FunctionInstance4perl(new_X, IncidenceMatrix< NonSymmetric >, perl::Canned< const pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > >);
   FunctionInstance4perl(new_X, IncidenceMatrix< NonSymmetric >, perl::Canned< const pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > > >, perl::Canned< const Set< int > >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::RowChain<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > const&, pm::SingleIncidenceRow<pm::Set_with_dim<pm::Set<int, pm::operations::cmp> const&> > > > >, perl::Canned< const Set< int > >);
   FunctionInstance4perl(new_X, IncidenceMatrix< NonSymmetric >, perl::Canned< const pm::MatrixMinor<pm::IncidenceMatrix<pm::NonSymmetric> const&, pm::Set<int, pm::operations::cmp> const&, pm::all_selector const&> >);
   OperatorInstance4perl(Binary_div, perl::Canned< const Wary< IncidenceMatrix< NonSymmetric > > >, perl::Canned< const Set< int > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
