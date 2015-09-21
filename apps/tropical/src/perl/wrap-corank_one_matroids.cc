/* Copyright (c) 1997-2015
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

namespace polymake { namespace tropical { namespace {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   FunctionWrapper4perl( pm::perl::ListReturn (int) ) {
      perl::Value arg0(stack[0]);
      IndirectWrapperReturnVoid( arg0 );
   }
   FunctionWrapperInstance4perl( pm::perl::ListReturn (int) );

   FunctionWrapper4perl( pm::Set<pm::Set<int, pm::operations::cmp>, pm::operations::cmp> (pm::Set<pm::Set<int, pm::operations::cmp>, pm::operations::cmp> const&, pm::Set<pm::Set<int, pm::operations::cmp>, pm::operations::cmp> const&) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0.get< perl::TryCanned< const Set< Set< int > > > >(), arg1.get< perl::TryCanned< const Set< Set< int > > > >() );
   }
   FunctionWrapperInstance4perl( pm::Set<pm::Set<int, pm::operations::cmp>, pm::operations::cmp> (pm::Set<pm::Set<int, pm::operations::cmp>, pm::operations::cmp> const&, pm::Set<pm::Set<int, pm::operations::cmp>, pm::operations::cmp> const&) );

   FunctionWrapper4perl( perl::Object (int, pm::Set<int, pm::operations::cmp>) ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      IndirectWrapperReturn( arg0, arg1.get< perl::TryCanned< const Set< int > > >() );
   }
   FunctionWrapperInstance4perl( perl::Object (int, pm::Set<int, pm::operations::cmp>) );

///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
