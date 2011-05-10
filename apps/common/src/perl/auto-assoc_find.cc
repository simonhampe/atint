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
#include "polymake/perl/assoc.h"
#include "polymake/Map.h"
#include "polymake/Vector.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( assoc_find_X_X, T0,T1 ) {
      perl::Value arg0(stack[0]), arg1(stack[1]);
      WrapperReturn( pm::perl::find_element(arg0.get<T0>(), arg1.get<T1>()) );
   };

   FunctionInstance4perl(assoc_find_X_X, perl::Canned< const Map< int, Vector< Integer > > >, int);
   FunctionInstance4perl(assoc_find_X_X, perl::Canned< const Map< int, Vector< Rational > > >, int);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
