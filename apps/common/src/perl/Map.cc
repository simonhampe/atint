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
#include "polymake/Map.h"
#include "polymake/Vector.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   template <typename T0>
   FunctionInterface4perl( new, T0 ) {
      WrapperReturnNew(T0, () );
   };

   Class4perl("Polymake::common::Map_A_Int_I_Map_A_Int_I_Vector__Integer_Z_Z", Map< int, Map< int, Vector< Integer > > >);
   Class4perl("Polymake::common::Map_A_Int_I_Vector__Integer_Z", Map< int, Vector< Integer > >);
   Class4perl("Polymake::common::Map_A_Int_I_Map_A_Int_I_Vector__Rational_Z_Z", Map< int, Map< int, Vector< Rational > > >);
   Class4perl("Polymake::common::Map_A_Int_I_Vector__Rational_Z", Map< int, Vector< Rational > >);
   OperatorInstance4perl(Binary_brk, perl::Canned< const Map< int, Map< int, Vector< Integer > > > >, int);
   Class4perl("Polymake::common::Map_A_Int_I_Map_A_Int_I_String_Z_Z", Map< int, Map< int, std::string > >);
   Class4perl("Polymake::common::Map_A_Int_I_String_Z", Map< int, std::string >);
   FunctionInstance4perl(new, Map< int, Map< int, std::string > >);
   FunctionInstance4perl(new_X, Map< int, Map< int, Vector< Integer > > >, perl::Canned< const Map< int, Map< int, Vector< Integer > > > >);
   FunctionInstance4perl(new_X, Map< int, Map< int, Vector< Rational > > >, perl::Canned< const Map< int, Map< int, Vector< Rational > > > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
