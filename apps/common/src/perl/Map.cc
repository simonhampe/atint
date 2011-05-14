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
   template <typename T0>
   FunctionInterface4perl( new, T0 ) {
      WrapperReturnNew(T0, () );
   };

   Class4perl("Polymake::common::Map_A_Pair_A_Int_I_Int_Z_I_Vector__Integer_Z", Map< std::pair< int, int >, Vector< Integer > >);
   OperatorInstance4perl(Binary_brk, perl::Canned< const Map< std::pair< int, int >, Vector< Integer > > >, perl::TryCanned< const std::pair< int, int > >);
   Class4perl("Polymake::common::Map_A_Int_I_Map_A_Int_I_Vector__Integer_Z_Z", Map< int, Map< int, Vector< Integer > > >);
   Class4perl("Polymake::common::Map_A_Int_I_Vector__Integer_Z", Map< int, Vector< Integer > >);
   FunctionInstance4perl(new, Map< int, Map< int, Vector< Integer > > >);
   FunctionInstance4perl(new, Map< int, Vector< Integer > >);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, Map< int, Vector< Integer > > > >, int);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, Vector< Integer > > >, int);
   Class4perl("Polymake::common::Map_A_Int_I_Map_A_Int_I_Vector__Rational_Z_Z", Map< int, Map< int, Vector< Rational > > >);
   Class4perl("Polymake::common::Map_A_Int_I_Vector__Rational_Z", Map< int, Vector< Rational > >);
   FunctionInstance4perl(new, Map< int, Map< int, Vector< Rational > > >);
   FunctionInstance4perl(new, Map< int, Vector< Rational > >);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, Vector< Rational > > >, int);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, Map< int, Vector< Rational > > > >, int);
   OperatorInstance4perl(Binary_brk, perl::Canned< const Map< int, Map< int, Vector< Rational > > > >, int);
   Class4perl("Polymake::common::Map_A_Int_I_Map_A_Int_I_String_Z_Z", Map< int, Map< int, std::string > >);
   Class4perl("Polymake::common::Map_A_Int_I_String_Z", Map< int, std::string >);
   FunctionInstance4perl(new, Map< int, Map< int, std::string > >);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, Map< int, std::string > > >, int);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, std::string > >, int);
   Class4perl("Polymake::common::Map_A_String_I_Rational_Z", Map< std::string, Rational >);
   Class4perl("Polymake::common::Map_A_Int_I_Map_A_String_I_Rational_Z_Z", Map< int, Map< std::string, Rational > >);
   FunctionInstance4perl(new, Map< int, Map< std::string, Rational > >);
   FunctionInstance4perl(new, Map< std::string, Rational >);
   OperatorInstance4perl(Binary_brk, perl::Canned< Map< int, Map< std::string, Rational > > >, int);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
