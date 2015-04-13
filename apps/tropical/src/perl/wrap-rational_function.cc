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

#include "polymake/TropicalNumber.h"

namespace polymake { namespace tropical { namespace {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0>
   FunctionInterface4perl( computeDomain_T_x_f16, T0 ) {
      perl::Value arg0(stack[0]);
      WrapperReturnVoid( (computeDomain<T0>(arg0)) );
   };

   template <typename T0>
   FunctionInterface4perl( computeGeometricFunctionData_T_x_f16, T0 ) {
      perl::Value arg0(stack[0]);
      WrapperReturnVoid( (computeGeometricFunctionData<T0>(arg0)) );
   };

   FunctionInstance4perl(computeGeometricFunctionData_T_x_f16, Min);
   FunctionInstance4perl(computeGeometricFunctionData_T_x_f16, Max);
   FunctionInstance4perl(computeDomain_T_x_f16, Min);
   FunctionInstance4perl(computeDomain_T_x_f16, Max);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} } }
