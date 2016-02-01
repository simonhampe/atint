/*
	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License
	as published by the Free Software Foundation; either version 2
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor,
	Boston, MA  02110-1301, USA.

	---
	Copyright (C) 2011 - 2015, Simon Hampe <simon.hampe@googlemail.com>

	Contains functions to compute the affine transform of a cycle 
	*/

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/IncidenceMatrix.h"

namespace polymake { namespace tropical {

   template <typename Addition>
      perl::Object dual_addition_version(perl::Object ring_cycle) {
         int n = ring_cycle.give("N_ELEMENTS");
         int r = ring_cycle.give("RANK");
         Array<IncidenceMatrix<> > pres = ring_cycle.give("NESTED_PRESENTATIONS");
         Array<int> coef = ring_cycle.give("NESTED_COEFFICIENTS");

         perl::Object result(perl::ObjectType::construct<typename Addition::dual>("MatroidRingCycle"));
            result.take("N_ELEMENTS") << n;
            result.take("RANK") << r;
            result.take("NESTED_PRESENTATIONS") << pres;
            result.take("NESTED_COEFFICIENTS") << coef;
         return result;
      }
  
   UserFunctionTemplate4perl("# @category Conversion of tropical addition"
         "# Takes a MatroidRingCycle and converts it to the dual tropical addition"
         "# @param MatroidRingCycle<Addition> M"
         "# @return MatroidRingCycle",
         "dual_addition_version<Addition>(MatroidRingCycle<Addition>)");

}}
