/*
 T his program is free s*oftware; you can redistribute it and/or
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
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 Defines basic operations on MinMaxFunctions
 */

#include "polymake/client.h"

#ifndef MINMAX_FUNCTIONS_H
#define MINMAX_FUNCTIONS_H

namespace polymake { namespace atint {
 
  /**
    @brief Computes the sum of two MinMaxFunctions that are supposed to have the same USES_MIN property
    @param MinMaxFunction f An arbitrary MinMaxFunction
    @param MinMaxFunction g A MinMaxFunction such that f->[[USES_MIN]] == g->[[USES_MIN]] and is defined on the same domain.
    @param Bool reduce Optional. False by default. If true, the function reduces the amount of terms by computing vertices of the newton polytope of the sum
    @return MinMaxFunction The sum of both functions.
    */
  perl::Object add_minmax_functions(perl::Object f, perl::Object g, bool reduce = false);
  
  /**
    @brief Scales a MinMaxFunction by a given Rational a
    @param MinMaxFunction f An arbitrary MinMaxFunction
    @param Rational a A scalar values
    @return MinMaxFunction The scaled function
    */
  perl::Object scale_minmax_function(perl::Object f, Rational a);
  
  
}
}

#endif // MINMAX_FUNCTIONS_H