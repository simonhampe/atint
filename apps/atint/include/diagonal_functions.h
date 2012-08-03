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
 */

#include "polymake/client.h"

#ifndef DIAGONAL_FUNCTIONS_H
#define DIAGONAL_FUNCTIONS_H

namespace polymake { namespace atint {
  
  /**
   @brief This function lifts a tropical variety to one dimension higher by adding a last coordinate = 0 and a lineality space generator (1,...,1)
   @param perl::Object X A WeightedComplex
   @return perl::Object The lifted WeightedComplex
  */
  perl::Object lift_variety(perl::Object X);

  
  /**
   @brief This function mods out the lineality space (1,..,1) by setting the last coordinate to be minus the sum of the remaining coordinates
   @param perl::Object X A WeightedComplex
   @return perl::Object The projected WeightedComplex
  */
  perl::Object project_variety(perl::Object X);
}}
#endif // DIAGONAL_FUNCTIONS_h