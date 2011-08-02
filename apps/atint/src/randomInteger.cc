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
 Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>
 */

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/RandomGenerators.h"

namespace polymake { namespace atint {

  Array<Integer> randomInteger(const int& max_arg, const int &n) {
    static Integer upperBound = 0;
    static UniformlyRandomRanged<Integer> rg(max_arg);
      if(max_arg != upperBound)  {
	rg = UniformlyRandomRanged<Integer>(max_arg);
	upperBound = max_arg;
      }
    Array<Integer> result(n);
    for(int i = 0; i < n; i++) {
	result[i] = rg.get();
    }
    return result;
  }
  
  
  
  UserFunction4perl("# @category Random number generators"
                  "# Returns n random integers in the range 0.. (max_arg-1),inclusive"
                  "# Note that this algorithm is not optimal for real randomness:"
                  "# If you change the range parameter and then change it back, you will"
                  "# usually get the exact same sequence"
                  "# @param int max_arg The upper bound for the random integers"
                  "# @param int n The number of integers to be created"
                  "# @return Array<Integer>",
                  &randomInteger,"randomInteger($, $)");                  
  
}}