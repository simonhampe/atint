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

#ifndef ATINT_MORPHISM_SPECIAL_H
#define ATINT_MORPHISM_SPECIAL_H


namespace polymake { namespace atint {
  
  
  //Documentation see perl wrapper
  perl::Object evaluation_map(int n, int r, Matrix<Rational> delta, int i);
  
  //Documentation see perl wrapper
  perl::Object evaluation_map_d(int n, int r, int d, int i);
  
  
  
}}

#endif // ATINT_MORPHISM_SPECIAL_H