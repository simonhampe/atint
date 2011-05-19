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

#ifndef ATINT_DIVISOR_H
#define ATINT_DIVISOR_H

namespace polymake { namespace tropical {

/**
@brief Takes two fans and computes the intersection of both. The function relies on the fact that the latter fan is complete (i.e. its support is the whole ambient space) to compute the intersection correctly.
@param fan An arbitrary polyhedral fan
@param completeFan A complete polyhedral fan
@return fan::PolyhedralFan The intersection of both fans (whose support is equal to the support of fan). The 
resulting fan uses homogeneous coordinates if and only fan does. If fan has a property TROPICAL_WEIGHTS, 
the tropical weights of the refinement are also computed. If fan is zero-dimensional (i.e. a point), fan is returned.
*/
perl::Object intersect_complete_fan(perl::Object fan, perl::Object completeFan);


}}

#endif