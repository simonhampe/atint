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
 Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"

#ifndef ATINT_REFINEMENT_H
#define ATINT_REFINEMENT_H

namespace polymake { namespace atint { 
  
/**
 @brief Takes a collection of weighted cones of the same dimension and refines them such that they form a weighted complex. Weights of cones lying one over the other add up.
 @param Matrix<Rational> rays A matrix of rays of the cones. Has to be irredundant and normalized.
 @param Vector<Set<int> > cones The maximal cones, given in terms of row indices of their rays
 @param Vector<Integer> weights The i-th element is the weight of the cone described by the i-th row of cones
 @param bool uses_homog Whether the rays are given in homogeneous coordinates
 @return perl::Object A WeightedComplex object, whose support is the union of the cones. The weight of a cone is the sum of the weights of the orginal cones containing it. The result is homogeneous iff uses_homog is true
 */
perl::Object complexify(Matrix<Rational> rays, Vector<Set<int> > max_cones, Vector<Integer> weights, bool uses_homog); 
}}
#endif // ATINT_REFINEMENT_H