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
 
 This file includes functionalities to compute properties of WeightedComplex objects
 that should be available on the c++ side
 */

#include "polymake/client.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"

#ifndef WEIGHTED_COMPLEX_RULES_H
#define WEIGHTED_COMPLEX_RULES_H

namespace polymake { namespace atint { 
  
struct CodimensionOneResult {
  IncidenceMatrix<> codimOneCones;
  IncidenceMatrix<> codimOneInMaximal;
};

/**
 @brief Computes the codimension one data of a polyhedral complex
 @param Matrix<Rational> rays The ray matrix of the complex
 @param IncidenceMatrix<> maximalCones The incidence matrix of the maximal cells
 @param bool uses_homog Whether the rays are given in homogeneous coordinates
 @param Matrix<Rational> linspace The generators of the lineality space
 @param Array<Set<int>> local_restriction The local restriction of the complex
 @return CodimensionOneResult A struct containing:
 1) An IncidenceMatrix<> called codimOneCones describing the codimension one cells in terms of the rays
 2) An IncidenceMatrix<> called codimOneInMaximal describing which codim one cell lies in which maximal cells
*/
CodimensionOneResult calculateCodimOneData(Matrix<Rational> rays, IncidenceMatrix<> maximalCones, bool uses_homog, Matrix<Rational> linspace, Array<Set<int> > local_restriction);

/**
 @brief Check whether a given cone set is compatible with a given set of local restrictions
 @param Set<int> cone A set of (ray) indices
 @param Array<Set<int>>  local_restriction A list of sets of ray indices
 @return true, if and only if cone is either contained in one of the sets of local_restriction or
 contains one of them
 */
 bool is_coneset_compatible(const Set<int> &cone, const Array<Set<int> > &local_restriction);

}}


#endif // WEIGHTED_COMPLEX_RULES_H