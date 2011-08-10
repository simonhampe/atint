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

Contains the definitions for basicoperations.cc
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"

#ifndef ATINT_BASICOP_H
#define ATINT_BASICOP_H

namespace polymake { namespace atint{
  
  /**
    @brief Takes a matrix and returns the row indices where the first coordinate is nonzero and where the first coordinate is zero in  two different sets
    @param Matrix<Rational> m The matrix whose rows we consider
    @param Set<int> affine A reference to a set where the indices of the affine rows are stored
    @param Set<int> directional A reference to a set where the indices of the directional rows are stored
    @param bool uses_homog When set to false, all rows are directional and affine only contains the index -1
    @return std::pair<Set<int>, Set<int> > The first set contains the row indices of rows that start with a nonzero entry, the second set is the complement
  */
  inline std::pair<Set<int>, Set<int> > separateRays(Matrix<Rational> m, Set<int> &affine, Set<int> &directional, bool uses_homog);
  
  /**
    @brief Takes a list of WeightedComplex objects (that may be weighted and that may use homogeneous coordinates) and computes the cartesian product of these. If any complex uses homogeneous coordinates, so will the result. If any complex has weights, all non-weighted complexes will be treated as having constant weight 1.
    @return The cartesian product of the complexes
  */
  perl::Object compute_product_complex(std::vector<perl::Object> complexes) ;
  
  /**
   @brief Takes a weighted complex and computes the refinement of the fan along the halfspace fans given by the rows of a matrix
   @param perl::Object fan The weighted complex to be refined, given in homogeneous coordinates.
   @param Matrix<Rational> facets A matrix whose rows define halfspace fans in the following way: The matrix should have fan->CMPLX_AMBIENT_DIM+1 columns. The function will then intersect with the halfspace complexes defined by (0,row) for each row.
   @return perl::Object The corr. refinement of fan
  */
  perl::Object facetRefinement(perl::Object fan, Matrix<Rational> facets);
  
  /**
   @brief Takes a polyhedral complex and returns a list of all the local vertex fans, i.e. for each affine ray r, the list contains the fan Star_complex(r) (in non-homogeneous coordinates)
   @param WeightedComplex complex A tropical variety
   @return perl::ListReturn A list of WeightedComplex objects in non-homogeneous coordinates. The i-th complex corresponds to the i-th affine ray ( vertex). If the complex is not in homogeneous coordinates, the list contains just the complex itself
   */
  perl::ListReturn fan_decomposition(perl::Object complex);
  
  /**
   @brief Takes a polyhedral complex and applies an affine linear transformation, given by a translate vector and a matrix. The method assumes the function is bijective and preserves cones, i.e. it just applies the transformation to the rays and lineality space and leaves the cones and weights unchanged.
   @param WeightedComplex complex The complex to be transformed, supposed to be in homogeneous coordinates (if translate != 0)
   @param Vector<Rational> translate A vector whose dimension should be equal to the column dimension of the transformation matrix
   @param Matrix<Integer> matrix An integer (square) invertible matrix. 
   @return WeightedComplex The transformed complex
   */
  perl::Object affineTransformation(perl::Object complex, Vector<Rational> translate, Matrix<Integer> matrix);
}}

#endif