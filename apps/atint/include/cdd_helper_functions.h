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
 
 This defines wrapper / helper functions for using the polymake-cdd-interface
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"

#ifndef CDD_HELPER_FUNCTIONS_H
#define CDD_HELPER_FUNCTIONS_H

namespace polymake { namespace atint { 

  struct fan_intersection_result {
    Matrix<Rational> rays;
    Matrix<Rational> lineality_space;
    IncidenceMatrix<> cones;
    IncidenceMatrix<> xcontainers;
    IncidenceMatrix<> ycontainers;
  };
  
  
  /**
   @brief Inserts a list of new rays into a matrix of existing rays, taking care of doubles. This method DOES change the matrix of existing rays.
   @param Matrix<Rational> rays A reference to a matrix of (normalized) rays or vertices. This matrix will potentially be changed
   @param Matrix<Rational> nrays A list of new rays or vertices (not necessarily normalized), that will be added
   @param bool is_normalized Whether the new rays are also already normalized
   @param bool uses_homog Whether both matrices are given in homogeneous coordinates or not
   @return Vector<int> At position i contains the row index of the new ray nrays[i] in the modified matrix rays
   */
  Vector<int> insert_rays(Matrix<Rational> &rays, Matrix<Rational> nrays, bool is_normalized, bool uses_homog);
  
  /**
   @brief Normalizes a ray matrix: Vertices begin with a 1 and the first non-zero coordinate of a ray is +-1
   @param Matrix<Rational> The row vectors to be normalized. This method modifies this matrix!
   @param bool uses_homog Whether the rays are in homog. coordinates. If so, then every row starting with a nonzero entry is a vertex
   */
  void cdd_normalize_rays(Matrix<Rational> &rays, bool uses_homog);
  
  /**
   @brief Computes the intersection of two rational polyhedra, x and y, given in terms of rays and lineality space
   @param Matrix<Rational> xrays The rays of x
   @param Matrix<Rational> xlin The lineality space of x
   @param Matrix<Rational> yrays The rays of y
   @param Matrix<Rational> ylin The lineality space of y
   @param bool uses_homog Whether all rays and lineality generators are given in homogeneous coordinates. 
   @return An std::pair of Matrix<Rational>, containing first rays, then lineality space of the intersection. The rays are normalized according to whether the input was homogeneous or not
   */
  std::pair<Matrix<Rational>, Matrix<Rational> > cdd_cone_intersection(
      const Matrix<Rational> &xrays, const Matrix<Rational> &xlin, 
      const Matrix<Rational> &yrays, const Matrix<Rational> &ylin, bool uses_homog);

  /**
   @brief Computes the set-theoretic intersection of two polyhedral complexes, x and y.
   @param Matrix<Rational> xrays The rays of the complex x
   @param Matrix<Rational> xlin The lineality space of the complex x
   @param IncidenceMatrix<> xcones The cones of the complex x, in terms of xrays
   @param Matrix<Rational> yrays The rays of the complex y
   @param Matrix<Rational> ylin The lineality space of the complex y
   @param IncidenceMatrix<> ycones The cones of the complex y, in terms of yrays
   @param bool uses_homog Whether the result should be given in homogeneous coordinates. Regardless of homogeneity all matrices above should have the same number of columns. This only influences how the rays are normalized
   @return The set-theoretic intersection of x and y in a struct containing:
   1) Matrix<Rational> rays The rays of the intersection complex. The result is in homogeneous coordinates, if and only if one of x or y is in homog. coordinates.
   2) Matrix<Rational> lineality The lineality space of the complex
   3) IncidenceMatrix<> cones The maximal cones of the complex. More precisely: All pairwise intersections s \cap t of cones s in x, t in y. In particular, cones in this matrix might be contained in one another (but never equal). This list will not contain the empty cone, if the result should be in homogeneous coordinates (but might contain it otherwise)
   4) IncidenceMatrix<> xcontainers For the i-th cone of the intersection, the i-th row gives the indices of all maximal cones in x containing it
   5) IncidenceMatrix<> ycontainers For the i-th cone of the intersection, the i-th row gives the indices of all maximal cones in y containing it
   */
  fan_intersection_result cdd_fan_intersection(	
      Matrix<Rational> xrays, Matrix<Rational> xlin, IncidenceMatrix<> xcones,
      Matrix<Rational> yrays, Matrix<Rational> ylin, IncidenceMatrix<> ycones,
      bool uses_homog);
  
  
}}
  
#endif // CDD_HELPER_FUNCTIONS_H