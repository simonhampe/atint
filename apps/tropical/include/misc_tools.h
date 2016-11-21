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

	Contains some miscellaneous tools.

*/


#ifndef POLYMAKE_ATINT_MISC_TOOLS_H
#define POLYMAKE_ATINT_MISC_TOOLS_H

#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
#include "polymake/RandomGenerators.h"
#include "polymake/linalg.h"
#include "polymake/tropical/homogeneous_convex_hull.h"

namespace polymake { namespace tropical {

	/**
	  @brief Takes a matrix and returns the row indices where the first coordinate is nonzero and where the first coordinate is zero in  two different sets
	  @param Matrix<Rational> m The matrix whose rows we consider
	  @return std::pair<Set<int>, Set<int> > The first set contains the row indices of rows that start with a zero entry, the second set is the complement
	  */
   template <typename MType> inline
      std::pair<Set<int>, Set<int> > far_and_nonfar_vertices(const GenericMatrix<MType> &m) {
         const auto& first_col_supp = support(m.col(0));
         return std::pair<Set<int>, Set<int> >( sequence(0, m.rows()) - first_col_supp, first_col_supp);
      }

	/*
	 * @brief Takes a polyhedral complex and returns [[CONES]] summarized into one single incidence matrix.
	 * @param PolyhedralComplex
	 * @return IncidenceMatrix<>
	 */
	IncidenceMatrix<> all_cones_as_incidence(perl::Object complex);

	/*
	 * @brief Converts an incidence matrix to a Vector<Set<int> >
	 */
	template <typename MType> inline
      Vector<Set<int> > incMatrixToVector(const GenericIncidenceMatrix<MType> &i) {
         return Vector<Set<int> >(i.rows(), entire(rows(i)));
      }

   inline Array<Integer> randomInteger(const int& max_arg, const int &n) {
      static UniformlyRandomRanged<Integer> rg(max_arg);
      return Array<Integer>(n, rg.begin());
   }

	/**
	  @brief Computes all vectors of dimension n with entries +1 and -1. 
	  They are sorted such that each vector v has the row index determined by the sum: 
	  sum_{i: v_i = 1} 2^i (where i runs from 0 to n-1)
	  @param int n The column dimension of the matrix
	  @return Matrix<Rational> A 2^n by n matrix containing all +-1-vectors of dimension n
	  */
	Matrix<Rational> binaryMatrix(int n);

	/**
	  @brief Assumes v is a vector with entries +1 and -1 only. 
	  Returns sum_{i: v_i = 1} 2^i (where i runs from 0 to n-1
	  */
	template <typename VType> inline
      int binaryIndex(const GenericVector<VType> &v) {
         int result = 0;
         for(const auto &i : indices(attach_selector( v.top(), operations::positive())))
               result += pow(2,i);
         return result;
      }

	/**
   	@brief Helper function for the refinement function. 
		Given a polyhedral cell in terms of rays and lineality space, it computes, whether a given ray 
		is contained in this cell (possibly modulo (1,..,1)). 
	   @param Matrix<Rational> rays The rays of the cell
   	@param Matrix<Rational> lineality The lineality space of the cell
   	@param Vector<Rational> ray The ray to be tested
		@param bool is_projective Whether coordinates are given as tropical projective coordinates. 
		(False means they're affine).
		@param solver A convex hull solver
   	@returns true, if and only if ray lies in the cone
   */
	template <typename ch_solver>
		bool is_ray_in_cone(const Matrix<Rational> &rays, const Matrix<Rational> &lineality, 
				const Vector<Rational> &ray, bool is_projective, ch_solver& sv) {
			std::pair<Matrix<Rational>, Matrix<Rational> > facets = 
				is_projective ? enumerate_homogeneous_facets(rays, lineality, sv) :
				sv.enumerate_facets(rays,lineality, false,false);
			//Check equations
			for(auto l = entire(rows(facets.second)); !l.at_end(); ++l) {
				if(*l * ray != 0) return false;
			}
			//Check facets
			for(auto f = entire(rows(facets.first)); !f.at_end(); ++f) {
				if(*f * ray < 0) return false;
			}
			return true;
		}//END is_ray_in_cone
}}

#endif
