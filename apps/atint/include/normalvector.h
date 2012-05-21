/*
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License version 2 as published by the Free Software Foundation.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.

You should have received a copy of the GNU Library General Public License
along with this library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
Boston, MA 02110-1301, USA.

---
Copyright (C) 2011, Simon Hampe <hampe@mathematik.uni-kl.de>
*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"

#ifndef ATINT_NORMALVECTOR_H
#define ATINT_NORMALVECTOR_H

namespace polymake { namespace atint {

/**
@brief Computes the gcd of a and b and returns it. s and t are set such that gcd(a,b) = s * a + t * b
 @param Integer a first argument of gcd(,)
 @param Integer b second argument of gcd(,)
 @param Integer s coefficient of a (will be set)
 @param Integer t coefficient of b (will be set)
 @return The gcd of a and b",
*/
Integer gcdext(Integer a, Integer b, Integer &s, Integer &t);


/**
 @brief Takes a rational vector and makes it an Integer vector by multiplying with the least common multiple of the denominators of the entries.
 */
Vector<Integer> makeInteger(const Vector<Rational> &v);

/**
 @brief Takes a rational matrix and makes it an Integer vector by multiplying each with the least common multiple of the denominators of the entries.
 */
Matrix<Integer> makeInteger(const Matrix<Rational> &v);

/**
@brief Takes a rational matrix and makes each row a primitive vector in Z^n by multiplying each row with the least common multiple of the denominators of the entries and dividing afterwards by the gcd of the entries.
*/
Matrix<Integer> makePrimitiveInteger(const Matrix<Rational> &m);

/**
@brief Takes a rational vector and transforms it into a primitive vector in Z^n by multiplying with the least common multiple of the denominators of the entries and dividing afterwards by the gcd of the entries.
*/
Vector<Integer> makePrimitiveInteger(const Vector<Rational> &v);

/**
@brief Takes two matrices whose rows define the dual of the linear span of cone tau and sigma. Assuming that tau is a codimension one face of sigma, computes a representative of the primitive lattice normal vector of sigma wrt tau. 
@param Matrix taumatrix a codimension one face of sigma, given as a matrix defining its linear span
@param Matrix sigmamatrix an arbitrary cone, given as a matrix defining its linear span
@param Vector additionalRay A ray that is contained in sigma, but not in tau. Used to calculate proper orientation of the normal vector. The orientation of the normal vector is determined by additionalRay, in the sense that it will point "towards" this ray. More precisely, if h is the hypersurface defining tau wrt sigma and h * (additionalRay) >= 0, then h* normal >= 0
@return The lattice normal of sigma with respect to tau.
*/
Vector<Integer> latticeNormal(const Matrix<Rational> &tmatrix, const Matrix<Rational> &smatrix, const Vector<Rational> &additionalRay);

/**
@brief Assuming that tau is a codimension one face of sigma computes a representative of the primitive lattice normal vector of sigma with respect to tau.
@param polytope::Cone tau a codimension one face of tau, given as a cone
@param polytope:: Cone sigma an arbitrary cone
@return The lattice normal vector of sigma wrt tau
*/
Vector<Integer> latticeNormalByCone(const perl::Object &tau, const perl::Object &sigma);

/**
@brief Computes a lattice basis for a cone spanned by a list of rays and with a given lineality space
@param Matrix<Rational> rays The rays of the cone
@param Matrix<Rational> linspace The lineality space generators of the cone
@param bool uses_homog Optional parameter, false by default. If true, the cone is assumed to be given in homogeneous coordinates and the lattice basis will be computed in such a way that the first, homogenizing coordinate is zero.
@returns Matrix<Integer> A basis for the lattice of the cone
*/
Matrix<Integer> latticeBasisFromRays(const Matrix<Rational> &rays, const Matrix<Rational> &linspace, bool uses_homog = false);

/**
@brief Takes a cone and computes a Z-basis of the vector space spanned by the cone, returned as row vectors of a matrix.
@param polytope::Cone cone A cone for which a Z-basis is to be computed
@return An integer matrix whose rows span the linear space spanned by cone
*/
Matrix<Integer> latticeBasis(const perl::Object &cone);

/**
  @brief Takes a vector v and a matrix with column dimension equal to the dimension of v. Assuming that v is in the row span of the matrix, it computes one(!) possible representation of v in these generators. It does this by performing a standard (partial) gaussian reduction algorithm
  @param v The vector supposed to be contained in the row span of the generators
  @param generators  A set of row vectors whose linear span should contain v
  @return A vector (a1,..,an) such that v = (a1,...,an) * generators. It returns a vector of dimension 0, if
  v is not in the span of the generators. An error is thrown if the dimensions of v and the generators mismatch
  
*/
Vector<Rational> linearRepresentation(const Vector<Rational> &v, const Matrix<Rational> &generators);

/**
@brief  This method takes a set of row indices for [[RAYS]] (or [[CMPLX_RAYS]] in the homogeneous case) and a vector that is supposed to be in the span of these row vectors and the lineality space (or their affine space, see description of [[LATTICE_NORMAL_FCT_VECTOR]]). It then computes the corresponding representation in these vectors
@param s a set of row indices of [[RAYS]] (or [[CMPLX_RAYS]])
@param v a vector supposed to lie in the span (or affine span) of [[RAYS]] (or [[CMPLX_RAYS]]) + [[LINEALITY_SPACE]]
@param ambient_dim The ambient dimension of the fan
@param uses_homog Whether the fan uses homogeneous coordinates
@param rays The matrix of [[RAYS]] or [[CMPLX_RAYS]]
@param linealitySpace A matrix of generators of the lineality space
@param lineality_dim The dimension of the lineality space
@return A vector of length [[N_RAYS]] + [[LINEALITY_DIM]] (or [[CMPLX_RAYS]]->rows() + [[LINEALITY_DIM]] in the homogeneous case) with linear coefficients of a representation in the generators chosen via s. The last elements always refer to the lineality space.
@throw std::runtime_error If the vector is not in the linear span of the give vectors
*/
Vector<Rational> functionRepresentationVector(const Set<int> &rayIndices, const Vector<Rational> &v,
					      int ambient_dim, bool uses_homog, 
					      const Matrix<Rational> &rays,
					      const Matrix<Rational> &linealitySpace,
					      int lineality_dim); 



}}

#endif // __h_