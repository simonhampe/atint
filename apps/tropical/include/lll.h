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
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"

#ifndef ATINT_LLL_H
#define ATINT_LLL_H

namespace polymake { namespace tropical {
  
  /**
    @brief Computes a row Hermite normal form (HNF) of matrix. It returns the normal form and stores the corresponding unimodular transformation matrix in tfMatrix and the dimension of the kernel of the transposed matrix in kernelDimension. This algorithm is the LLL based HNF algorithm described in "Extended gcd and Hermite normal form algorithms via lattice base reduction" by Havas, Majewski,Matthews. That algorithm actually computes the HNF in row-reversed order, so we reverse it by hand before returning.
    @param matrix The matrix for which the normal form should be computed
    @param tfmatrix The unimodular transformation matrix is stored in this, i.e. lllHNF = tfmatrix * matrix
    @param kernelDimension The dimension of the kernel of transpose(matrix), or matrix.cols() - rank(matrix) is stored in this
    @return A row Hermite normal form of matrix. If matrix.cols() <= matrix.rows(), this means that lllHNF = T / 0_(kernelDimension x matrix.cols), where T is a regular, upper triangular square matrix. Hence the last k rows of tfmatrix form a Z-basis of Ker(transpose(matrix))
  */
  Matrix<Integer> lllHNF(const Matrix<Integer> &matrix, Matrix<Integer> &tfmatrix, Integer &kernelDimension);
  
  
}}

#endif