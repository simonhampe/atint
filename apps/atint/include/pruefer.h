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
 
 Contains the signatures for pruefer.cc
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"


#ifndef PRUEFER_H
#define PRUEFER_H

/**
 @brief This computes the set of all Pruefer sequences of order n fulfilling one of a list of certain valency condition. These conditions are given as a matrix of integers, seen as a list of row vectors. Each row has length k+1, where k is the number of bounded edges. Column c_i, i= 0,..,k stands for the interior vector labelled n+i and indicates what valence it should have. That means that in a sequence corresponding to row r the vertex n+i occurs valences(r,i)-1 times.
 @param int n The number of leaves of rational curves for which we compute Pruefer sequences.
 @param Matrix<int> valences. Each row prescribes a valence for each interior vertex
 @return Matrix<int> A list of all Pruefer sequences fulfilling one of the valency conditions (as row vectors).
 */
Matrix<int> prueferSequenceFromValences(int n, Matrix<int> valences);

/**
 @brief This computes the set of all Pruefer sequences corresponding to k-dimensional combinatorial types in M_0,n
 @param int n The number of leaves of rational curves
 @param int k The number of  bounded edges in rational curves
 @return Matrix<int> A list of all Pruefer sequences of combinatorial types of curves with n leaves and k bounded edges (as row vectors).
 */
Matrix<int> dimension_k_prueferSequence(int n, int k) ;

/**
 @brief Takes a list of Pruefer sequences and decodes them into a WeightedComplex containing for each sequence the cone that corresponds to it
 @param int n The parameter n of the M_0,n on which the sequence is defined.
 @param Matrix<int> A list of Pruefer sequences (as row vectors)
 @return perl::Object A WeightedComplex (but without any weights)
 */
perl::Object complex_from_prueferSequences(int n, Matrix<int> pseq);

#endif // PRUEFER_H