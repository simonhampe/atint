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
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint {

    using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
  
    ///////////////////////////////////////////////////////////////////////////////////////
    
    Matrix<Integer> lllHNF(const Matrix<Integer> &matrix, Matrix<Integer> &tfmatrix, Integer &kernelDimension) {
	
	
      
	//Initialize variables
	int rows = matrix.rows();
	int cols = matrix.cols();
	Matrix<Integer> A(matrix);
	tfmatrix = unit_matrix<Integer>(rows);
	
	Vector<Integer> D(rows+1); //D(i) is the product of the norms squared of the first i Gram (row) vectors of 
				   //tfmatrix, initially all 1
				   //D(0) = 1 throughout
	  for(int i = 0; i <= rows; i++) D[i] = 1;
	Matrix<Integer> Lambda(rows,rows); //The LLL - Lambda-coefficients, intially all 0
	  
	int m1 = 3; int n1 = 4; //The efficiency coefficient of the LLL algorithm, alpha = m1/n1
				//can be in (1/4,1]. Higher = Slower, but smaller results
	int k = 1; //The current row to be reduced
	
	
	
	//If there is exactly one nonzero element in the first nonzero column of A and it lies in the last row
	//and is negative, then we multiply the last row by -1
	//This is a special case not covered by the algorithm below. If we omit it, then A(0,0) will be negative
	Vector<Integer> zvector(rows);
	for(int c = 0; c < cols; c++) {
	    if(A.col(c) != zvector) {
	      for(int r = 0; r < rows; r++) {
		  if(A(r,c) != 0) {
		      if(r == rows-1 && A(r,c) < 0) {
			A.row(r) = - A.row(r);
			tfmatrix(r,r) = - 1;
		      }
		      break;
		  }
	      }
	      break;
	    }
	}
	
	
	dbglog << "Reducing a " << rows << " x " << cols << " matrix" << endl;
	while(k < rows) {
	    dbglog << "Starting reduction of row k = " << k << endl;
	    for(int i = k-1; i >= 0; i--) {
		dbglog << "Reducing wrt i = " << i << endl;
		//Reduce row k, using row i --------------------------------------------------------------------------------
		//First find the first columns col1,col2, such that ai,col1 != 0, ak,col2 != 0 (if none is found, set to n+1)
		int col1 = cols;
		int col2 = cols;
		for(int j = 0; j < cols; j++) {
		    if(col1 == cols && A(i,j) != 0) col1 = j;
		    if(col2 == cols && A(k,j) != 0) col2 = j;
		    if(col1 < cols && col2 < cols) break;
		}
		dbglog << "col1 = " << col1 << ", col2 = " << col2 << endl;
		//Check for positivity of this first nonzero element. If it is negative, invert the appropriate 
		//LLL coefficients Lambda and the rows of the transformation matrix
		if(col1 < cols && A(i,col1) < 0) {
		    A.row(i) = - A.row(i);
		    tfmatrix.row(i) = - tfmatrix.row(i);
		    for(int r = 1; r < rows; r++) { 
			for(int s = 0; s < r; s++) {
			    if( (r == i) || (s == i) ) Lambda(r,s) = - Lambda(r,s);
			}
		    }
		}
		//Now find the appropriate reduction multiplier q
		Integer q(0);
		if(col1 < cols) {
		    q = A(k,col1) / A(i,col1);
		}
		else {
		   Integer plambda = Lambda(k,i);
		      plambda = plambda > 0? plambda : - plambda;
		   if(2* plambda > D[i+1]) {
		      dbglog << "Reducing only B" << endl;
		      dbglog << "L_k,i = " << Lambda(k,i) << ", Di+1 = " << D[i+1] << endl;
		      q = (2 * plambda + D[i+1]) / (2*D[i+1]);// = nearest integer to |Lambda(k,i)|/ D[i+1];
		      if(Lambda(k,i) < 0) q = -q;
		   }
		}
		dbglog << "Setting q = " << q << endl;
		//If q != 0, apply reduction step
		if(q != 0) {
		    A.row(k) = A.row(k) - q* A.row(i);
		    tfmatrix.row(k) = tfmatrix.row(k) - q* tfmatrix.row(i);
		    Lambda(k,i) = Lambda(k,i) - q* D[i+1];
		    for(int j = 0; j <= i-1; j++) {
			Lambda(k,j) = Lambda(k,j) - q* Lambda(i,j);
		    }
		}
		dbglog << "Now A =\n" << A << endl;
		dbglog << "B = \n" << tfmatrix << endl;
		//if i=k-1, check for the Lovasz condition ---------------------------------------------------------------
		//If it is fulfilled, reduce wrt the other rows, otherwise step down one row
		if(i == k-1) {
		    if( ((col1 <= col2) && (col1 < cols)) ||
			((col1 == col2) && (col1 == cols) && (n1* (D[k-1]*D[k+1]) + (Lambda(k,k-1)*Lambda(k,k-1))) < m1* D[k]* D[k])) {
			dbglog << "Lovasz condition not fulfilled, swapping " << k << " with " << k-1 << endl;
			//If it is not fulfilled, swap rows k and k-1
			for(int c = 0; c < cols; c++) A(k,c).swap(A(k-1,c));
			for(int c = 0; c < rows; c++) tfmatrix(k,c).swap(tfmatrix(k-1,c));			
			for( int j = 0; j <= k-2; j++) {
			    Lambda(k,j).swap(Lambda(k-1,j));
			}
			dbglog << "Recalculating D's and Lambdas" << endl;
			for( int i = k+1; i < rows; i++) {
			    Integer t = Lambda(i,k-1)* D[k+1] - Lambda(i,k-1) * Lambda(k,k-1);
			    dbglog << "Setting t = " << t << endl;
			    dbglog << "Li,k-1 * Lk,k-1  = " << (Lambda(i,k-1)* Lambda(k,k-1)) << endl;
			    dbglog << "Li,k * Dk-1 = " << (Lambda(i,k) * D[k-1]) << endl;
			    Lambda(i,k-1) = ( (Lambda(i,k-1) * Lambda(k,k-1)) + (Lambda(i,k) * D[k-1])  ) / D[k];
			    dbglog << "Setting Lambda_i,k-1 = " << Lambda(i,k-1) << endl;
			    Lambda(i,k) = t / D[k];
			    dbglog << "Setting Lambda_i,k = " << Lambda(i,k) << endl;
			}
			D[k] = (  (D[k-1] * D[k+1]) + (Lambda(k,k-1) * Lambda(k,k-1))  ) / D[k]; 
			dbglog << "Setting Dk-1 = " << D[k] << endl;
			//After swapping, step down one row
			if(k > 1) k--;
			dbglog << "Stepping down to k = " << k << endl;
			break;
		    }
		}
		//After row k has been reduced wrt all lower rows, go up one row -----------------------------------------
		if(i == 0) {
		  k++;
		  dbglog << "Stepping up to k = " << k << endl;
		}
		
		
	    }
	}
	//Return results--------------------------------
	
	//Invert row order of matrix A (and tfmatrix)
	for(int swapr = 0; 2*swapr < rows; swapr++) {
	    for(int swapc = 0; swapc < cols; swapc++) A(swapr,swapc).swap(A(rows-swapr-1,swapc));
 	    for(int swapc = 0; swapc < rows; swapc++) tfmatrix(swapr,swapc).swap(tfmatrix(rows-swapr-1,swapc));
	}
	
	//Find first nonzero row, searching from below
	int fnzero = -1;
	for(int r = rows-1; r >= 0; r--) {
	    if(A(r,cols-1) != 0) {
		fnzero = r; break;
	    }
	}
	kernelDimension = fnzero == -1? 0 : (rows -1 ) - fnzero;
	return A;
    }
    
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
    
 UserFunction4perl("# @category Linear algebra"
                  "# Computes a row HNF of an integer matrix A"
                  "# It returns the normal form and stores the unimodular transformation matrix and the kernel dimension  of the transposed matrix in"
                  "# the last two parameters. The algorithm is the LLL-based HNF alg. by Havas, Majevski, Matthews"
                  "# @param Matrix<Integer> matrix the matrix for which the transformation is computed"
                  "# @param Matrix<Integer> tfmatrix The matrix that will contain the transformation matrix"
                  "# @param Integer kdim This will be set to dim Ker(A)"                  
                  "# @return Matrix",
                  &lllHNF,"lllHNF(Matrix<Integer>,Matrix<Integer>,Integer)");           

  
}
}