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
#include "polymake/atint/lll.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
	
      using namespace atintlog::donotlog;
    //using namespace atintlog::dolog;
    //using namespace atintlog::dotrace;
  
      
      ///////////////////////////////////////////////////////////////////////////////////////
      
	/*
		Wrapper for the GMP's gcdext method. 
	*/
	Integer gcdext(Integer a, Integer b, Integer &s, Integer &t) {
		//Copy values of a and b and initialize GMP variables for s and t and the gcd as well
		mpz_t at; mpz_t bt; mpz_t st; mpz_t tt; mpz_t gcd;
		mpz_init(at); mpz_init(bt); mpz_init(st); mpz_init(tt); mpz_init(gcd);
		
		mpz_set_si(at,a); mpz_set_si(bt,b); 
			
		//Compute gcdext
		mpz_gcdext(gcd,st,tt,at,bt);
				
		//Copy back values of s and t
		s = Integer(st);
		t = Integer(tt);
		
		return Integer(gcd);	
	}	
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	/*
		Computes a determinant +-1 matrix U with Z entries, such that AU = (0 | T), where T is a regular lower
		triangular matrix (and A = matrix). k will be set to the dimension of the kernel of A (= the number of 0 columns of
		AU), so that the first k columns of U are a Z-basis of Ker(A). This algorithm is based on the 
		algorithm found in \"Cohen: A course in computational algebraic number theory, p. 74\"
		This algorithm has the problem that the coefficients created by the gcdext algorithm explode quickly, even
		if the initial entries of matrix are small. This makes it useless for practical applications.
	*/
	
	/*Matrix<Integer> znormaltransform(const Matrix<Integer>& matrix, Integer &kdim) {
		//Iteration variables		
		int m = matrix.rows();
		int n = matrix.cols();
		int i = m;
		int j = n;
		int k = n;		
		int l = (m <= n ? 1 : m-n+1);
		Matrix<Integer> A(matrix); 		
		
		//Helper variables
		Integer u;
		Integer v;
		Integer g;		
		Vector<Integer> B;	
		
		//Result variable
		Matrix<Integer> U = unit_matrix<Integer>(n);
		
		//Algorithm to compute kernel
		//cout << "Reducing A = \n" << A << endl;
		while(true) {
			while(j > 1) {
				j--;
				//cout << "i = " << i << ",j = " << j << ", aij = " << A(i-1,j-1) << ", aik = " << A(i-1,k-1) << endl;
				if(A(i-1,j-1) != 0) {
					 
					 g = gcdext(A(i-1,k-1),A(i-1,j-1),u,v);
					 
					 
					 //cout << "gcd(" << A(i-1,k-1) << "," << A(i-1,j-1) << ") = " << g << "; u = " << u << ", v = " << v << endl;
					 B = u * A.col(k-1) + v * A.col(j-1);
					 
					 Integer aikbyg = (div_exact(A(i-1,k-1),g));
					 Integer aijbyg = (div_exact(A(i-1,j-1),g));					 
					 
					 A.col(j-1) = aikbyg * A.col(j-1)
					 				- aijbyg * A.col(k-1);
					 A.col(k-1) = B;
					 B = u * U.col(k-1) + v * U.col(j-1);
					 U.col(j-1) = aikbyg * U.col(j-1)
					 				- aijbyg * U.col(k-1);
					 U.col(k-1) = B;
				} 
			}		
			if(A(i-1,k-1) == 0) {
				k++;			
			}
			if(i != l) {
				i--; k--;
				j = k;
			}
			else {	
				kdim = Integer(k-1);
				return U;
			}			
		}
	}*/

	///////////////////////////////////////////////////////////////////////////////////////
	
	/*
		Takes a rational nxm matrix and multiplies each row with a minimal integer such that it becomes 
		a primitive vector in Z^m 
	*/
	Matrix<Integer> makePrimitiveInteger(const Matrix<Rational> &m) {
		Matrix<Integer> result(m.rows(), m.cols());
		for(int r = 0; r < m.rows(); r++) { 
			Integer lc = 1;
			for(int c = 0; c < m.cols(); c++) {
				lc = lcm(lc,denominator(m(r,c)));
			}
			result.row(r) = lc* m.row(r);
		}
		return result;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	Vector<Integer> makePrimitiveInteger(const Vector<Rational> &v) {
	  Vector<Integer> result(v.dim());
	  Integer lc = 1;
	  for(int c = 0; c < v.dim(); c++) {
		lc = lcm(lc,denominator(v[c]));
	  }
	  result = lc* v;
	  
	  return result;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	/*
		Assuming that tauamtrix is a codimension one subspace of sigmamatrix (considered as dual spaces of the row spaces), computes a lattice normal vector
		The orientation of the normal vector is determined by additionalRay, in the sense that it will point "towards" this ray.
		More precisely, if h is the hypersurface defining tau wrt sigma and h * (additionalRay) >= 0, then h* normal >= 0
	*/
	Vector<Integer> latticeNormal(const Matrix<Rational> &tmatrix, const Matrix<Rational> &smatrix, const Vector<Rational> &additionalRay) {
	      dbgtrace << "Making the matrices integer" << endl;
	      Matrix<Integer> taumatrix = makePrimitiveInteger(tmatrix);
	      Matrix<Integer> sigmamatrix = makePrimitiveInteger(smatrix);
	      dbgtrace << "Taumatrix = \n" << taumatrix << "\nSigmamatrix = \n" << sigmamatrix << endl;
	      int rk = rank(sigmamatrix);
	      int rowcount = taumatrix.rows();
	      //Find the row of taumatrix that is not in the span of sigmamatrix
	      //and append it to the bottom of sigmamatrix		
	      for(int r = 0; r < rowcount; r++) {
		      Vector<Integer> testedRow = taumatrix.row(r);
		      if(rank((sigmamatrix / testedRow)) > rk) {
			      sigmamatrix = (sigmamatrix / testedRow);
			      break;				
		      }
		      //If we arrive at the end of the loop, we haven't found anything and an error
		      //is thrown			
		      if(r == taumatrix.rows()-1) {
			      throw std::runtime_error("latticeNormal: tau is not a face of sigma");
		      }
	      }
	      dbgtrace << "Transformed sigmamatrix = \n" << sigmamatrix << endl;
	      Integer k;
	      //Matrix<Integer> tfmatrix = znormaltransform(sigmamatrix, k); // --> Dies on large integers
	      Matrix<Integer> tfmatrix;
	      lllHNF( T(sigmamatrix), tfmatrix,k);//Compute the HNF of sigmamatrix-transposed
	      //Now the normal vector is the (n-k)-th row of tfmatrix (with n = ambient dimension)
	      dbgtrace << "The transformation matrix is \n" << tfmatrix << "\nKernel dimension is " << k << endl;
	      Vector<Integer> lnormal = tfmatrix.row	(sigmamatrix.cols()- 1 - k);
	      
	      
	      //Determine orientation
	      dbgtrace << "Determining orientation" << endl;
	      Vector<Integer> hyper = sigmamatrix.row(sigmamatrix.rows()-1);
	      dbgtrace << "Defining hyperplane = " << hyper << endl;
	      dbgtrace << "Normal vector = " << lnormal << endl;
	      dbgtrace << "Additional ray = " << additionalRay << endl;
	      hyper = (hyper * additionalRay >= 0 ? hyper : (-1)*hyper);
	      dbgtrace << "Setting hyperplane orientation to " << hyper << endl;
	      return (hyper * lnormal >= 0? lnormal : (-1) * lnormal);	      
	}
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	/*
		Assuming that tau is a codimension one face of sigma, computes a lattice normal vector
	*/	
	Vector<Integer> latticeNormalByCone(const perl::Object &tau, const perl::Object &sigma) {
		//Read out the matrices whose rows are the dual basis to V_tau (or V_sigma) and make them integer 		
		Matrix<Rational> taumatrix = tau.give("LINEAR_SPAN");
		Matrix<Rational> sigmamatrix = sigma.give("LINEAR_SPAN");
		dbgtrace << "Linspan of tau = \n" << taumatrix << endl;
		dbgtrace << "Linspan of sigma = \n" << sigmamatrix << endl;
		
		//Find the additional ray of sigma
		Matrix<Rational> sigmarays =sigma.give("RAYS");
		
		//We have to check manually, whether tau is just the origin, since the give("RAYS") would
		//break in this case (bug?)
		Vector<Rational> additionalRay;
		if(tau.give("CONE_DIM") == 0) {
		  additionalRay = sigmarays.row(0);
		}
		else {
		   Matrix<Rational> taurays = tau.give("RAYS");
		   int rk = rank(taurays);
		    for(int row = 0; row < sigmarays.rows(); row++) {
			if(rank((taurays / sigmarays.row(row))) > rk) {
			    additionalRay = sigmarays.row(row);
			    break;
			}
		    }
		}
				
		dbgtrace << "additional ray = " << additionalRay << endl;
		return latticeNormal(taumatrix, sigmamatrix,additionalRay);		
	}
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	/*
	  Compute a lattice basis of the vector space spanned by cone (more precisely it computes Lambda_cone)
	*/
	Matrix<Integer> latticeBasis(const perl::Object &cone) {
	  Matrix<Rational> prematrix = cone.give("LINEAR_SPAN");
	  Matrix<Integer> conematrix = makePrimitiveInteger(prematrix);
	  conematrix = T(conematrix);
	  Matrix<Integer> tfmatrix;
	  Integer k;
	  lllHNF(conematrix, tfmatrix,k);
	  //Copy the last k rows of tfmatrix 
	  Matrix<Integer> latticeB(k,conematrix.rows());
	  for(int r = 0; r < latticeB.rows(); r++) {
	      latticeB.row(r) = tfmatrix.row(tfmatrix.rows() - 1 - r);
	  }
	  return latticeB;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////
	
	Vector<Rational> linearRepresentation(const Vector<Rational> &v, const Matrix<Rational> &generators) {
	  Vector<Rational> solution(generators.rows());
	  //Copy arguments
	  Matrix<Rational> A(generators);
	  Vector<Rational> w(v);
	  Matrix<Rational> U = unit_matrix<Rational>(generators.rows()); //The transformation matrix for A
	  
	  if(v.dim() != generators.cols()) {
	    throw std::runtime_error("Dimension mismatch of generating set and vector");
	  }
	  
	  int usedRows = 0; //The number of times we actually used a row for reducing
	  Rational coeff = 0;
	  Vector<Rational> zv = zero_vector<Rational>(w.dim());
	  if(w == zv) return solution;
	  
	  //Go through each column of w / A and try to reduce it using a row of A
	  for(int c = 0; c < v.dim(); c++) {
	    dbgtrace << "Reducing column " << c+1 << endl;
	    //Find the first row of A such that A(row,c) != 0 and use this row to reduce the column c.
	    //Then move it to the end (i.e. above the row we used the last time)
	    for(int r = 0; r < A.rows() - usedRows; r++) {
	      if(A(r,c) != 0) {
		dbgtrace << "Reducing with row " << r+1 << endl;
		//First reduce w, if necessary
		if(w[c] != 0) {
		  coeff = w[c] / A(r,c);
		  solution += coeff * U.row(r);
		  w -= coeff * A.row(r);
		  if(w == zv) return solution;
		}
		//Actually we first move the row to the end
		for(int sc = 0; sc < A.cols(); sc++) {
		    A(r,sc).swap(A(A.rows()-usedRows-1,sc));
		}
		for(int uc = 0; uc < U.cols(); uc++) {
		    U(r,uc).swap(U(A.rows()-usedRows-1,uc));
		}
		//Now we reduce all rows below it (the ones above, we don't need anymore)
		for(int s = 0; s < A.rows() - usedRows - 1; s++) {
		  if(A(s,c) != 0) {
		    coeff = A(s,c) / A(A.rows() - usedRows -1,c);
		    A.row(s) -= coeff * A.row(A.rows() - usedRows - 1);
		    U.row(s) -= coeff * U.row(A.rows() - usedRows - 1);
		  }
		}
		usedRows++;
		dbgtrace << "Now A = \n" << A << endl;
		dbgtrace << "Now w = \n" << w << endl;
		break;
	      }
	      //If we arrive at this point, we cant reduce w, so its not in the linear span
	      if(r == A.rows() - usedRows - 1 && w[c] != 0) {
		dbgtrace << "Not in linear span" << endl;
		return Vector<Rational>(0);
	      }
	    }
	  }
	  
	  return solution;
	}
	
// ------------------------- PERL WRAPPERS ---------------------------------------------------
	
/*UserFunction4perl("# @category Linear algebra"
                  "# Computes a matrix //U// with entries in //Z// with determinatn +-1, such that"
                  "# //A// times //U// = (0 | T), with T in lower triangular form. It sets //kdim// to"
                  "# be the dimension of Ker(A), such that the first k columns of U are a Z-basis of"
                  "# Ker(A)"
                  "# @param Matrix matrix the matrix for which the transformation is computed"
                  "# @param kdim This will be set to dim Ker(A)"
                  "# @return Matrix",
                  &znormaltransform,"znormaltransform(Matrix<Integer>,Integer)");   */               

UserFunction4perl("# @category Arithmetic"
					  "# Computes the gcd of //a// and //b// and returns it. //s// and //t// are set"
					  "# such that gcd(//a//,//b//) = //s// * //a// + //t// * //b//"
					  "# @param Integer a first argument of gcd(,)"
					  "# @param Integer b second argument of gcd(,)"
					  "# @param Integer s coefficient of a (will be set)"
					  "# @param Integer t coefficient of b (will be set)"
					  "# @return Integer",
					  &gcdext, "gcdext(Integer, Integer, Integer, Integer)");               

UserFunction4perl("# @category Tropical geometry"
						"# Assuming that tau is a codimension one face of sigma, computes a representative of"
						"# the primitive lattice normal vector of sigma with respect to tau"
						"# @param polytope::Cone tau a codimension one face of tau, given as a cone"
						"# @param polytope::Cone sigma an arbitrary cone"
						"# @return Vector",
						&latticeNormalByCone, "latticeNormalByCone(polytope::Cone,polytope::Cone)");
					

UserFunction4perl("# @category Tropical geometry"
						"# Takes two matrices whose rows define the dual of the linear span of cone tau and sigma."
						"# Assuming that tau is a codimension one face of sigma, computes a representative of"
						"# the primitive lattice normal vector of sigma with respect to tau"
						"# @param Matrix taumatrix a codimension one face of sigma, given as a matrix defining its linear span"
						"# @param Matrix sigmamatrix an arbitrary cone, given as a matrix defining its linear span"
						"# @param Vector additionalRay A ray that is contained in sigma, but not in tau. Used to calculate proper orientation of the normal vector."
						"# @return Vector",
						&latticeNormal, "latticeNormal(Matrix<Rational>,Matrix<Rational>, Vector<Rational>)");

UserFunction4perl("# @category Tropical geometry"
					  "# Takes a cone and computes a Z-basis of the vector space spanned by the cone,"
					  "# returned as row vectors of a matrix"
					  "# @param polytope::Cone cone A cone for which a Z-basis is to be computed"
					  "# @return Matrix<Integer>",
					  &latticeBasis,"latticeBasis(polytope::Cone)");
					  
UserFunction4perl("# @category Linear algebra"
		  "# Takes a vector v and a matrix with column dimension equal to the dimension of v. Assuming that "
		  "# v is in the row span of the matrix, it computes one(!) possible representation of v in these "
		  "# generators"
		  "# @param Vector v The vector supposed to be contained in the row span of the generators"
		  "# @param Matrix generators  A set of row vectors whose linear span should contain v"
		  "# @return Vector A vector (a1,..,an) such that v = (a1,...,an) * generators. It returns a vector of" 
		  "# dimension 0, if v is not in the span of the generators. An error is thrown if the dimensions of "
		  "# v and the generators mismatch",
		  &linearRepresentation,"linearRepresentation(Vector<Rational>,Matrix<Rational>)");
					  
}
}
