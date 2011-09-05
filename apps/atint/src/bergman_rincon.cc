/*
 T his *program is free software; you can redistribute it and/or
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
 
 This file contains the algorithms to compute bergman fans as described
 by Felipe Rincón in his paper "Computing tropical linear spaces".
 See also http://math.berkeley.edu/~felipe/bergman
 */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Map.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  //using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
  
  /**
   @brief Takes a matrix and computes all sets of column indices, s.t. the corresponding column vectors form a basis of the column space of the matrix.
   @param Matrix<Rational> m 
   @return IncidenceMatrix An incidence matrix whose rows correspond to the bases and where each row contains the colum indices contained in that basis.
   */
  IncidenceMatrix<> computeMatrixBases(Matrix<Rational> m) {
    //First we determine the rank of the matrix
    int r = rank(m);
    
    Vector<Set<int> > result;
    Array<Set<int> > rsets = pm::Subsets_of_k<Set<int> > ( sequence(0,m.cols()),r );
    for(int i = 0; i < rsets.size(); i++) {
      if(rank(m.minor(All,rsets[i])) == r) {
	result |= rsets[i];
      }
    }
    return IncidenceMatrix<>(result);
  }
  
  /**
   @brief Takes a matrix and finds all columns that are coloops in the matrix matroid, i.e. all column vectors that are contained in every basis of the column space. In other words, the rank of the matrix decreases, if I remove such a column
   @param Matrix<Rational> m
   @return Set<int> A set of column indices of columns that are coloops
   */
  Set<int> computeMatrixColoops(Matrix<Rational> m) {
    //First we determine the rank of the matrix
    int r = rank(m);
    
    Set<int> coloops;
    for(int c = 0; c < m.cols(); c++) {
      if(rank(m.minor(All,~scalar2set(c))) < r) {
	coloops += c;
      }
    }
    return coloops;
  }
  
  /**
    @brief For a given base B of a matroid and an element not in B, this method computes the set C(k,B)-{k}, where C(k,B) is the fundamental circuit of k over B
    @param IncidenceMatrix<> bases The set of all bases of the matroid
    @param int mybase The row index of the relevant base in the matrix bases
    @param int k An element not in mybase
    @return Set<int> The set C(k,B) - {k}
  */
  inline Set<int> computeFk(const IncidenceMatrix<> &bases, int mybase, int k) {
    Set<int> Fk;
    //We go through all elements i of mybase and check if B - i + k is independent, i.e. equal to a basis
    Set<int> B = bases.row(mybase);
    for(Entire<Set<int> >::iterator i = entire(B); !i.at_end(); i++) {
      Set<int> C = B; C -= (*i); C += k;
      bool independent = false;
      //Check if the set C is contained in any basis
      for(int row = 0; row < bases.rows(); row++) {
	if(row != mybase) {
	  if(C == bases.row(row)) {
	    independent = true;
	    break;
	  }
	}
      }
      if(independent) {
	Fk += (*i);
      }
    }
    return Fk;
  }
  
  /**
    @brief For a given base B of a linear matroid and an element not in B, this method computes the set C(k,B)-{k}, where C(k,B) is the fundamental circuit of k over B. This is much faster then computeFk for general matroids
    @param IncidenceMatrix<> bases The set of all bases of the matroid
    @param int mybase The row index of the relevant base in the matrix bases
    @param int k An element not in mybase
    @param Matrix<Rational> m The matrix representing the matroid
    @return Set<int> The set C(k,B) - {k}
   */
  inline Set<int> computeFkLinear(const IncidenceMatrix<> &bases, int mybase, int k, const Matrix<Rational> &m) {
    Set<int> Fk;
    Vector<int> B(bases.row(mybase));
    //Compute the row reduced matrix s.t. the minor corresponding to mybase is the identity
    Matrix<Rational> I = inv(m.minor(All,bases.row(mybase)));
    Matrix<Rational> A = I * m;
    for(int i = 0; i < A.rows(); i++) {
      if(A(i,k) != 0) {
	Fk += B[i];
      }
    }
    return Fk;
  }
  
  /**
    @brief Computes the bergman fan of a matroid that is supposed to be loop- and coloop-free
    @param int n The number of elements of the ground set of the matroid
    @param IncidenceMatrix<> bases The bases of the matroid
    @param bool is_linear Whether the matroid is represented by a Q-matrix
    @param Matrix<Rational> m if is_linear is true, this contains the matrix representing the matroid
    @return WeightedComplex The bergman fan of the matroid, including the (1,..,1)-lineality space.
  */
  perl::Object bergman_fan(int n, const IncidenceMatrix<> &bases, bool is_linear, const Matrix<Rational> &m) {
      //Prepare result variables
      Matrix<Rational> rays(0,n);
      Vector<Set<int> > cones;
      Matrix<Rational> lineality(0,n);
	lineality /= ones_vector<Rational>(n);
      Vector<Integer> weights;
      
      //Now compute cones for each basis
      Set<int> complete = sequence(0,n);
      int count = 0;
      for(int B = 0; B < bases.rows(); B++) {
	//Initialize fundamental circuits
	Vector<int> complement(complete - bases.row(B));
	Vector<Set<int> > Fksets (complement.dim());
	Vector<Vector<int > > Fklists(complement.dim());
	dbgtrace << "Computing for basis " << bases.row(B) << endl;
	for(int k = 0; k < Fksets.dim(); k++) {
	  Fksets[k] = is_linear? computeFkLinear(bases,B,complement[k],m) :
				 computeFk(bases,B,complement[k]);
	  Fklists[k] = Vector<int>(Fksets[k]);
	  dbgtrace << "For k = " << complement[k] << " Fk = " << Fklists[k] << endl;
	}
	
	//Now we go through all the regressive compatible pairs (p,L)
	Vector<int> p(complement.dim()); //contains the images p(k)
	Vector<int> L; //the ordering on B (first element = smallest, etc)
	int k = 0; //The k we currently define an image for (actually: index in complement)
	//For given k, iterations[0]..iterations[k] describe the iterations for all k' <= k
	// Semantics: (-1,-1): We just started with this k, check F_k cap L
	// (0,-1): Start adding elements from F_k \ L
	// (b,x): In the last iteration we took the element at index b in F_k and inserted it
	// before position x. Find the next valid position
	Vector<std::pair<int,int> > iterations(complement.dim());
	for(int t = 0; t < iterations.dim(); t++) {
	  iterations[t] = std::make_pair(-1,-1);
	}	
	while(k != -1) {
	    //If we're through all k's, make a cone -------------------------
	    if(k >= complement.dim()) {
	      //TODO: Make a cone
	      dbgtrace << "Compatible pair: \n: Order: " << L <<"\np: " << p << endl;
	      count++;
	      k--;
	      continue;
	    }
	    //Check F_k \cap L for elements --------------------------------
	    if(iterations[k].first == -1) {
	      dbgtrace << "Iteration for " << complement[k] << " at " << iterations[k] << endl;
	      iterations[k] = std::make_pair(0,-1); 
	      int newimage = -1;
	      for(int l = 0; l < L.dim(); l++) {
		if(Fksets[k].contains(L[l])) {
		    newimage = L[l];break;
		}
	      }
	      if(newimage >= 0) {
		p[k] = newimage;
		k++;
		continue;
	      }
	    }
	    //Go through the basis elements -------------------------------
	    //Remove the last b we added
	    if(iterations[k].second != -1) {
	      L = L.slice(~scalar2set(iterations[k].second));
	    }
	    dbgtrace << "Iteration for " << complement[k] << " at " << iterations[k] << endl;
	    dbgtrace << "L is " << L << endl;
	    //Find the next valid (b,x)
	    int b = iterations[k].first;
	    int x = iterations[k].second;
	    bool found = false;
	    do {
	      x++;
	      //If we arrive after the last position, take the next b
	      if(x > L.dim()) {
		x = 0; b++;
	      }
	      //Check that b is not in L
	      bool binL = false;
	      for(int l = 0; l < L.dim(); l++) {
		if(Fklists[k][b] == L[l]) {
		  binL = true; break;
		}
	      }
	      if(binL) {
		x = 0; b++;
	      }
	      //Check that b is smaller than k (and the list size)
	      if(b >= Fklists[k].dim()) {
		break;
	      }
	      if(Fklists[k][b] >= complement[k]) {
		break;
	      }
	      //Check validity of b
	      bool isvalid = true;
	      for(int l = 0; isvalid && l < complement.dim() && complement[l] < complement[k]; l++) {
		if(Fksets[l].contains(Fklists[k][b])) {
		    //check if (p(complement[l]) > b), i.e. p(c[l]) occurs at position x or later
		    for(int y = x; y < L.dim(); y++) {
		      if(L[y] == p[l]) {
			isvalid = false; break;
		      }
		    }
		}
	      }
	      if(isvalid) {
		//If we arrive here, we found a valid position
		found = true; break;
	      }
	    } while(true);
	    //If we found a position, continue the iteration
	    if(found) {
	      p[k] = Fklists[k][b];
	      iterations[k] = std::make_pair(b,x);
	      Vector<int> newL(L.dim()+1);
		for(int y = 0; y < x; y++) newL[y] = L[y];
		newL[x] = Fklists[k][b];
		for(int y = x+1; y < newL.dim(); y++) newL[y] = L[y-1];
	      L = newL;
	      k++;
	      continue;
	    }
	    //Otherwise step down
	    else {
	      iterations[k] = std::make_pair(-1,-1);
	      k--;
	      continue;
	    }
	    
	} //End compute cones
	
      } //End for all bases
      dbgtrace << "Found " << count << " pairs " << endl;
      return perl::Object("WeightedComplex");
      
  }
  
  
  
//   void test(IncidenceMatrix<> bases, int n) {
//     //Go through all bases
//     Set<int> complete = sequence(0,n);
//     for(int i = 0; i < bases.rows(); i++) {
//       //Compute complement
//       Set<int> B = bases.row(i);
//       Set<int> C = complete - B;
//       for(Entire<Set<int> >::iterator k = entire(C); !k.at_end(); k++) {
// 	computeFk(bases, i, *k);
//       }
//     }
//   }
//   
//   void test2(IncidenceMatrix<> bases, int n, Matrix<Rational> m) {
//     //Go through all bases
//     Set<int> complete = sequence(0,n);
//     for(int i = 0; i < bases.rows(); i++) {
//       //Compute complement
//       Set<int> B = bases.row(i);
//       Set<int> C = complete - B;
//       for(Entire<Set<int> >::iterator k = entire(C); !k.at_end(); k++) {
// 	computeFkLinear(bases, i, *k,m);
//       }
//     }
//   }
  
  Function4perl(&computeMatrixBases,"computeMatrixBases(Matrix<Rational>)");
  Function4perl(&computeMatrixColoops,"computeMatrixColoops(Matrix<Rational>)");
  Function4perl(&computeFkLinear,"computeFk(IncidenceMatrix, $,$,Matrix<Rational>)");
  Function4perl(&bergman_fan,"computeBergmanFan($,IncidenceMatrix,$,Matrix<Rational>)");
//   Function4perl(&test,"test(IncidenceMatrix, $)");
//   Function4perl(&test2,"test2(IncidenceMatrix, $,Matrix<Rational>)");
  
}}

