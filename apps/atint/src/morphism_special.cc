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
 Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 
 This file contains creation functions for special types of morphisms
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/moduli.h"

namespace polymake { namespace atint { 
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object evaluation_map(int n, int r, Matrix<Rational> delta, int i) {
    if(n <= 0 || r <= 0 || delta.rows() <= 0 || i <= 0 || i > n) {
      throw std::runtime_error("Cannot create evaluation map: Invalid parameters");
    }
    
    //The matrix of ev_i is pr(R^r) + G*pr(M_0,N), where
    //G is the matrix evaluating the differences between the marked points
    //Since we take ev_n as the special evaluation, the map reads under moduli coordinates:
    // R^((n-1) over 2) -> R^r, x |-> - sum_(k=1)^|delta| x_{ik} v_k (where v_k is the k-th
    // row of delta). 
    //Modding out the lineality space is realized by simply forgetting the last column of the
    //corresponding matrix
    
    //Projection matrices
    int N = n + delta.rows();
    int modulidim = (N*(N-1))/2 - N;
    
    Matrix<Rational> projR = Matrix<Rational>(r,modulidim) | unit_matrix<Rational>(r);
    
    Matrix<Rational> projM = unit_matrix<Rational>(modulidim) | Matrix<Rational>(modulidim,r);
    
    Matrix<Rational> G(r,modulidim);
    if(i < n) {    //If i is the special leaf, G is the zero map
      //Create index map (i,j) -> index(i,j)
      Matrix<int> E = pair_index_map(N-1);
      int evalIndex = delta.rows() + i-1;
      //Set all rows corr. to coordinate (n+i,k) to -  v_k
      for(int k = 0; k < delta.rows(); k++) {
	if(E(evalIndex,k) < G.cols()) {
	  G.col(E(evalIndex,k)) = - delta.row(k);
	}
      }
    }
    
    Matrix<Rational> map_matrix = projR + G*projM;
    
    perl::Object morphism("Morphism");
      morphism.take("MATRIX") << map_matrix;
      
    return morphism;
    
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object evaluation_map_d(int n, int r, int d, int i) {
    if(n <= 0 || r <= 0 || d <= 0 || i <= 0 || i > n) {
      throw std::runtime_error("Cannot create evaluation map: Invalid parameters");
    }
    
    //Create standard d-fold direction matrix
    Matrix<Rational> delta(0,r);
    for(int x = 0; x <= r; x++) {
      Vector<Rational> append;
      if(x == 0) append = ones_vector<Rational>(r);
      else append = -unit_vector<Rational>(r,x-1);
      for(int j = 1; j <= d; j++) {
	delta /= append;
      }
    }
    
    dbgtrace << "Delta: " << delta << endl;
    
    return evaluation_map(n,r,delta,i);
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object projection_map(int n, Set<int> coords) {
    
    //Create matrix
    Matrix<Rational> proj_matrix(coords.size(), n);
    int image_index = 0;
    for(Entire<Set<int> >::iterator c = entire(coords); !c.at_end(); c++) {
      if(*c >= n) {
	throw std::runtime_error("Cannot create projection: Image dimension larger than domain dimension");
      }
      proj_matrix.col(*c) = unit_vector<Rational>(n,image_index);
      image_index++;
    }
    
    perl::Object result("Morphism");
      result.take("MATRIX") << proj_matrix;
      
    return result;
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object projection_map_default(int n, int m) {
    if(m > n) {
      throw std::runtime_error("Cannot create projection: Image dimension larger than domain dimension");
    }
    return projection_map(n, sequence(0,m));
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry/Morphisms"
		    "# This creates the i-th evaluation function on M_0,n^(lab)(R^r,Delta)"
		    "# (which is actually realized as M_0,(n+|Delta|) x R^r)"
		    "# @param Int n The number of marked points"
		    "# @param Int r The dimension of the embedding space"
		    "# @param Matrix<Rational> Delta The directions of the unbounded edges"
		    "# @param Int i The index of the marked point that should be evaluated. i "
		    "# should lie in between 1 and n"
		    "# @return Morphism ev_i. Its domain is the ambient space of the moduli space "
		    "# in matroid coordinates cross R^r",
		    &evaluation_map,"evaluation_map($,$,Matrix<Rational>,$)");
  
  UserFunction4perl("# @category Tropical geometry/Morphisms"
		    "# This creates the i-th evaluation function on M_0,n^(lab)(R^r,d)"
		    "# (which is actually realized as M_0,(n+d(r+1)) x R^r)"
		    "# This is the same as calling the function"
		    "# evaluation_map(Int,Int,Matrix<Rational>,Int) with the standard d-fold"
		    "# degree as matrix"
		    "# @param Int n The number of marked points"
		    "# @param Int r The dimension of the embedding space"
		    "# @param Int d The degree of the embedding. The direction matrix will be"
		    "# the standard d-fold directions (first (1,...1), then the inverted unit vectors)"
		    "# @param Int i The index of the marked point that should be evaluated. i "
		    "# should lie in between 1 and n"
		    "# @return Morphism ev_i. Its domain is the ambient space of the moduli space "
		    "# in matroid coordinates cross R^r",
		    &evaluation_map_d,"evaluation_map($,$,$,$)"); 
  
  UserFunction4perl("# @category Tropical geometry/Morphisms"
		    "# This creates a linear projection from R^n to a given set of coordinates"
		    "# @param Int n The dimension of the domain"
		    "# @param Set<Int> The set of coordinates to which this map should project (starting "
		    "# the count at 0)"
		    "# @return Morphism The corresponding projection as a global linear map",
		    &projection_map, "projection_map($, Set<Int>)");
  
  UserFunction4perl("# @category Tropical geometry/Morphisms"
		    "# This computes the projection from R^n to R^m (for m < n) onto the first m coordinates"
		    "# @param Int n Dimension of domain"
		    "# @param Int m Dimension of image"
		    "# @return Morphism The corresponding projection",
		    &projection_map_default, "projection_map($,$)");
  
  
}}