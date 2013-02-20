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
#include "polymake/Map.h"
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
    
    //dbgtrace << "Delta: " << delta << endl;
    
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
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object forgetful_map(int n, Set<int> leaves_to_forget) {
    
    //First check, that the leaves are in (1,..,n)
    if((leaves_to_forget * sequence(1,n)).size() < leaves_to_forget.size()) {
      throw std::runtime_error("Cannot compute forgetful map: The forgotten leaves should be in {1,..,n}");
    }
    
    //Compute domain and image dimension
    int domain_dim = (n*n - 3*n)/2;
    int small_n = n - leaves_to_forget.size();
    int image_dim = (small_n * small_n - 3*small_n)/2;
    
    
    
    //Check if we forget so many leaves that we get the zero map
    if(small_n <= 3) {
      perl::Object result("Morphism");
	result.take("MATRIX") << Matrix<Rational>(0,domain_dim);
      return result;
    }
    //Check if we don't forget anything at all
    if(leaves_to_forget.size() == 0) {
      perl::Object result("Morphism");
	Matrix<Rational> um = unit_matrix<Rational>(domain_dim);
	result.take("MATRIX") << um;
      return result;
    }
    
    //Prepare map mapping remaining leaves to {1,..,n-|leaves_to_forget|}
    Map<int,int> remaining_map;
    int next_index = 1;
    for(int i = 1; i <= n; i++) {
      if(!leaves_to_forget.contains(i)) {
	remaining_map[i] = next_index;
	next_index++;
      }
    }
    
    //Compute the unit vectors corresponding to the 1-edge-flats in the target M_0,small_n
    Map<int, Map<int,int> > unit_edges;
    int unit_index = 0;
    for(int i = 1; i < small_n-2; i++) {
      unit_edges[i] = Map<int,int>();
      for(int j = i+1; j < small_n; j++) {
	(unit_edges[i])[j] = unit_index;
	unit_index++;
      }
    }
    
    //Prepare matrix representing forgetful map
    Matrix<Rational> ffm(image_dim,0);
    
    //Compute image of each neagtive unit vector -e_i, correpsonding to all v_{k,l}
    //with 1 <= k < l < n and {k,l} != {n-2,n-1}
    for(int k = 1; k < n-2; k++) {
      for(int l = k+1; l < n; l++) {
	Set<int> klset; klset += k; klset +=l;
	//If klset contains forgotten leaves, the ray is mapped to zero
	if( (leaves_to_forget * klset).size() > 0) {
	  ffm |= zero_vector<Rational>(image_dim);
	}
	else {
	  //First we compute the new number of the leaves in the image
	  int newk = remaining_map[k];
	  int newl = remaining_map[l];
	  
	  //The easy case is that newl is not the new maximal leaf
	  if(newl < small_n) {
	    //If (newk,newl) != (small_n - 2, small_n-1), the image is just the corr. unit vector,
	    //otherwise its minus the sum of all unit_vectors
	    if(newk < small_n -2) {
	      ffm |= unit_vector<Rational>(image_dim, (unit_edges[newk])[newl]);
	    }
	    else {
	      ffm |= (- ones_vector<Rational>(image_dim));
	    }
	  }//END if newl < small_n
	  else {
	    //If the new set contains the maximal leaf, we have to take the complement to 
	    //compute the moduli coordinates
	    Set<int> complement = (sequence(1,small_n) - newk) - newl;
	    //Compute moduli coordinates as usual: Sum up over all flat vectors of the 
	    //1-edge flats, i.e. all pairs contained in complement
	    Vector<Rational> ray_image(image_dim);
	    Vector<int> slist(complement);
	    for(int i = 0; i < slist.dim(); i++) {
	      for(int j = i+1; j < slist.dim(); j++) {
		//If its the last edge (small_n-2, small_n-1), we add the ones_vector, otherwise, 
		//we  substract the corr. unit_vector
		if(slist[i] == small_n-2) {
		  ray_image += ones_vector<Rational>(image_dim);
		}
		else {
		  ray_image[(unit_edges[slist[i]])[slist[j]]] -= 1;
		}
	      }
	    }//END compute ray image
	    
	    //Recall that we computed the image of -e_i, so we have - ray_image as image
	    ffm |= (- ray_image);
	    
	  }//END if newl == small_n
	  
	}//END compute image of ray
      }//END iterate l
    }//END iterate k
    
    
    //Return result
    perl::Object result("Morphism");
      result.take("MATRIX") << ffm;
      
    return result;
    
  }//END function forgetful_map
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category The moduli space M_0,n"
		    "# This creates the i-th evaluation function on M_0,n^(lab)(R^r,Delta)"
		    "# (which is actually realized as M_0,(n+|Delta|) x R^r)"
		    "# @param Int n The number of marked points"
		    "# @param Int r The dimension of the embedding space"
		    "# @param Matrix<Rational> Delta The directions of the unbounded edges (given as row vectors)"
		    "# @param Int i The index of the marked point that should be evaluated. i "
		    "# should lie in between 1 and n"
		    "# Note that the i-th marked point is realized as the |Delta|+i-th leaf in M_0,(n+|Delta|)"
		    "# and that the R^r - coordinate is interpreted as the position of the n-th leaf. "
		    "# In particular, ev_n is just the projection to the R^r-coordinates"
		    "# @return Morphism ev_i. Its domain is the ambient space of the moduli space "
		    "# in matroid coordinates cross R^r",
		    &evaluation_map,"evaluation_map($,$,Matrix<Rational>,$)");
  
  UserFunction4perl("# @category The moduli space M_0,n"
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
  
  UserFunction4perl("# @category Morphisms and functions"
		    "# This creates a linear projection from R^n to a given set of coordinates"
		    "# @param Int n The dimension of the domain"
		    "# @param Set<Int> s The set of coordinates to which this map should project (starting "
		    "# the count at 0)"
		    "# @return Morphism The corresponding projection as a global linear map",
		    &projection_map, "projection_map($, Set<Int>)");
  
  UserFunction4perl("# @category Morphisms and functions"
		    "# This computes the projection from R^n to R^m (for m < n) onto the first m coordinates"
		    "# @param Int n Dimension of domain"
		    "# @param Int m Dimension of image"
		    "# @return Morphism The corresponding projection",
		    &projection_map_default, "projection_map($,$)");
  
  UserFunction4perl("# @category The moduli space M_0,n"
		    "# This computes the forgetful map from the moduli space M_0,n to M_0,(n-|S|)"
		    "# @param Int n The number of leaves in the moduli space M_0,n"
		    "# @param Set<Int> S The set of leaves to be forgotten. Should be a subset of (1,..,n)"
		    "# @return Morphism The forgetful map. It will identify the remaining leaves "
		    "# i_1,..,i_(n-|S|) with the leaves of M_0,(n-|S|) in order of their size",
		    &forgetful_map,"forgetful_map($,Set<Int>)");
  
  
}}