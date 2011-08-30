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
 
 This file includes functionality for the type RationalCurve
 */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/Map.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/PowerSet.h"
#include "polymake/Array.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;

  //Documentation see perl wrapper
  int moduliDimensionFromLength(int length) {
    Rational s = sqrt(1 + 8*length);
    return (1+s) / 2;
  }
  
  //Documentation see perl wrapper
  perl::Object curveFromMetric(Vector<Rational> metric) {
    // We prepare the metric by making sure, all entries are > 0
    // and by adding Phi(sum of unit vectors) to ensure all leaves have
    // positive distance
    int n = moduliDimensionFromLength(metric.dim());
    metric += 2*ones_vector<Rational>(metric.dim());
    //For simplicity we ignore the first row and column and start counting at 1
    Matrix<Rational> d(n+1,n+1); 
    int mindex = 0;
    for(int i = 1; i < n; i++) {
      for(int j = i+1; j <= n; j++) {
	d(i,j) = metric[mindex];
	d(j,i) = metric[mindex];
	mindex++;
      }
    }
    //Now check for nonpositive entries
    for(int i = 1; i < n; i++) {
      for(int j = i+1; j <= n; j++) {
	if(d(i,j) <= 0) {
	  Vector<Rational> add = (d(i,j) +1) * ones_vector<Rational>(d.cols());
	  d.row(i) += add;
	  d.col(i) += add;
	  d(i,i) = 0;
	}
      }
    }
      
    dbgtrace << "Starting with metric matrix\n" << d << endl;
    
    
    //Result variable
    Vector<Rational> coeffs;
    Vector<Set<int> > sets;
    //Prepare vertex set, leaf map 
    dbgtrace << "Moduli dimension is " << n << endl;
    Set<int> V = sequence(1,n);
    Map<int,Set<int> > leaves;
    for(int i = 1; i <=n; i++) {
      Set<int> singleset; singleset += i;
      leaves[i] = singleset;
    }
    dbgtrace << "Starting with leaf map " << leaves << endl;
    
    //Now inductively remove pairs of vertices until only 3 are left
    while(V.size() > 3) {
      dbgtrace << "Have " << V.size() << " vertices. Reducing..." << endl;
      //Find the triple (p,q,r) that maximizes the Buneman term
      int p,q,r;
	p = q = r = 0;
      bool init = false;
      Rational max = 0;
      for(Entire<Set<int> >::iterator a = entire(V); !a.at_end(); a++) {
	for(Entire<Set<int> >::iterator b = entire(V); !b.at_end(); b++) {
	  if(*b != *a) { 
	    for(Entire<Set<int> >::iterator c = entire(V); !c.at_end(); c++) {
	      if(*c != *b && *c != *a) {
		Rational newd = d(*a,*c) + d(*b,*c) - d(*a,*b);
		dbgtrace << *a <<","<<*b<<","<<*c<<": " << newd << endl;
		if(newd > max || !init) {
		    max = newd;
		    p = *a; q = *b; r = *c;
		    init = true;
		}
	      }
	    }
	  }
	}
      } //End find maximal triple
      dbglog << "Maximal triple is (" << p << "," << q << "," << r << ")" << endl;
      
      //Compute distances to the new virtual element t
      Rational dtp = (d(p,q) + d(p,r) - d(q,r));
	dtp /= Rational(2);
	dbgtrace << "dtp: " << dtp << endl;
      Vector<Rational> dtx(d.cols()); //Again, start counting from 1
      int x = 0;
      for(Entire<Set<int> >::iterator i = entire(V); !i.at_end(); i++) {
	if(*i != p) {
	    dbgtrace << "Setting distance d(t," << *i << ")" << endl;
	    dtx[*i] = d(*i,p) - dtp;
	    if(*i != q && dtx[*i] == 0) {
	      x = *i;
	    }
	}
      }
      V = V - p; 
      V = V - q;
      dbgtrace << "Computed new distances" << endl;
      
      //Now 'add' the new vertex
      if(x> 0) {
	dbgtrace << "Attaching to vertex " << x << endl;
	leaves[x] = leaves[x] + leaves[p] + leaves[q];
	if(leaves[p].size() > 1 && leaves[p].size() < n-1) {
	  if(d(p,x) != 0) {
	    coeffs |= d(p,x);
	    sets |= leaves[p];
	  }
	}
	if(leaves[q].size() > 1 && leaves[q].size() < n-1) {
	  if(d(q,x) != 0) {
	    coeffs |= d(q,x);
	    sets |= leaves[q];
	  }
	}
      }
      else {
	dbgtrace << "Creating new vertex" << endl;
	//We update the distance matrix, since we add a new element
	d = d | zero_vector<Rational>();
	d = d / zero_vector<Rational>();
	int t = d.cols() -1;
	for(Entire<Set<int> >::iterator i = entire(V); !i.at_end(); i++) {
	    d(*i,t) = d(t,*i) = dtx[*i];
	}
	if(leaves[p].size() > 1 && leaves[p].size() < n-1) {
	  if(dtp != 0) {
	    coeffs |= dtp;
	    sets |= leaves[p];
	  }
	}
	if(leaves[q].size() > 1 && leaves[q].size() < n-1) {
	  if(dtx[q] != 0) {
	    coeffs |= dtx[q];
	    sets |= leaves[q];
	  }
	}
	//Now add the new vertex
	V += t;
	leaves[t] = leaves[p] + leaves[q];
      }
      dbgtrace << "Distance matrix\n" << d << endl;
      dbgtrace << "Leaf map\n" << leaves << endl;
    } //End while(>3)
    
    //Now treat the basic cases of size 2 and 3
    Vector<int> vAsList(V);
    if(V.size() == 3) {
      dbgtrace << "Remaining: " << vAsList << endl;
      //Solve the linear system given by the pairwise distances
      Matrix<Rational> A(3,3);
	//Create the inverse matrix of the distance relation and multiply it with the distance vectors
	A(0,0) = A(0,1) = A(1,0) = A(1,2) = A(2,1) = A(2,2) = 0.5;
	A(0,2)  = A (1,1) = A(2,0) = -0.5;
      Vector<Rational> B(3);
	B[0] = d(vAsList[0],vAsList[1]); 
	B[1] = d(vAsList[0],vAsList[2]);
	B[2] = d(vAsList[1],vAsList[2]);
      dbgtrace << "Solving " << A << "," << B << endl;
      Vector<Rational> a = A * B;
      dbgtrace << "Result: " << a << endl;
      for(int i = 0; i < 3; i++) {
	if(a[i] != 0) {
	    if(leaves[vAsList[i]].size() > 1 && leaves[vAsList[i]].size() < n-1) {
	      coeffs |= a[i];
	      sets |= leaves[vAsList[i]];
	    }
	}
      }
    }//End case size == 3
    if(V.size() == 2) {
      if(leaves[vAsList[0]].size() > 1 && leaves[vAsList[0]].size() < n-1) {
	if(d(vAsList[0],vAsList[1]) != 0) {
	  coeffs |= d(vAsList[0],vAsList[1]);
	  sets |= leaves[vAsList[0]];
	}
      }
    }
    
    //Now we're done, so we create the rational curve
    perl::Object curve("RationalCurve");
      curve.take("SETS") << sets;
      curve.take("COEFFS") << coeffs;
      curve.take("N_LEAVES") << n;
    return curve;
  }
  
  /**
    @brief Takes a linear combination of abstract n-marked curves with 1 bounded edge, described by their partitions and the corresponding edge length and computes the resulting metric
    @param IncidenceMatrix<> sets A list of partitions of {1,..,n}. May be redundant.
    @param Vector<Set<int> > coeffs A list of arbitrary rational coefficients. Superfluous coefficients are ignored, missing ones replaced by 0.
    @param int n The size of the leaf set
    @return Vector<Rational> A curve metric of length (n over 2)
  */
  Vector<Rational> metricFromCurve(IncidenceMatrix<> sets, Vector<Rational> coeffs, int n) {
    //Create distance matrix (we count from 1 for simplicity)
    Matrix<Rational> d(n+1,n+1);
    Set<int> completeSet = sequence(1,n);
    //Go through all sets
    for(int s = 0; s < sets.rows(); s++) {
      //If we have no more coefficients, stop calculating
      if(s >= coeffs.dim()) break;
      //Otherwise add the coefficients to the appropriate distances
      Rational c = coeffs[s];
      Set<int> sset = sets.row(s);
      Set<int> complement = completeSet - sset;
      for(Entire<Set<int> >::iterator selt = entire(sset); !selt.at_end(); selt++) {
	for(Entire<Set<int> >::iterator celt = entire(complement); !celt.at_end(); celt++) {
	    dbgtrace << "Adding " << c << " to distance of " << *selt << ", " << *celt << endl;
	    d(*selt,*celt) += c;
	    d(*celt, *selt) += c;
	}
      }
    }
    
    //Now convert to a vector
    Vector<Rational> result;
    for(int i = 1; i < n; i++) {
      for(int j = i+1; j <= n; j++) {
	result |= d(i,j);
      }
    }
    
    return result;
  }
  
  //Documentation see perl wrapper
  perl::Object curveFromModuli(Vector<Rational> moduliVector) {
    //Insert projection 0
    moduliVector |= 0;
    
    //Convert vector to a map
    int n = moduliDimensionFromLength(moduliVector.dim())+1;
    Matrix<Rational> d(n,n);
    int index = 0;
    for(int i = 1; i < n-1; i++) {
      for(int j = i+1; j <= n-1; j++) {
	d(i,j) = moduliVector[index];
	index++;
      }
    }
    
    //Now apply mapping
    Vector<Rational> metric;
    for(int i = 1; i < n; i++) {
      for(int j = i+1; j <= n; j++) {
	if(j == n) {
	    metric |= 0;
	}
	else {
	    metric |= (2* d(i,j));
	}
      }
    }
    
    return curveFromMetric(metric); 
  }
  
  /**
   @brief Takes a rational n-marked abstract curve and computes its representation in the matroid coordinates of the moduli space
   @param perl::Object The rational curve
   @return Vector<Rational>
   */
  Vector<Rational> moduliFromCurve(perl::Object curve) {
    //Extract values
    IncidenceMatrix<> sets = curve.give("SETS");
    Vector<Rational> coeffs = curve.give("COEFFS");
    int n = curve.give("N_LEAVES");
    
    //Create edge index map (i,j) -> vector index
    Matrix<Rational> E(n,n); int index = 0;
    for(int i = 1; i < n-1; i++) {
      for(int j = i+1; j <= n-1; j++) {
	E(i,j) = E(j,i) = index;
	index++;
      }
    }
    
    //Compute ambient dimension of moduli space
    int raydim = (n*(n-1))/2 - n;
    Set<int> completeSet = sequence(1,n);
    
    Vector<Rational> result(raydim);
    Vector<Rational> onlyones(ones_vector<Rational>(raydim));
    
    //Map each set to matroid coordinates with appropriate coefficient
    for(int s = 0; s < sets.rows(); s++) {
      Set<int> sset = sets.row(s);
      //Make sure the set does not contain n
      if(sset.contains(n)) sset = completeSet - sset;
      //Now create the flat vector for the complete graph on vertices in sset
      Vector<Rational> slist(sset);
      for(int i = 0; i < slist.dim(); i++) {
	for(int j = i+1; j < slist.dim(); j++) {
	  int edgeindex = E(slist[i],slist[j]);
	  if(edgeindex < raydim) {
	    result[E(slist[i],slist[j])] -= coeffs[s];
	  }
	  else { //If edgeindex = raydim, this is the coordinate that we modded out
	    result += coeffs[s] * onlyones;
	  }
	}
      }
    }
    
    return result;
  }
  
  //Documentation see perl wrapper
  perl::ListReturn curveFromModuliMatrix(Matrix<Rational> m) {
    perl::ListReturn result;
    
    for(int i = 0; i < m.rows(); i++) {
      result << curveFromModuli(m.row(i));
    }
    
    return result;
  }
  
  /**
   @brief Takes three integer values and checks whether two of them are equal and >= than the third
   */
  inline bool fpcCheck(int a, int b, int c) {
    if(a == b && a >= c) return true;
    if(a == c && a >= b) return true;
    if(b == c && b >= a) return true;
    return false;
  }
  
  //Documentation see perl wrapper
  perl::ListReturn curveFromMetricMatrix(Matrix<Rational> m) {
    perl::ListReturn result;
    
    for(int i = 0; i < m.rows(); i++) {
      result << curveFromMetric(m.row(i));
    }
    
    return result;
  }
  
  //Documentation see perl wrapper
  perl::ListReturn testFourPointCondition(Vector<Rational> v) {
    //Convert metric into map
    int n = moduliDimensionFromLength(v.dim());
    Matrix<Rational> d(n+1,n+1);
    
    int mindex = 0;
    for(int i = 1; i < n; i++) {
      for(int j = i+1; j <= n; j++) {
	d(i,j) = d(j,i) = v[mindex];
	mindex++;
      }
    }
    
    //First we test all 4-element subsets
    Set<int> complete = sequence(1,n);
    Array<Set<int> > fours =  pm::Subsets_of_k<Set<int> > ( complete,4 );
    for(int f = 0; f < fours.size(); f++) {
      Vector<int> l(fours[f]);
      int a = d(l[0],l[1]) + d(l[2],l[3]);
      int b = d(l[0],l[2]) + d(l[1],l[3]);
      int c = d(l[0],l[3]) + d(l[1],l[2]);
      //Check that two of a,b,c are equal and not less than the third
      if(!fpcCheck(a,b,c)) {
	perl::ListReturn fault;
	  fault << l[0] << l[1] << l[2] << l[3];
	return fault;
      }
    }
    //Now we check all 3-element subsets
    Array<Set<int> > threes = pm::Subsets_of_k<Set<int> >(complete, 3);
    for(int f = 0; f < threes.size(); f++) {
      Vector<int> l(threes[f]);
      //Now check the three possibilities, where the fourth element is equal to any of the three
      for(int t = 0; t < l.size(); t++) {
	int a = d(l[0],l[1]) + d(l[2],l[t]);
	int b = d(l[0],l[2]) + d(l[1],l[t]);
	int c = d(l[0],l[t]) + d(l[1],l[2]);
	//Check that two of a,b,c are equal and not less than the third
	if(!fpcCheck(a,b,c)) {
	  perl::ListReturn fault;
	    fault << l[0] << l[1] << l[2] << l[t];
	  return fault;
	}
      }
    }
    //Now we check all 2-element subsets
    Array<Set<int> > twos = pm::Subsets_of_k<Set<int> >(complete, 2);
    for(int f = 0; f < twos.size(); f++) {
      Vector<int> l(twos[f]);
      //We have three possibilites for the other two z,t: t=x,z=y or t=z=x or t=z=y
      for(int p = 1; p <= 3; p++) {
	int t = p < 3? l[0] : l[1];
	int z = p != 2? l[1] : l[0];
	int a = d(l[0],l[1]) + d(z,t);
	int b = d(l[0],z) + d(l[1],t);
	int c = d(l[0],t) + d(l[1],z);
	//Check that two of a,b,c are equal and not less than the third
	if(!fpcCheck(a,b,c)) {
	  perl::ListReturn fault;
	    fault << l[0] << l[1] << z << t;
	  return fault;
	}
      }
    }
    perl::ListReturn result;
    return result;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------

  UserFunction4perl("# @category Tropical geometry"
		    "# Takes a positive length (of a vector) and assumes it is of the form (n over 2)"
		    "# It then computes n"
		    "# @param Int k The length = (n over 2)"
		    "# @return Int n. If k is not of the form (n over 2), the result is arbitrary",
		    &moduliDimensionFromLength, "moduliDimensionFromLength($)");
		     
  UserFunction4perl("# @category Tropial geometry"
		    "# Takes a vector from Q^(n over 2) that describes an n-marked rational abstract"
		    "# curve as a distance vector between its leaves. It then computes the "
		    "# curve corresponding to this vector."
		    "# @param Vector<Rational> v A vector of length (n over 2). Its entries are "
		    "# interpreted as the distances d(i,j) ordered lexicographically according to i,j. However, they need not be positive, as long as v is equivalent to a proper "
		    "# metric modulo leaf lengths."
		    "# @return RationalCurve",
		    &curveFromMetric,"rational_curve_from_metric(Vector<Rational>)");
  
  UserFunction4perl("# @category Tropical geometry"
		    "# Takes a vector from Q^((n over 2) - n) that lies in M_0,n (in its matroid coordinates "
		    "# and computes the corresponding rational curve."
		    "# @param Vector<Rational> v A vector in the moduli space"
		    "# @return RationalCurve",
		    &curveFromModuli,"rational_curve_from_moduli(Vector<Rational>)");
  
  UserFunction4perl("# @category Tropical geometry"
		    "# Takes a matrix whose rows are elements in the moduli space M_0,n in matroid "
		    "# coordinates. Returns a list, where the i-th element is the curve corr. to "
		    "# the i-th row in the matrix"
		    "# @param Matrix<Rational> m"
		    "# @return RationalCurve : An array of RationalCurves",
		    &curveFromModuliMatrix, "rational_curve_list_from_moduli(Matrix<Rational>)");
  
  UserFunction4perl("# @category Tropical geometry"
		    "# Takes a matrix whose rows are metrics of rational n-marked curves."
		    "# Returns a list, where the i-th element is the curve corr. to "
		    "# the i-th row in the matrix"
		    "# @param Matrix<Rational> m"
		    "# @return RationalCurve : An array of RationalCurves",
		    &curveFromMetricMatrix, "rational_curve_list_from_metric(Matrix<Rational>)");
  
  UserFunction4perl("# @category Tropical geometry" 
		    "# Takes a metric vector in Q^{(n over 2)} and checks whether it fulfills "
		    "# the four-point condition, i.e. whether it lies in M_0,n. More precisely "
		    "# it only needs to be equivalent to such a vector"
		    "# @param Vector<Rational> v The vector to be checked"
		    "# @return Int A quadruple (array) of indices, where the four-point condition "
		    "# is violated or an empty list, if the vector is indeed in M_0,n",
		    &testFourPointCondition, "testFourPointCondition(Vector<Rational>)");
		    
  
  Function4perl(&metricFromCurve, "metric_from_curve(IncidenceMatrix, Vector<Rational>, $)");
  Function4perl(&moduliFromCurve, "moduli_from_curve(RationalCurve)");
  
}}