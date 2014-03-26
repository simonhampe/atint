/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 * 
 * ---
 * Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 * 
 * 
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Array.h"
#include "polymake/Set.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/fan_piecewise_divisor.h"
#include "polymake/atint/normalvector.h"

namespace polymake { namespace atint { 
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  void cln(perl::Object fan) {
    
    
    
    //dbgtrace << "Extracting properties" << endl;
    
    //Extract basic properties of fan
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    IncidenceMatrix<> codimInc = fan.give("CODIM_1_IN_MAXIMAL_CONES");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    
    //This will contain the result
    Map<int, Map<int, Vector<Integer> > > latticeNormals;
    if(codimOneCones.rows() == 0) {
//       fan.take("LATTICE_NORMALS") << latticeNormals;
      return;
    }
       
    //This equation is added to all cone linear span matrices (and stands for intersecting
    // with (x0 = 1), if we use homogeneous coordinates, to ensure that the lattice normal
    // is of the form (0,...)
    Vector<Rational> intereq(unit_vector<Rational>(ambient_dim,0));
    
    //dbgtrace << "Computing linear spans" << endl;
    
    //Compute all the linear spans of the cones before, so we don't do it several times
    Vector<Matrix<Rational> > codimone;
    Vector<Matrix<Rational> > maximal;
    for(int facet = 0; facet < codimOneCones.rows(); facet++) {
      codimone |= null_space(rays.minor(codimOneCones.row(facet),All) / linspace);
    }
    for(int mcone = 0; mcone < maximalCones.rows(); mcone++) {
      maximal |=  null_space(rays.minor(maximalCones.row(mcone),All) / linspace);
    }
    
    //dbgtrace << "Computing lattice normals" << endl;
    
    //Go through all codim one faces
    for(int facet = 0; facet < codimone.dim(); facet++) {
      //Create map for this facet
      latticeNormals[facet] = Map<int, Vector<Integer> >();
      
      //Construct the dual of the linear span of the cone:
      //= kernel of generators
      Matrix<Rational> facetmatrix = codimone[facet];
      if(uses_homog) {
	facetmatrix = facetmatrix / intereq;
      }
      Set<int> adjacentCones = codimInc.row(facet);
      //Go through all adjacent cones
      for(Entire<Set<int> >::iterator e=entire(adjacentCones); !e.at_end(); ++e) {
	Matrix<Rational> maxmatrix = maximal[*e];
	if(uses_homog) {
	  maxmatrix = maxmatrix / intereq;
	}
	//Extract the additional ray by taking any index in the maximal cone
	// and not in the codim 1 cone
	int additionalRayIndex = *(maximalCones.row(*e) - codimOneCones.row(facet)).begin();
	Vector<Rational> additionalRay(rays.row(additionalRayIndex));
	
	//Compute normalvector
	(latticeNormals[facet])[*e] = latticeNormal(facetmatrix,maxmatrix,additionalRay);
      }
    }
    
//     fan.take("LATTICE_NORMALS") << latticeNormals;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  void cpl(perl::Object fan) {
    
    //dbgtrace << "Extracting properties" << endl;
    
    //Extract basic properties of fan
    IncidenceMatrix<> codimInc = fan.give("CODIM_1_IN_MAXIMAL_CONES");
      IncidenceMatrix<> codimInMaximal = T(codimInc);
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    IncidenceMatrix<> codimOneCones = fan.give("CODIM_1_FACES");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linspace = fan.give("LINEALITY_SPACE");
    Matrix<Integer> lattice_generators = fan.give("LATTICE_GENERATORS");;
    IncidenceMatrix<> lattice_bases = fan.give("LATTICE_BASES");
    bool uses_homog = fan.give("USES_HOMOGENEOUS_C");
    
    //This will contain the result
    Map<int, Map<int, Vector<Integer> > > latticeNormals;
    if(codimOneCones.rows() == 0) {
//       fan.take("LATTICE_NORMALS") << latticeNormals;
      return;
    }
    
    //Create maps for facets
    for(int co = 0; co < codimOneCones.rows(); co++) {
      latticeNormals[co] = Map<int,Vector<Integer> >();
    }
    
    //Compute normal vectors
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      //Compute the projection
      //dbgtrace << "Computing projection matrix for cone " << mc << endl;
      Matrix<Integer> lb = lattice_generators.minor(lattice_bases.row(mc),All) ;
      //This matrix is used to compute the preimage of the projection normal
      Matrix<Integer> inverse = T(lb);
	if(uses_homog) inverse = zero_vector<Integer>() | inverse;
      Set<int> I = basis_cols(lb);
      Matrix<Rational> S = inv(T(lb.minor(All,I)));
      Matrix<Rational> P(S.rows(),lattice_generators.cols());
      P.minor(All,I) = S;
      
      //dbgtrace << "Projection matrix reads: " << P << endl;
      //dbgtrace << "lb: " << lb << endl;
      
      //Project maximal cone rays and lineality and compute span
      Matrix<Rational> mc_rays(rays.minor(maximalCones.row(mc),All));
	//Identify vertices
	Set<int> vertices;
	if(uses_homog) {
	  for(int r = 0; r < mc_rays.rows(); r++) {
	    if(mc_rays(r,0) == 1) vertices += r;
	  }
	}
	//Project rays
	mc_rays = mc_rays * T(P);
	//Add homog. "1" where necessary
	if(uses_homog) {
	  mc_rays = zero_vector<Rational>() | mc_rays;
	  mc_rays.col(0).slice(vertices) = ones_vector<Rational>(vertices.size());
	}
	
      Matrix<Rational> proj_lineality = linspace * T(P);
      if(uses_homog && proj_lineality.rows() > 0) {
	proj_lineality = zero_vector<Rational>() | proj_lineality;
      }
      
      Matrix<Rational> mc_span = null_space(mc_rays / proj_lineality);
      
      if(uses_homog) {
	//Intersect with x0 = 0 in homog. coordinates
	mc_span /= unit_vector<Rational>(lb.rows()+1,0);
      }      
      
      //Now iterate the facets
      Set<int> facets = codimInMaximal.row(mc);
      for(Entire<Set<int> >::iterator co = entire(facets); !co.at_end(); co++) {
	//dbgtrace << "Computing for facet " << *co << endl;
	//Project facet
	Matrix<Rational> co_rays = rays.minor(codimOneCones.row(*co),All);
	  //Identify vertices
	  Set<int> co_vertices;
	  if(uses_homog) {
	    for(int r = 0; r < co_rays.rows(); r++) {
	      if(co_rays(r,0) == 1) co_vertices += r;
	    }
	  }
	  //Project rays
	  co_rays = co_rays * T(P);
	  //Add homog. "1" where necessary
	  if(uses_homog) {
	    co_rays = zero_vector<Rational>() | co_rays;
	    co_rays.col(0).slice(co_vertices) = ones_vector<Rational>(co_vertices.size());
	  }
	  
	Matrix<Rational> co_span = null_space(co_rays / proj_lineality);
	if(uses_homog) {
	  //Intersect with x0 = 0 in homog. coordinates
	  co_span /= unit_vector<Rational>(co_span.cols(),0);
	}
	
	//Extract the additional ray by taking any index in the maximal cone
	// and not in the codim 1 cone
	int additionalRayIndex = *(maximalCones.row(mc) - codimOneCones.row(*co)).begin();
	Vector<Rational> additionalRay(rays.row(additionalRayIndex) * T(P));
	if(uses_homog) {
	  Rational add = rays(additionalRayIndex,0) == 1? 1 : 0;
	  additionalRay = add | additionalRay;
	}
	
	//dbgtrace << "Mc: " << mc_rays << endl;
	//dbgtrace << "Mc-span" << mc_span << endl;
	//dbgtrace << "Co: " << co_rays << endl;
	//dbgtrace << "co-span: " << co_span << endl;
	//dbgtrace << "ar: " << additionalRay << endl;
	
	//Compute facet normal in projection
	Vector<Integer> normal = latticeNormal(co_span, mc_span,additionalRay);
	
	//dbgtrace << "ln: " << normal << endl;
	
	//Compute preimage
	(latticeNormals[*co])[mc] = inverse * normal;
	//dbgtrace << "result: " << (latticeNormals[*co])[mc] << endl;
	
      }//END iterate facets
    }//END iterate maximal cones
    
//     fan.take("LATTICE_NORMALS") << latticeNormals; 
    
    /*int dim = sigma.give("CONE_DIM");
    Set<int> I = basis_cols(latticeBasis);
    Matrix<Rational> S = inv(T(latticeBasis.minor(All,I)));
    Matrix<Rational> P(dim,latticeBasis.cols());
    P.minor(All,I) = S;
    
    pm::cout << "Trafo " << P << endl;
    
    Matrix<Rational> taurays = tau.give("RAYS");
      taurays = T(P * T(taurays));
    Matrix<Rational> sigmarays = sigma.give("RAYS");
      sigmarays = T(P * T(sigmarays));
      
    pm::cout << "Tau: " << taurays << endl;
    pm::cout << "sigma: " << sigmarays << endl;
      
    perl::Object trtau("polytope::Cone");
      trtau.take("RAYS") << taurays;
    perl::Object trsigma("polytope::Cone");
      trsigma.take("RAYS") << sigmarays;
      
    Vector<Integer> ln = latticeNormalByCone(trtau,trsigma);
    pm::cout << "ln: " << ln << endl;
    
    return T(latticeBasis)*ln;*/
    
    
  }
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&cln,"cln(WeightedComplex)");
  Function4perl(&cpl,"cpl(WeightedComplex)");
  
//   using namespace atintlog::donotlog;
//   //using namespace atintlog::dolog;
// //   using namespace atintlog::dotrace;
//   
//   //Returns all possible vectors of length "length", filled with elements in "content"
//   Matrix<int> fill_value_vector(int length, Vector<int> content) {
//     Matrix<int> result(0,length);
//     if(content.dim() == 1) {
//       result /= (content[0] * ones_vector<int>(length));
//       return result;
//     }
//     //Go through all subsets
//     Array<Set<int> > subsets = all_subsets(sequence(0,length));
//     for(int s = 0; s < subsets.size(); s++) {
//       //We want each entry at least once, so we leave enough space
// //       if(length - subsets[s].size() >= content.dim() -1 && subsets[s].size() > 0) {
// 	//Fill subset s with first entry
// 	Vector<int> v(length);
// 	v.slice(subsets[s]) = content[0] * ones_vector<int>(subsets[s].size());
// 	//Continue recursively:
// 
// 	Matrix<int> recursive = 
// 	  fill_value_vector(length - subsets[s].size(), content.slice(~scalar2set(0)));
// 	for(int r = 0; r < recursive.rows(); r++) {
// 	  v.slice(~subsets[s]) = recursive.row(r);
// 	  result /= v;
// 	}
//       
// //       }
//     }
//     return result;
//   }
// 
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   Matrix<Rational> eq_matrix(perl::Object div) {
//     
//     Set<int> first_half = sequence(0,3);
//     Set<int> second_half = sequence(3,3);
//     //Extract values
//       Matrix<Rational> rays = div.give("RAYS");
//       IncidenceMatrix<> codim = div.give("CODIM_1_FACES");
//       IncidenceMatrix<> codimInMax = div.give("CODIM_1_IN_MAXIMAL_CONES");
//       Map<int, Map<int, Vector<Rational> > > lnFunctionVector = div.give("LATTICE_NORMAL_FCT_VECTOR");
//       Matrix<Rational> lsumFunctionVector = div.give("LATTICE_NORMAL_SUM_FCT_VECTOR");
//       Vector<Integer> weights = div.give("TROPICAL_WEIGHTS");
//       
//       //Find the four diagonal rays (if they exist)
//       Set<int> diag_rays;
//       for(int r = 0; r < rays.rows(); r++) {
// 	if(rays.row(r).slice(first_half) == rays.row(r).slice(second_half)) diag_rays += r;
//       }
//       
//       //Find the six cones that consist of diagonal rays
//       Set<int> diag_cones;
//       for(int c = 0; c < codim.rows(); c++) {
// 	if((codim.row(c) * diag_rays).size() == 2) diag_cones += c;
//       }
//       
// //       dbgtrace << "Diag rays " << diag_rays << ", " << "Diag cones " << diag_cones << endl;
//       
//       //Now check if we can cut out the diagonal with the appropriate linear equation system
//       Matrix<Rational> eq(0,12);
//       for(int co = 0; co < codim.rows(); co++) {
// 	Vector<Rational> v = lsumFunctionVector.row(co);
// 	Set<int> adjacent = codimInMax.row(co);
// 	for(Entire<Set<int> >::iterator mc = entire(adjacent); !mc.at_end(); mc++) {
// 	    v += weights[*mc] * (lnFunctionVector[co])[*mc];
// 	}
// 	eq /= v;
//       }
// //       dbgtrace << "Matrix " << eq << endl;
//       Vector<Rational> desired = zero_vector<Rational>(eq.rows());
// 	desired.slice(diag_cones) = ones_vector<Rational>(6);
// 	
// 	return (eq | desired);
//   }
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   Matrix<Rational> test_diagonal_combinations(perl::Object D) {
//     
// //     Vector<Rational> first_index;
// //     Vector<Rational> second_index;
// //     Vector<Rational> third_index;
// //     for(int i = 0; i < 4; i++) {
// //       for(int j = i+1; j < 4; j++) {
// // 	for(int k = j+1; k < 4; k++) {
// // 	  first_index |= i;
// // 	  second_index |= j;
// // 	  third_index |= k;
// // 	}
// //       }
// //     }
//     
//     int translate = 4;
//     
//     pm::cout << "Creating possibilities" << endl;
//     Matrix<int> possibilities = fill_value_vector(12, Vector<int>(sequence(0,3)));
//     
//     
//     Matrix<Rational> result(0,12);
//     for(int p = 0; p < possibilities.rows(); p++) {
//       pm::cout << "Testing " << p << " of " << possibilities.rows() << ", " << result.rows() << " found " << endl;
//       
//       int i = 0;
//       Vector<Set<int> > piecewise;
//       Vector<Integer> signs;
//       for(int a = 0; a < 4; a++) {
// 	for(int b = a+1; b < 4; b++) {
// 	   for(int c = b+1; c < 4; c++) {
// 	    Set<int> cone;
// 	    cone += (a + translate * possibilities(p,i));
// 	    cone += (b + translate * possibilities(p,i+1));
// 	    cone += (c + translate * possibilities(p,i+2));
// 	    piecewise |= cone;
// 	    signs |= (possibilities(p,i) == 2? -1 : 1) * (possibilities(p,i+1) == 2? -1 : 1) * (possibilities(p,i+2) == 2? -1 : 1);
// 	    i+=3;
// 	   }
// 	}
//       }
//       
//       perl::Object div  = piecewise_divisor(D, piecewise, signs);
//       
//       //Extract values
//       Vector<Integer> weights = div.give("TROPICAL_WEIGHTS");
//       Matrix<Rational> rays = div.give("RAYS");
//       if(weights == ones_vector<Integer>(4)) {
// 	if(rays.minor(All,sequence(0,3)) == rays.minor(All,sequence(3,3))) {
// 	    result /= possibilities.row(p);
// 	    pm::cout << "WORKS" << endl;
// 	}
//       }
//       
//       
//       
//       
//     }//END go through possibilities
//     
//     
//     return result;
//   }
//   
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   Matrix<Rational> test_sign_combinations(perl::Object d, IncidenceMatrix<> cones) {
//     
//     Vector<int> signs; signs |= 1; signs |= -1;
//     Matrix<int> possibilities = fill_value_vector(10,signs);
//     
//     Matrix<Rational> result(0,10);
//     
//     for(int s = 0; s < possibilities.rows(); s++) {
//       pm::cout << "Testing " << s << " of " << possibilities.rows() << endl;
//       
//       perl::Object div  = piecewise_divisor(d, cones, Vector<Integer>(possibilities.row(s)));
//       
//       //Extract values
//       Vector<Integer> weights = div.give("TROPICAL_WEIGHTS");
//       Matrix<Rational> rays = div.give("RAYS");
//       if(weights == ones_vector<Integer>(6)) {
// 	if(rays.minor(All,sequence(0,3)) == rays.minor(All,sequence(3,3))) {
// 	    result /= possibilities.row(s);
// 	    pm::cout << "WORKS" << endl;
// 	}
//       }
//       
//     }
//     
//     return result;
//     
//   }

  
//   
  // ------------------------- PERL WRAPPERS ---------------------------------------------------

 
//   Function4perl(&eq_matrix, "l32em(WeightedComplex)");
//   Function4perl(&test_diagonal_combinations, "l32td(WeightedComplex)");
//   Function4perl(&test_sign_combinations, "l43ts(WeightedComplex,IncidenceMatrix)");
  
}}