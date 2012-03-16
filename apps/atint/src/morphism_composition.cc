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
 * Contains a function to compute the composition of two morphisms (or actually: of any two 
 * piecewise affine linear functions)
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/linalg.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/morphism_pullback.h"

namespace polymake { namespace atint { 
    
  using polymake::polytope::cdd_interface::solver;
  
//   using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
    /**
   @brief Computes the composition g(f) of two morphisms f and g (as in f:X->Y, g:Y->Z). Actually, f and g can also be piecewise linear maps on their domains, the method will work equally well. It naturally requires that the image of f is contained in the domain of g (except if f is global and surjective, see documentation of third parameter)
   @param perl::Object f A Morphism
   @param perl::Object g A Morphism, whose DOMAIN contains the image of f (see also next parameter)
   @param bool is_global_and_surjective Optional, false by default. If set to true, f is assumed to be a globally affine linear surjective map. In this case, g's DOMAIN needn't be the whole space. The composition will then be defined on the preimage of g's DOMAIN
   @return perl::Object A Morphism object, the composition "g after f"
   */
   perl::Object morphism_composition(perl::Object f, perl::Object g, bool is_global_and_surjective = false) {
      //First we homogenize to make sure everything is compatible
      f = f.CallPolymakeMethod("homogenize");
      g = g.CallPolymakeMethod("homogenize");
    
      //Extract values of f
      perl::Object f_domain = f.give("DOMAIN");
	Matrix<Rational> f_rays = f_domain.give("CMPLX_RAYS");
	Matrix<Rational> f_lin = f_domain.give("LINEALITY_SPACE");
	int f_ambient_dim = f_domain.give("CMPLX_AMBIENT_DIM");
	IncidenceMatrix<> f_cones = f_domain.give("CMPLX_MAXIMAL_CONES");
      bool is_global = f.give("IS_GLOBAL");
      Matrix<Rational> f_on_rays = f.give("RAY_VALUES");
      Matrix<Rational> f_on_lin = f.give("LIN_VALUES");
      Matrix<Rational> matrix; Vector<Rational> translate;
      if(is_global) {
	Matrix<Rational> m = f.give("MATRIX");
	Vector<Rational> t = f.give("TRANSLATE");
	matrix = m; translate = t;
      }
      
      //Extract values of g
      perl::Object g_domain = g.give("DOMAIN");
	Matrix<Rational> g_rays = g_domain.give("CMPLX_RAYS");
	Matrix<Rational> g_lin = g_domain.give("LINEALITY_SPACE");
	IncidenceMatrix<> g_cones = g_domain.give("CMPLX_MAXIMAL_CONES");
	int g_dim = g_domain.give("CMPLX_DIM");
      Matrix<Rational> g_on_rays = g.give("RAY_VALUES");
      Matrix<Rational> g_on_lin = g.give("LIN_VALUES");
      Vector< Matrix<Rational> > g_hreps_ineq;
      Vector< Matrix<Rational> > g_hreps_eq;
      
      //Prepare result variables
      Matrix<Rational> pullback_rays(0, f_ambient_dim);
      Matrix<Rational> pullback_lineality(0, f_ambient_dim);
	bool lineality_computed = false;
      Vector<Set<int> > pullback_cones;
      //The following two variables contain for each cone of the pullback domain the representation 
      //as an affine linear function on this cone
      Vector<Matrix<Rational> > pullback_matrices; 
      Vector<Vector<Rational> > pullback_translates;
      
      //Compute H-representations of all cones in g_domain
      solver<Rational> sv;
      for(int gcone = 0; gcone < g_cones.rows(); gcone++) {
	std::pair<Matrix<Rational>, Matrix<Rational> > p = sv.enumerate_facets(
	  zero_vector<Rational>() | g_rays.minor(g_cones.row(gcone),All), 
	  zero_vector<Rational>() | g_lin, true,false);
	g_hreps_ineq |= p.first;
	g_hreps_eq |= p.second;
      }//END compute fcone-H-rep
      
      //Now iterate all cones of f's domain
      for(int fcone = 0; fcone < f_cones.rows(); fcone++) {
	dbglog << "Computing on function cone " << fcone << endl;
	//Compute H-representation of the image of the cone: We have to convert the image values to homogeneous
	//coordinates
	Vector<Rational> homogenizer(f_cones.row(fcone).size());
	int homrow = 0;
	Set<int> fcone_set = f_cones.row(fcone);
	for(Entire<Set<int> >::iterator r = entire(fcone_set); !r.at_end(); r++) {
	  if(f_rays(*r,0) == 1) {
	    homogenizer[homrow] == 1; homrow++;
	  }
	}
	int image_dim = rank(f_on_rays.minor(f_cones.row(fcone),All) / f_on_lin);
	std::pair<Matrix<Rational>, Matrix<Rational> > image_rep = sv.enumerate_facets(
	      zero_vector<Rational>() | (homogenizer | f_on_rays.minor(f_cones.row(fcone),All)),
	      zero_vector<Rational>() | (zero_vector<Rational>() | f_on_lin), true, false); 
	
	//Compute representation of morphism on current cone
	Matrix<Rational> fmatrix;
	Vector<Rational> ftranslate;
	if(is_global) {
	    fmatrix = matrix;
	    ftranslate = translate;
	}
	else {
	    computeConeFunction(f_rays.minor(f_cones.row(fcone),All), f_lin,
				true, f_on_rays.minor(f_cones.row(fcone),All),
				f_on_lin, ftranslate, fmatrix);
	}
	
	dbgtrace << "Local representation is " << ftranslate << " and " << fmatrix << endl;
	
	//Iterate all cones of the function
	for(int gcone = 0; gcone < g_cones.rows(); gcone++) {
	  dbglog << "Intersecting with g cone " << gcone << endl;
	  //Compute intersection (or take the whole cone if the morphism is global and surjective)
	  Matrix<Rational> intersection_ineq;
	  Matrix<Rational> intersection_eq;
	  
	  if(is_global_and_surjective) {
	    intersection_ineq = g_hreps_ineq[gcone];
	    intersection_eq = g_hreps_eq[gcone];
	  }
	  else {
	    //TODO: Maybe just canonicalize h-rep and count equations?
	    intersection_ineq = image_rep.first / g_hreps_ineq[gcone];
	    intersection_eq = image_rep.second / g_hreps_eq[gcone];
	    
	    std::pair<Matrix<Rational>, Matrix<Rational> > p = sv.enumerate_vertices(
		image_rep.first / g_hreps_ineq[gcone],
		image_rep.second / g_hreps_eq[gcone],true,true);
	    
	    //Check dimension of intersection - if its not the correct one, take the next g cone
	    int interdim = rank(p.first.minor(All,~scalar2set(0))/ p.second.minor(All,~scalar2set(0))) - 1;
	    if( (is_global_and_surjective && interdim != g_dim) || 
		(!is_global_and_surjective && interdim != image_dim)) continue;
	    
	  }
	  
	  dbglog << "Computing representation on g cone" << endl;
	  
	  //Compute g's representation on the current cone
	  Vector<Rational> gtranslate;
	  Matrix<Rational> gmatrix;
	  computeConeFunction(g_rays.minor(g_cones.row(gcone),All), g_lin, true,
			      g_on_rays.minor(g_cones.row(gcone),All), g_on_lin, gtranslate, 
			      gmatrix);
	  
	  //Compute preimage of the intersection cone
	  //If (b,-A) is the representation of (in)equalities of the cone
	  // and x |-> v + Mx is the representation of the morphism, then
	  //(b - Av, -AM) is the representation of the preimage
	  //Mind the additional zero's we have to add for cddlib!
	  
	  dbglog << "Computing preimage and pullback" << endl;
	  
	  Set<int> homog_cols; homog_cols += 0; homog_cols += 1; //The first two columns are for homogenizing
	  Matrix<Rational> preimage_eq = intersection_eq.minor(All,~homog_cols) * fmatrix;
	    preimage_eq = zero_vector<Rational>() | 
			  ((intersection_eq.col(1) + intersection_eq.minor(All,~homog_cols) * ftranslate) 
			    | preimage_eq);
	  Matrix<Rational> preimage_ineq = intersection_ineq.minor(All,~homog_cols) * fmatrix;
	    preimage_ineq = zero_vector<Rational>() |
			  ((intersection_ineq.col(1) + intersection_ineq.minor(All,~homog_cols) * ftranslate)
			    | preimage_ineq);

	  std::pair<Matrix<Rational>, Matrix<Rational> > preimage_cone = sv.enumerate_vertices(
			    preimage_ineq, preimage_eq, true,true);
	  
	  dbgtrace << "Preimage has rays " << preimage_cone.first << " and lin " << preimage_cone.second << endl;
	
	  //Canonicalize rays and add the cone
// 	  if(!lineality_computed) {
// 	    pullback_lineality = 
// 	  }
	  
	  
	  //Now we compute the representation of h = g after f 
	  Matrix<Rational> hmatrix = gmatrix * fmatrix;
	  Vector<Rational> htranslate = gmatrix * ftranslate + gtranslate;
	  
	  
	  
	  
	  
	}//END iterate all function cones
	
	
	
	
      }//END iterate cones of morphism
      
      return perl::Object("Morphism");
      
  }//END morphism_composition
  
  void test(Matrix<Rational> rays, Matrix<Rational> lin, perl::Object morphism) {
    solver<Rational> sv;
    std::pair<Matrix<Rational>, Matrix<Rational> > eq = sv.enumerate_facets(zero_vector<Rational>() | rays,
									    zero_vector<Rational>() | lin,
									    true,false);
    
    pm::cout << "Ineq: " << eq.first << ", eq: " << eq.second << endl;
    
    Matrix<Rational> matrix = morphism.give("MATRIX");
    Vector<Rational> translate = morphism.give("TRANSLATE");
    
    pm::cout << "matrix " << matrix << endl;
    pm::cout << "translate " << translate << endl;
    
    Set<int> homog_cols; homog_cols += 0; homog_cols += 1;
    
    Matrix<Rational> newineq = eq.first.minor(All,~homog_cols) * matrix;
      newineq = (eq.first.col(1) + eq.first.minor(All,~homog_cols) * translate) | newineq;
      newineq = zero_vector<Rational>() | newineq;
      
    Matrix<Rational> neweq = eq.second.minor(All,~homog_cols) * matrix;
      neweq = (eq.second.col(1) + eq.second.minor(All,~homog_cols) * translate) | neweq;
      neweq = zero_vector<Rational>() | neweq;
    
    pm::cout << "New ineq " << newineq << endl;
    pm::cout << "new eq " << neweq << endl;
    
//     Matrix<Rational> transineq = 
    
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&test, "pitest(Matrix<Rational>, Matrix<Rational>,Morphism)");
  Function4perl(&morphism_composition, "morphism_composition(Morphism, Morphism; $=0)");
  
}}