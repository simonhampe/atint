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
#include "polymake/atint/morphism_composition.h"


namespace polymake { namespace atint { 
    
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
   
   //Documentation see header
   perl::Object morphism_composition(perl::Object f, perl::Object g) {
      //First we homogenize to make sure everything is compatible
      f = f.CallPolymakeMethod("homogenize");
      g = g.CallPolymakeMethod("homogenize");
    
      // -------------------------- PREPARATIONS ----------------------------------- //
      
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
      Vector< Matrix<Rational> > f_hreps_ineq;
      Vector< Matrix<Rational> > f_hreps_eq;
      
      //Extract values of g
      perl::Object g_domain = g.give("DOMAIN");
	Matrix<Rational> g_rays = g_domain.give("CMPLX_RAYS");
	Matrix<Rational> g_lin = g_domain.give("LINEALITY_SPACE");
	IncidenceMatrix<> g_cones = g_domain.give("CMPLX_MAXIMAL_CONES");
// 	int g_dim = g_domain.give("CMPLX_DIM");
      Matrix<Rational> g_on_rays = g.give("RAY_VALUES");
      Matrix<Rational> g_on_lin = g.give("LIN_VALUES");
      Vector< Matrix<Rational> > g_hreps_ineq;
      Vector< Matrix<Rational> > g_hreps_eq;
      
      //Prepare result variables
      Matrix<Rational> pullback_rays(0, f_ambient_dim);
      Matrix<Rational> pullback_lineality(0, f_ambient_dim);
	bool lineality_computed = false;
      Vector<Set<int> > pullback_cones;
      Set<Set<int> > pullback_cones_set; //Used to check for doubles
      Matrix<Rational> pullback_ray_values(0,g_on_rays.cols());
      Matrix<Rational> pullback_lin_values(0,g_on_lin.cols());
      //The following two variables contain for each cone of the pullback domain the representation 
      //as an affine linear function on this cone
      Vector<Matrix<Rational> > pullback_matrices; 
      Vector<Vector<Rational> > pullback_translates;
      
      //Compute H-representations of all cones in f_domain and g_domain
      solver<Rational> sv;
      for(int fcone = 0; fcone < f_cones.rows(); fcone++) {
	std::pair<Matrix<Rational>, Matrix<Rational> > p = sv.enumerate_facets(
	  zero_vector<Rational>() | f_rays.minor(f_cones.row(fcone),All), 
	  zero_vector<Rational>() | f_lin, true,false);
	f_hreps_ineq |= p.first;
	f_hreps_eq |= p.second;
      }//END compute fcone-H-rep
      for(int gcone = 0; gcone < g_cones.rows(); gcone++) {
	std::pair<Matrix<Rational>, Matrix<Rational> > p = sv.enumerate_facets(
	  zero_vector<Rational>() | g_rays.minor(g_cones.row(gcone),All), 
	  zero_vector<Rational>() | g_lin, true,false);
	g_hreps_ineq |= p.first;
	g_hreps_eq |= p.second;
      }//END compute gcone-H-rep
      
      //Compute a homogenized version of the image of f
      Vector<Rational> homogenizer(f_rays.rows());
      for(int fr = 0; fr < f_rays.rows(); fr++) {
	if(f_rays(fr,0) == 1) {
	    homogenizer[fr] = 1;
	}
      }
      Matrix<Rational> f_on_rays_homog = homogenizer | f_on_rays;
      Matrix<Rational> f_on_lin_homog = zero_vector<Rational>() | f_on_lin;
      
      // --------------------------- COMPUTE GEOMETRY ---------------------------------- //
      
      //Now iterate all cones of f's domain
      for(int fcone = 0; fcone < f_cones.rows(); fcone++) {
	dbglog << "Computing on function cone " << fcone << endl;
	//Compute H-representation of the image of the cone
	std::pair<Matrix<Rational>, Matrix<Rational> > image_rep = sv.enumerate_facets(
	      zero_vector<Rational>() | (f_on_rays_homog.minor(f_cones.row(fcone),All)),
	      zero_vector<Rational>() | f_on_lin_homog, true, false); 
	int image_dim = image_rep.first.cols() - image_rep.second.rows() - 2;
	dbgtrace << "Image has dimension " << image_dim << endl;
	
	//Find out if f is global and surjective 
	bool is_global_and_surjective = is_global && image_dim == f_on_rays.cols();
	
	
	dbgtrace << "Image of cone has H-rep " << image_rep << endl;
	
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
	
// 	dbgtrace << "Local representation is " << ftranslate << " and " << fmatrix << endl;
	
	//Iterate all cones of the function
	for(int gcone = 0; gcone < g_cones.rows(); gcone++) {
	  dbglog << "Intersecting with g cone " << gcone << endl;
	  //Compute intersection (or take the whole cone if the morphism is global and surjective)
	  Matrix<Rational> intersection_ineq;
	  Matrix<Rational> intersection_eq;
	  
	  if(is_global_and_surjective) {
	    intersection_ineq = g_hreps_ineq[gcone];
	    intersection_eq = g_hreps_eq[gcone];
	    dbgtrace << "global, hence valid " << endl;
	  }
	  else {
	    //Compute an irredundant H-rep of the intersection
	    intersection_ineq = image_rep.first / g_hreps_ineq[gcone];
	    intersection_eq = image_rep.second / g_hreps_eq[gcone];
	    Matrix<Rational> isMatrix = intersection_ineq / intersection_eq;
	    
	    std::pair<Bitset,Bitset> isection = 
		  sv.canonicalize( intersection_ineq,intersection_eq,1);
	    
		  
	      
	    intersection_ineq = isMatrix.minor(isection.first,All);
	    intersection_eq = isMatrix.minor(isection.second,All);
		  
	    int interdim = isMatrix.cols()  - isection.second.size() - 2 ;
	    dbgtrace << "intersection dimension is " << interdim << endl;
	    
	    //Check dimension of intersection - if its not the correct one, take the next g cone
	    if( interdim != image_dim) continue;
	    
	    dbgtrace << "Cone is valid " << endl;	    
	  }
	  
	  dbgtrace << "Intersection cone has H-rep " << intersection_ineq << "\n" << intersection_eq << endl;
	  
	  dbglog << "Computing representation on g cone" << endl;
	  
	  //Compute g's representation on the current cone
	  Vector<Rational> gtranslate;
	  Matrix<Rational> gmatrix;
	  computeConeFunction(g_rays.minor(g_cones.row(gcone),All), g_lin, true,
			      g_on_rays.minor(g_cones.row(gcone),All), g_on_lin, gtranslate, 
			      gmatrix);

// 	  dbgtrace << "g's representation on this cone " << gmatrix << " and " << gtranslate << endl;
	  
	  //Compute preimage of the intersection cone
	  //If (b,-A) is the representation of (in)equalities of the cone
	  // and x |-> v + Mx is the representation of the morphism, then
	  //(b - Av, -AM) is the representation of the preimage
	  //Mind the additional zero's we have to add for cddlib!
	  
	  dbglog << "Computing preimage and pullback" << endl;
	  
	  Set<int> homog_cols; homog_cols += 0; homog_cols += 1; //The first two columns are for homogenizing
	  Matrix<Rational> preimage_ineq = intersection_ineq.minor(All,~homog_cols) * fmatrix;
	    preimage_ineq = zero_vector<Rational>() |
			  ((intersection_ineq.col(1) + intersection_ineq.minor(All,~homog_cols) * ftranslate)
			    | preimage_ineq);
		
	  Matrix<Rational> preimage_eq(0,preimage_ineq.cols());
	  if(intersection_eq.rows() > 0) { //For the equalities consider the special case that there are none
	      preimage_eq = intersection_eq.minor(All,~homog_cols) * fmatrix;
	      preimage_eq = zero_vector<Rational>() | 
			    ((intersection_eq.col(1) + intersection_eq.minor(All,~homog_cols) * ftranslate) 
			    | preimage_eq);
	  }
	  
	  //Intersect with the fcone
	  preimage_ineq /= f_hreps_ineq[fcone];
	  preimage_eq /= f_hreps_eq[fcone];
			  

	  dbgtrace << "Preimage ineq " << preimage_ineq << "\n eq " << preimage_eq << endl;
	  
			  
	  std::pair<Matrix<Rational>, Matrix<Rational> > preimage_cone = sv.enumerate_vertices(
			    preimage_ineq, preimage_eq, true,true);
	  
// 	  dbgtrace << "Preimage has rays " << preimage_cone.first << " and lin " << preimage_cone.second << endl;

	  //Dehomogenize
	  Matrix<Rational> preimage_rays = preimage_cone.first.minor(All,~scalar2set(0));
	  Matrix<Rational> preimage_lin = preimage_cone.second.minor(All,~scalar2set(0));
	  
	  dbgtrace << "Canonicalizing rays" << endl;
	  
	  //Canonicalize rays and create cone
	  if(!lineality_computed) {
	    pullback_lineality = preimage_lin;
	    lineality_computed = true;
	    dbgtrace << "Setting lineality to " << pullback_lineality << endl;
	  }
	  Set<int> pcone; 
	  for(int r = 0; r < preimage_rays.rows(); r++) {
	    //Canonicalize ray
	    if(preimage_rays(r,0) != 0) {
	      preimage_rays.row(r) *= (1/preimage_rays(r,0));
	    }
	    else {
	      for(int c = 1; c < preimage_rays.cols(); c++) {
		if(preimage_rays(r,c) != 0) {
		  preimage_rays.row(r) *= (1/abs(preimage_rays(r,c)));
		  break;
		}
	      }
	    }
	    //Find correct ray index
	    int ray_index = -1;
	    for(int oray = 0; oray < pullback_rays.rows(); oray++) {
	      if(pullback_rays.row(oray) == preimage_rays.row(r)) {
		ray_index = oray; 
		break;
	      }
	    }
	    //Insert ray if necessary and add index to set
	    if(ray_index == -1) {
	      pullback_rays /= preimage_rays.row(r);
	      ray_index = pullback_rays.rows() -1;
	    }
	    pcone += ray_index;
	  }
	  dbgtrace << "Ray set is " << pcone << endl;
	  //Add cone if it doesn't exist yet
	  if(!pullback_cones_set.contains(pcone)) {
	    pullback_cones |= pcone;
	    pullback_cones_set += pcone;
	    dbgtrace << "Adding cone " << pcone << endl;
	    
	    //Now we compute the representation of h = g after f 
	    Matrix<Rational> hmatrix = gmatrix * fmatrix;
	    Vector<Rational> htranslate = gmatrix * ftranslate + gtranslate;
	    
	    dbgtrace << "Composition on preimage: " << hmatrix << " and " << htranslate << endl;
	    
	    pullback_matrices |= hmatrix;
	    pullback_translates |= htranslate;
	    
	  }
	  dbgtrace << "Rays now read " << pullback_rays << endl;
	  
	  
	}//END iterate all function cones
	
      }//END iterate cones of morphism
      
      // ------------------------- COMPUTE VALUES --------------------------------- //
      
      dbglog << "Computing values on CMPLX_RAYS " << endl;
      
      //Compute CMPLX_RAYS / CMPLX_MAXIMAL_CONES
      
      perl::Object pullback_domain("WeightedComplex");
	pullback_domain.take("RAYS") << pullback_rays;
	pullback_domain.take("MAXIMAL_CONES") << pullback_cones;
	pullback_domain.take("LINEALITY_SPACE") << pullback_lineality;
	pullback_domain.take("USES_HOMOGENEOUS_C") << true;
	
      Matrix<Rational> pb_cmplx_rays = pullback_domain.give("CMPLX_RAYS");
	Matrix<Rational> pb_crays_dehomog = pb_cmplx_rays.minor(All,~scalar2set(0));
      IncidenceMatrix<> pb_cmplx_cones = pullback_domain.give("CMPLX_MAXIMAL_CONES");
	IncidenceMatrix<> pb_cones_by_rays = T(pb_cmplx_cones);
      
      //Go trough all rays
      int basepoint = -1; //Save the first vertex
      for(int cr = 0; cr < pb_cmplx_rays.rows(); cr++) {
	//Take any cone containing this ray
	int cone_index = *(pb_cones_by_rays.row(cr).begin());
	Matrix<Rational> cone_matrix = pullback_matrices[cone_index];
	Vector<Rational> cone_translate = pullback_translates[cone_index];
	//If its a vertex, just compute the value
	if(pb_cmplx_rays(cr,0) == 1) {
	    pullback_ray_values /= 
		(cone_matrix* pb_crays_dehomog.row(cr) + cone_translate);
	    if(basepoint == -1) basepoint = cr;
	}
	//Otherwise find an associated vertex
	else {
	    Set<int> rays_in_cone = pb_cmplx_cones.row(cone_index);
	    for(Entire<Set<int> >::iterator aRay = entire(rays_in_cone); !aRay.at_end(); aRay++) {
	      if(pb_cmplx_rays(*aRay,0) == 1) {
		Vector<Rational> sum = pb_crays_dehomog.row(*aRay) + pb_crays_dehomog.row(cr);
		pullback_ray_values /= 
		  (
		  (cone_matrix* sum + cone_translate) - 
		  (cone_matrix* pb_crays_dehomog.row(*aRay) + cone_translate)
		  );
		  break;
	      }
	    }
	}
      }
      
      //Now compute lineality values
      Matrix<Rational> pb_lin_dehomog = pullback_lineality.minor(All,~scalar2set(0));
      for(int l = 0; l < pullback_lineality.rows(); l++) {
	Vector<Rational> sum = pb_crays_dehomog.row(basepoint) + pb_lin_dehomog.row(l);
	pullback_lin_values /= 
	  (
	  (pullback_matrices[0] * sum + pullback_translates[0]) - 
	  pullback_ray_values.row(basepoint)
	  );
      }
      
      // ------------------------------- RETURN RESULT ------------------------------- //
      
      dbglog << "Done. Preparing result" << endl;
      
      perl::Object result("Morphism");
	result.take("DOMAIN") << pullback_domain;
	result.take("RAY_VALUES") << pullback_ray_values;
	result.take("LIN_VALUES") << pullback_lin_values;
      
      return result;
      
  }//END morphism_composition
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&morphism_composition, "morphism_composition(Morphism, Morphism)");
  
}}