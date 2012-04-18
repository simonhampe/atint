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
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/cdd_helper_functions.h"
#include "polymake/atint/WeightedComplexRules.h"

namespace polymake { namespace atint { 
  
  using polymake::polytope::cdd_interface::solver;
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;

  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  void cdd_normalize_rays(Matrix<Rational> &rays, bool uses_homog) {
    for(int r = 0; r < rays.rows(); r++) {
      //Vertices are normalized to have first coordinate 1
      if(rays(r,0) != 0 && uses_homog) {
	rays.row(r) *= (1/rays(r,0));
      }
      //Rays are normalized to have first nonzero coordinate +-1
      else {
	for(int c = 0; c < rays.cols(); c++) {
	  if(rays(r,c) != 0) {
	    rays.row(r) *= abs(1/rays(r,c));  
	    break;
	  }
	}
      }
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  std::pair<Matrix<Rational>, Matrix<Rational> > cdd_cone_intersection(
      const Matrix<Rational> &xrays, const Matrix<Rational> &xlin, 
      const Matrix<Rational> &yrays, const Matrix<Rational> &ylin, bool uses_homog) {
    
    solver<Rational> sv;
    
	
    //Compute facets
    std::pair<Matrix<Rational>, Matrix<Rational> > x_eq =
      sv.enumerate_facets(zero_vector<Rational>() | xrays, zero_vector<Rational>() | xlin,true,false);
    std::pair<Matrix<Rational>, Matrix<Rational> > y_eq =
      sv.enumerate_facets(zero_vector<Rational>() | yrays, zero_vector<Rational>() | ylin,true,false);
    
    //Compute intersection rays
    std::pair<Matrix<Rational>, Matrix<Rational> > inter =
      sv.enumerate_vertices(x_eq.first / y_eq.first, x_eq.second / y_eq.second, true,true);
      
    //Truncate and normalize
    inter.first = inter.first.minor(All,~scalar2set(0));
    inter.second = inter.second.minor(All,~scalar2set(0));
    cdd_normalize_rays(inter.first,uses_homog);
    
    return inter;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see header
  fan_intersection_result cdd_fan_intersection(	
      Matrix<Rational> xrays, Matrix<Rational> xlin, IncidenceMatrix<> xcones,
      Matrix<Rational> yrays, Matrix<Rational> ylin, IncidenceMatrix<> ycones,
      bool uses_homog) {
      
    solver<Rational> sv;
    
    //Precompute h-representations of the x-cones and y-cones
    Vector<std::pair<Matrix<Rational>, Matrix<Rational> > > xequations;
    for(int xc = 0; xc < xcones.rows(); xc++) {
      xequations |= sv.enumerate_facets(zero_vector<Rational>() | xrays.minor(xcones.row(xc),All),
					zero_vector<Rational>() | xlin,true,false);
    }
    Vector<std::pair<Matrix<Rational>, Matrix<Rational> > > yequations;
    for(int yc = 0; yc < ycones.rows(); yc++) {
      yequations |= sv.enumerate_facets(zero_vector<Rational>() | yrays.minor(ycones.row(yc),All),
					zero_vector<Rational>() | ylin,true,false);
    }
    
    //Now compute intersections
    Matrix<Rational> interrays;
    Matrix<Rational> interlineality;
      bool lineality_computed = false;
    Vector<Set<int> > intercones;
    
    Vector<Set<int> > xcontainers;
    Vector<Set<int> > ycontainers;
    
    for(int xc = 0; xc < xcones.rows(); xc++) {
      for(int yc = 0; yc < ycones.rows(); yc++) {	
	//Compute intersection
	std::pair<Matrix<Rational>, Matrix<Rational> > inter = 
	  sv.enumerate_vertices(xequations[xc].first / yequations[yc].first,
				xequations[xc].second / yequations[yc].second,
				true,true);
	  
	if(!lineality_computed) {
	    interlineality = inter.second.rows() > 0 ? 
			  inter.second.minor(All,~scalar2set(0)) :
			  Matrix<Rational>();
	    lineality_computed = true;
	}
	
	
	//Truncate and normalize
	inter.first = inter.first.minor(All,~scalar2set(0));
	
	//The empty cone will not be included, if uses_homog = true
	if(uses_homog && inter.first.rows() == 0) {
	  continue;
	}
	
	//If we are in homog. coordinates and the cone contains no vertices (i.e. the intersection is actually 
	//empty), we leave it out
	if(uses_homog) {
	  if(inter.first.col(0) == zero_vector<Rational>(inter.first.rows())) continue;
	}
	
	cdd_normalize_rays(inter.first,uses_homog);
	
	//Insert rays into ray list and create cone
	Set<int> new_cone_set;
	bool new_ray_added = false;
	for(int r = 0; r < inter.first.rows(); r++) {
	    int ray_index = -1;
	    for(int oray = 0; oray < interrays.rows(); oray++) {
	      if(inter.first.row(r) == interrays.row(oray)) {
		ray_index = oray; 
		new_cone_set += ray_index;
		break;
	      }
	    }
	    if(ray_index == -1) {
	      interrays /= inter.first.row(r);
	      new_cone_set += (interrays.rows()-1);
	      new_ray_added = true;
	    }
	}
	
	//Make sure we don't add a cone twice
	//Also: Remember intersections that are contained in this one or contain this one
	Set<int> containedCones;
	Set<int> containerCones;
	int new_cone_index = -1;
	for(int j = 0; j < intercones.dim(); j++) {
	    if(intercones[j] == new_cone_set) {
	      new_cone_index = j; 
	    }
	    else {
	      int sz = (intercones[j] * new_cone_set).size();
	      if(sz == intercones[j].size()) containedCones += j;
	      if(sz == new_cone_set.size()) containerCones += j;
	    }
	}
	if(new_cone_index == -1) {
	    intercones |= new_cone_set;
	    new_cone_index = intercones.dim()-1;
	    xcontainers |= Set<int>();
	    ycontainers |= Set<int>();
	}
	
	//First add all containers from the containing cones
	for(Entire<Set<int> >::iterator sup = entire(containerCones); !sup.at_end(); sup++) {
	  xcontainers[new_cone_index] += xcontainers[*sup];
	  ycontainers[new_cone_index] += ycontainers[*sup];
	}
	//Add xc and yc as containers
	xcontainers[new_cone_index] += xc;
	ycontainers[new_cone_index] += yc;
	//Add all current containers to the contained cones
	for(Entire<Set<int> >::iterator sub = entire(containedCones); !sub.at_end(); sub++) {
	    xcontainers[*sub] += xcontainers[new_cone_index];
	    ycontainers[*sub] += ycontainers[new_cone_index];
	}
	
	
	
      }//END iterate ycones
    }//END iterate xcones
    
    //Create result:
    fan_intersection_result f;
      f.rays = interrays;
	if(interlineality.rows() == 0) interlineality = Matrix<Rational>(0,interrays.cols());
      f.lineality_space = interlineality;
      f.cones = IncidenceMatrix<>(intercones);
      f.xcontainers = IncidenceMatrix<>(xcontainers);
      f.ycontainers = IncidenceMatrix<>(ycontainers);
    return f;
    
  }
  
  perl::Object test_fan_intersection(perl::Object X, perl::Object Y) {
    Matrix<Rational> xrays = X.give("RAYS");
    Matrix<Rational> xlin = X.give("LINEALITY_SPACE");
    IncidenceMatrix<> xcones = X.give("MAXIMAL_CONES");
    bool x_homog = X.give("USES_HOMOGENEOUS_C");
    
    Matrix<Rational> yrays = Y.give("RAYS");
    Matrix<Rational> ylin = Y.give("LINEALITY_SPACE");
    IncidenceMatrix<> ycones = Y.give("MAXIMAL_CONES");
    bool y_homog = Y.give("USES_HOMOGENEOUS_C");
    
    fan_intersection_result f = cdd_fan_intersection(xrays,xlin, xcones, yrays,ylin,ycones, x_homog || y_homog);
    
    perl::Object p("WeightedComplex");
      p.take("RAYS") << f.rays;
      p.take("LINEALITY_SPACE") << f.lineality_space;
      p.take("MAXIMAL_CONES") << f.cones;
      
      pm::cout << "X-Containers: " << f.xcontainers << endl;
      pm::cout << "Y-Containers: " << f.ycontainers << endl;
      
      return p;
  }
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  Function4perl(&cdd_cone_intersection, "cdd_cone_intersection(Matrix<Rational>,Matrix<Rational>,Matrix<Rational>,Matrix<Rational>,$)");
  
  Function4perl(&cdd_normalize_rays, "cdd_normalize_rays(Matrix<Rational>, $)");
  
  Function4perl(&test_fan_intersection, "tfi(WeightedComplex, WeightedComplex)");
  
}}