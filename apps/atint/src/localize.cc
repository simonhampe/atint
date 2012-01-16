/*
 This *program is free software; you can redistribute it and/or
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
 
 This file provides convenience methods for creating locally restricted
 tropical varieties from given varieties
 */


#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/WeightedComplexRules.h"
#include "polymake/atint/refine.h"
#include "polymake/atint/specialvarieties.h"

namespace polymake { namespace atint { 
    
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object local_restrict(perl::Object complex, IncidenceMatrix<> cones) {
    //Extract values
    Vector<Set<int> > maximalCones = complex.give("MAXIMAL_CONES");
    Matrix<Rational> rays = complex.give("RAYS");
    Matrix<Rational> linspace = complex.give("LINEALITY_SPACE");
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    Vector<Integer> weights = complex.give("TROPICAL_WEIGHTS");
    
    //Find out which cones are no longe compatible
    Set<int> remainingCones;
    for(int c = 0; c < maximalCones.dim(); c++) {
      if(is_coneset_compatible(maximalCones[c], cones)) {
	remainingCones += c;
      }
    }
    
    //Adapt cone description and ray indices
    maximalCones = maximalCones.slice(remainingCones);
    weights = weights.slice(remainingCones);
    Set<int> usedRays = accumulate(maximalCones,operations::add());
    rays = rays.minor(usedRays,All);
    IncidenceMatrix<> newMaximalCones(maximalCones);
      newMaximalCones = newMaximalCones.minor(All,usedRays);
    cones = cones.minor(All,usedRays);
    
    perl::Object result("WeightedComplex");
      result.take("RAYS") << rays;
      result.take("MAXIMAL_CONES") << newMaximalCones;
      result.take("LINEALITY_SPACE") << linspace;
      result.take("USES_HOMOGENEOUS_C") << uses_homog;
      result.take("TROPICAL_WEIGHTS") << weights;
      result.take("LOCAL_RESTRICTION") << cones;
    
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object local_ray(perl::Object complex, int vertex) {
    //Convert vertex to incidence matrix
    Vector<Set<int> > matrix;
    Set<int> set;
      set += vertex;
      matrix |= set;
    return local_restrict(complex, IncidenceMatrix<>(matrix));
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object local_codim_1(perl::Object complex, int face) {
    //Convert codim face index to incidence matrix
    IncidenceMatrix<> codim = complex.give("CODIM_1_FACES");
    if(face >= codim.rows()) {
      throw std::runtime_error("Cannot localize at codim one face: Index is out of bounds.");
    }
    Vector<Set<int> > matrix;
      matrix |= codim.row(face);
    return local_restrict(complex, IncidenceMatrix<>(matrix));
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  perl::Object local_point(perl::Object complex, Vector<Rational> point) {
    //Normalize the vertex
    if(point.dim() <= 1) {
      throw std::runtime_error("Cannot localize at point: Point dimension is too low");
    }
    if(point[0] == 0) {
      throw std::runtime_error("Cannot localize at point: Point is not a vertex (or not given in homogeneous coordinates");
    }
    point /= point[0];
    
    //Homogenize the complex if necessary
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    if(!uses_homog) {
      complex = complex.CallPolymakeMethod("homogenize");
    }
    
    //First we refine the complex
    RefinementResult r = refinement(complex, point_variety(point),false,false,false,true);
    perl::Object refinedComplex = r.complex;
    
    //Then we look for the vertex
    Matrix<Rational> rays = refinedComplex.give("RAYS");
    Set<int> vertices = refinedComplex.give("VERTICES");
    int pointindex = -1;
    for(Entire<Set<int> >::iterator v = entire(vertices); !v.at_end(); v++) {
      if(rays.row(*v) == point) {
	pointindex = *v; break;
      }
    }
    
    //If we didn't find it, throw an error
    if(pointindex == -1) throw std::runtime_error("Cannot localize at point: Is not contained in support of complex.");
    
    //Otherwise localize
    return local_ray(refinedComplex,pointindex);
    
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  //Documentation see perl wrapper
  bool contains_point(perl::Object complex, Vector<Rational> point) {
    //Normalize the vertex
    if(point.dim() <= 1) {
      return false;
    }
    if(point[0] == 0) {
      return false;
    }
    point /= point[0];
    
    //Homogenize the complex if necessary
    bool uses_homog = complex.give("USES_HOMOGENEOUS_C");
    if(!uses_homog) {
      complex = complex.CallPolymakeMethod("homogenize");
    }
    
    //First we refine the complex
    RefinementResult r = refinement(complex, point_variety(point),false,false,false,true);
    perl::Object refinedComplex = r.complex;
    
    //Then we look for the vertex
    Matrix<Rational> rays = refinedComplex.give("RAYS");
    Set<int> vertices = refinedComplex.give("VERTICES");
    int pointindex = -1;
    for(Entire<Set<int> >::iterator v = entire(vertices); !v.at_end(); v++) {
      if(rays.row(*v) == point) {
	pointindex = *v; break;
      }
    }
    return pointindex != -1;
  }
  
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  UserFunction4perl("# @category Tropical geometry /Local geometry"
		    "# This takes a tropical variety and an IncidenceMatrix describing a set"
		    "# of cones (not necessarily maximal ones) of this variety. It will then"
		    "# create a variety that contains all compatible maximal cones and is"
		    "# locally restricted to the given cone set."
		    "# @param WeightedComplex complex An arbitrary weighted complex"
		    "# @param IncidenceMatrix cones A set of cones, indices refer to RAYS"
		    "# @return WeightedComplex The same complex, locally restricted to the given"
		    "# cones",
		    &local_restrict, "local_restrict(WeightedComplex,IncidenceMatrix)");
  
  UserFunction4perl("#@category Tropical geometry / Local geometry"
		    "# This takes a weighted complex and an index of one of its vertices or rays "
		    "# (the index is to be understood in RAYS)"
		    "# It then localizes the variety at this vertex / ray. The index should never"
		    "# correspond to a directional ray in a complex, since this would not be a cone"
		    "# @param WeightedComplex complex An arbitrary weighted complex"
		    "# @param Int ray The index of a ray/vertex in RAYS"
		    "# @return WeightedComplex The complex locally restricted to the given ray/vertex",
		    &local_ray, "local_ray(WeightedComplex,$)");
  
  UserFunction4perl("# @category Tropical geometry / Local geometry"
		    "# This takes a weighted complex and an index of one of its codimension one faces"
		    "# (The index is in CODIM_1_FACES) and computes the complex locally restricted"
		    "# to that face"
		    "# @param WeightedComplex complex An arbitrary weighted complex"
		    "# @param Int face An index of a face in CODIM_1_FACES"
		    "# @return WeightedComplex The complex locally restricted to the given face",
		    &local_codim_1, "local_codim_1(WeightedComplex,$)");
  
  UserFunction4perl("# @category Tropical geometry / Local geometry"
		    "# This takes a weighted complex and an arbitrary vertex in homogeneous "
		    "# coordinates that is supposed to lie in the support of the complex"
		    "# It then refines the complex such that the vertex is a cell in the polyhedral "
		    "# structure and returns the complex localized at this vertex"
		    "# @param WeightedComplex complex An arbitrary weighted complex"
		    "# @param Vector<Rational> v A vertex in homogeneous coordinates. It should lie"
		    "# in the support of the complex (otherwise an error is thrown)"
		    "# @return WeightedComplex The complex localized at the vertex",
		    &local_point, "local_point(WeightedComplex,Vector<Rational>)");
  
  UserFunction4perl("# @category Tropical geometry / Local geometry"
		    "# Takes a weighted complex and a point and computed whether that point lies in "
		    "# the complex"
		    "# @param perl::Object complex A weighted complex"
		    "# @param Vector<Rational> point An arbitrary vector in the same ambient"
		    "# dimension as complex"
		    "# @return bool Whether the point lies in the support of complex",
		    &contains_point,"contains_point(WeightedComplex,Vector<Rational>)");
		    
  
}}


