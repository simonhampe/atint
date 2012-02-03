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
 * This file contains the functions used to compute visual data
 * for a WeightedComplex
 */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  using polymake::polytope::cdd_interface::solver;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief Computes the data necessary to visualize the polyhedral complex. Will be called by CMPLX_VISUAL.
   @param fan The weighted polyhedral complex (in homogeneous coordinates!)
   @param scale The scaling factor for the directional rays that are added to the affine rays
   @param showWeights Whether the weights of the cells should be displayed in the center
   @return A perl array containg:
   A list of rational polytopes representing the complex.
   For each cell of the complex, there is a polytope obtained by adding all directional rays (with a scaling factor) to all
   affine rays.
   2)  A polytope::PointConfiguration that will contain the center of each cell as vertex, labelled with the corresponding weight. This is only computed if showWeights is true.
  */
  perl::ListReturn computeVisualPolyhedra(const perl::Object &fan, const Rational &scale, bool showWeights) { 
    //Extract values
    bool weightsExist = fan.exists("TROPICAL_WEIGHTS");
    Array<Integer> weights;
      if(weightsExist) {
	weights = fan.give("TROPICAL_WEIGHTS");
      }
    int ambient_dim = fan.give("FAN_AMBIENT_DIM");
    Matrix<Rational> rays = fan.give("RAYS");
    Matrix<Rational> linealitySpace = fan.give("LINEALITY_SPACE");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    
    dbgtrace << "Extracted values" << endl;
    
    //First separate affine and directional rays
    Set<int> affineRays;
    Set<int> directionalRays;
    for(int r = 0; r < rays.rows(); r++) {
      if(rays.row(r)[0] == 0) {
	directionalRays = directionalRays + r;
      }
      else {
	affineRays = affineRays + r;
      }
    }
    
    dbgtrace << "Separated rays" << endl;
    
    //Create a polytope for each cone
    perl::ListReturn result;
    
    //This will contain the cell centers with the weight labels
    perl::Object weightCenters("polytope::PointConfiguration");
    Matrix<Rational> centermatrix(0,ambient_dim);
    Vector<std::string> centerlabels;
      
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      
      dbgtrace << "Computing geometry for cone " << mc << endl;
      
      //First, create a point matrix by adding all directional rays to all affine rays
      Matrix<Rational> v(0,ambient_dim);
      
      Set<int> maxAffine = maximalCones.row(mc) * affineRays;
      Set<int> maxDirectional = maximalCones.row(mc) * directionalRays;
      
      for(Entire<Set<int> >::iterator aRay = entire(maxAffine); !aRay.at_end(); ++aRay) {
	v = v / rays.row(*aRay);
	for(Entire<Set<int> >::iterator dRay = entire(maxDirectional); !dRay.at_end(); ++dRay) {
	    v = v / (rays.row(*aRay) + scale * rays.row(*dRay));
	}
	for(int linrow = 0; linrow < linealitySpace.rows(); linrow++) {
	    v = v / (rays.row(*aRay) + scale * linealitySpace.row(linrow));
	    v = v / (rays.row(*aRay) - scale * linealitySpace.row(linrow));
	}
      }
      
      //Then create a rational polytope for labelling
      dbgtrace << "Points " << v << endl;
      perl::Object ratPolytope("polytope::Polytope<Rational>");
	ratPolytope.take("POINTS") << v;
      
      result << ratPolytope;
            
      //If necessary, compute centers and label with weights
      if(showWeights && weightsExist) {
	Vector<Rational> center = ratPolytope.give("VERTEX_BARYCENTER");
	centermatrix = centermatrix / center;
	std::ostringstream wlabel;
	wlabel << "# " << mc << ": " << weights[mc];
	centerlabels = centerlabels | wlabel.str();
      }      
      
      
      dbgtrace << "Done." << endl;
      
    }
    
    if(showWeights && weightsExist) {
      dbgtrace << "Computed weight labels, inserting them" << endl;
      weightCenters.take("POINTS") << centermatrix;
      weightCenters.take("LABELS") << centerlabels;
      dbgtrace << "Done" << endl;
    }
    
    result << weightCenters;
     
    dbgtrace << "Done." << endl;
      
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Computes the polyhedral data necessary for visualization with a bounding box
    @param perl::Object fan The polyhedral complex to be visualized, in homogeneous coordinates
    @param bool isRelative true, iff the bounding box is given relative to the complex
    @param bool showWeights If true, the barycenters of the polytopes are computed for weight labelling
    @param Rational bbDistance The relative distance of the border of the bounding box to the affine part of the complex (should be a positive number).
    @param bool onlyBoundingBox If true, only the relative bounding box with respect to the given bbDistance is computed and returned
    @param Matrix<Rational> bBox The absolute bounding box needed if isRelative is false. Given by two row vectors indicating the extreme points of the box
    @return A perl::ListReturn containing 
    1) the list of polytopes to be rendered
    2) A polytope::PointConfiguration that will contain the center of each cell as vertex, labelled with the corresponding weight. This is only computed if showWeights is true, but is contained in the ListReturn in any case.
    If however, onlyBoundingBox is true, the ListReturn will only contain a Matrix<Rational> specifying the relative  bounding box.
  */
  perl::ListReturn computeBoundedVisual(perl::Object fan, bool isRelative, bool showWeights,Rational bbDistance, bool onlyBoundingBox, Matrix<Rational> bBox) {
    //Extract values
    int ambient_dim = fan.give("CMPLX_AMBIENT_DIM");
    Matrix<Rational> rays = fan.give("RAYS");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    Matrix<Rational> facetNormals = fan.give("FACET_NORMALS");
    Matrix<Rational> facetNormalsInCones = fan.give("MAXIMAL_CONES_FACETS");
    Matrix<Rational> linearSpan = fan.give("LINEAR_SPAN_NORMALS");
    IncidenceMatrix<> linearSpanInCones = fan.give("MAXIMAL_CONES_LINEAR_SPAN_NORMALS");
    int fan_dim = fan.give("CMPLX_DIM");
    
     //First separate affine and directional rays
    Set<int> affineRays;
    Set<int> directionalRays;
    for(int r = 0; r < rays.rows(); r++) {
      if(rays.row(r)[0] == 0) {
	directionalRays = directionalRays + r;
      }
      else {
	affineRays = affineRays + r;
      }
    }  
    
    dbgtrace << "Computing bounding box..." << endl;
    
    //Compute facets of the bounding box
    Matrix<Rational> bbFacets(0,ambient_dim+1);
    
    //For each coordinate, contain minimum and maximum
    Vector<Rational> minCoord(ambient_dim);
    Vector<Rational> maxCoord(ambient_dim);
    //If bounding mode is relative, determine the maximal/minimal coordinates of the affine rays
    if(isRelative) {
      for(Entire<Set<int> >::iterator aff = entire(affineRays); !aff.at_end(); aff++) {
	for(int i = 0; i < ambient_dim; i++) {
	  Rational val = rays(*aff,i+1);
	  if(val > maxCoord[i]) maxCoord[i] = val;
	  if(val < minCoord[i]) minCoord[i] = val;
	}
      }
      //Now add the bbDistance to all values
      for(int i = 0; i < ambient_dim; i++) {
	maxCoord[i] += bbDistance;
	minCoord[i] -= bbDistance;
      }
      if(onlyBoundingBox) {
	Matrix<Rational> bb(0,ambient_dim);
	bb /= minCoord; 
	bb /= maxCoord;
	perl::ListReturn smallResult;
	  smallResult << bb;
	return smallResult;
      }
    }
    //otherwise take min and max from the given bounding box
    else {
      for(int i = 0; i < ambient_dim; i++) {
	maxCoord[i] = bBox(0,i) > bBox(1,i)? bBox(0,i) : bBox(1,i);
	minCoord[i] = bBox(0,i) < bBox(1,i)? bBox(0,i) : bBox(1,i);
      }
    }
    //Now make these coordinates into facets
    for(int i = 0; i < ambient_dim; i++) {
      Vector<Rational> facetVector = unit_vector<Rational>(ambient_dim,i);
      bbFacets /= (maxCoord[i] | -facetVector);
      bbFacets /= (-minCoord[i] | facetVector);
    }
    
    dbgtrace << "Done." << endl;
    
    perl::ListReturn result;
    
    //This will contain the cell centers with the weight labels
    perl::Object weightCenters("polytope::PointConfiguration");
    Matrix<Rational> centermatrix(0,ambient_dim);
    Vector<std::string> centerlabels;
    Array<Integer> weights;
    if(showWeights) {
      weights = fan.give("TROPICAL_WEIGHTS");
    }
    
    //Now compute all polyhedra to be rendered
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      dbgtrace << "Computing polytope of cone " << mc << endl;
      //Compute the facets ans equalities of the current cone and add the bbox facets
      Matrix<Rational> facets(0,ambient_dim+1);
      Matrix<Rational> linspan = linearSpan.minor(linearSpanInCones.row(mc),All);
	linspan = linspan;
      for(int fn = 0; fn < facetNormalsInCones.cols(); fn++) {
	if(facetNormalsInCones(mc,fn) == 1) {
	    facets /= facetNormals.row(fn);
	}
	if(facetNormalsInCones(mc,fn) == -1) {
	    facets /= (-facetNormals.row(fn));
	}
      }
      facets /= bbFacets;
      //facets = facets;
      
      dbgtrace << "Facets are " << facets << "Equalities are " << linspan << endl;
      
      //Compute the polytope vertices from that
      Matrix<Rational> polyRays = solver<Rational>().enumerate_vertices(zero_vector<Rational>()| facets, zero_vector<Rational>() | linspan, true,true).first;
      polyRays = polyRays.minor(All,~scalar2set(0));
      //Normalize
      for(int r = 0; r < polyRays.rows(); r++) {
	if(polyRays(r,0) != 0) polyRays.row(r) /= polyRays(r,0);
      }
      //We have to make sure that the polytope has
      //at least dim +1 vertices after cutting, otherwise its a point set or graph to the
      //visualization and all the Facet options don't work
      if(polyRays.rows() >= fan_dim+1) {
	perl::Object polytope("polytope::Polytope<Rational>");
	  polytope.take("VERTICES") << polyRays; //The polytope shouldn't have a lineality space
	result << polytope;
	
	//If weight labels should be displayed, compute the vertex barycenter of the polytope and
	// label it
	if(showWeights) {
	  Vector<Rational> barycenter = average(rows(polyRays));
	    //barycenter /= barycenter[0];
	  centermatrix = centermatrix / barycenter;
	  std::ostringstream wlabel;
	  wlabel << "# " << mc << ": " << weights[mc];
	  centerlabels = centerlabels | wlabel.str();
	}
      }
//       }
//       catch(...) { //An error should only occur if the polytope is empty. Then just omit it
// 	dbgtrace << "Cone " << mc << " not in bounding box. Omitting." << endl;
// 	pm::cout << "Test" << endl;
//       }
      
    }//END iterate rendered polyhedra
    
    if(showWeights) {
      weightCenters.take("POINTS") << centermatrix;
      weightCenters.take("LABELS") << centerlabels;
    }
    result << weightCenters;
    
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  
  Function4perl(&computeVisualPolyhedra, "computeVisualPolyhedra(WeightedComplex, Rational, $)");

  Function4perl(&computeBoundedVisual, "computeBoundedVisual(WeightedComplex, $, $, Rational,$, Matrix<Rational>)");
  
}}