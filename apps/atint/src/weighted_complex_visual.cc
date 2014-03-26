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
    
    //dbgtrace << "Extracted values" << endl;
    
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
    
    //dbgtrace << "Separated rays" << endl;
    
    //Create a polytope for each cone
    perl::ListReturn result;
    
    //This will contain the cell centers with the weight labels
    perl::Object weightCenters("polytope::PointConfiguration");
    Matrix<Rational> centermatrix(0,ambient_dim);
    Vector<std::string> centerlabels;
      
    for(int mc = 0; mc < maximalCones.rows(); mc++) {
      
      //dbgtrace << "Computing geometry for cone " << mc << endl;
      
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
      //dbgtrace << "Points " << v << endl;
      perl::Object ratPolytope("polytope::Polytope<Rational>");
	ratPolytope.take("POINTS") << v;
      
      result << ratPolytope;
            
      //If necessary, compute centers and label with weights
      if(showWeights && weightsExist) {
	Vector<Rational> center = ratPolytope.give("VERTEX_BARYCENTER");
	centermatrix = centermatrix / center;
	std::ostringstream wlabel;
	wlabel << "# " << mc << ": " << weights[mc];
	centerlabels |= wlabel.str();
      }      
      
      
      //dbgtrace << "Done." << endl;
      
    }
    
    if(showWeights && weightsExist) {
      //dbgtrace << "Computed weight labels, inserting them" << endl;
      weightCenters.take("POINTS") << centermatrix;
      weightCenters.take("LABELS") << centerlabels;
      //dbgtrace << "Done" << endl;
    }
    
    result << weightCenters;
     
    //dbgtrace << "Done." << endl;
      
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
   @brief This takes a list of rays (as row vectors) and computes the corresponding bounding box, i.e. it computes the minima /maxima over all coordinates, subtracts/adds a given distance and returns the resulting 2x(no of coordinates)-matrix. The first row contains the min-coords, the second the max-coords
   */
  Matrix<Rational> boundingBox(Matrix<Rational> rays, Rational distance) {
    if(rays.rows() == 0) return Matrix<Rational>(2,rays.cols());
    Vector<Rational> min_values = rays.row(0);
    Vector<Rational> max_values = rays.row(0);
    for(int r = 1; r < rays.rows(); r++) {
      for(int c = 0; c < rays.cols(); c++) {
	if(rays(r,c) < min_values[c]) min_values[c] = rays(r,c);
	if(rays(r,c) > max_values[c]) max_values[c] = rays(r,c);
      }
    }
    //Add distance
    min_values -= distance* ones_vector<Rational>(rays.cols());
    max_values += distance* ones_vector<Rational>(rays.cols());
    return min_values / max_values;
  }
  

  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  /**
    @brief Computes the polyhedral data necessary for visualization with a bounding box
    @param perl::Object fan The polyhedral complex to be visualized, in homogeneous coordinates
    @param bool isRelative true, iff the bounding box is given relative to the complex
    @param bool showWeights If true, the barycenters of the polytopes are computed for weight labelling
    @param Rational bbDistance The relative distance of the border of the bounding box to the affine part of the complex (should be a positive number).
    @param Matrix<Rational> bBox The absolute bounding box needed if isRelative is false. Given by two row vectors indicating the extreme points of the box
    @param Array<String> clabels If showWeights is false and this array has positive length, these strings will be used to label the maximal cones (missing labels are replaced by the emtpy string)
    @return A perl::ListReturn containing 
    1) the list of polytopes to be rendered
    2) A polytope::PointConfiguration that will contain the center of each cell as vertex, labelled with the corresponding weight. This is only computed if showWeights is true, but is contained in the ListReturn in any case.
  */
  perl::ListReturn computeBoundedVisual(perl::Object fan, bool isRelative, bool showWeights,Rational bbDistance, Matrix<Rational> bBox, Array<std::string> clabels) {
    //Extract values
    int ambient_dim = fan.give("CMPLX_AMBIENT_DIM");
    Matrix<Rational> rays = fan.give("RAYS");
    Set<int> vertices = fan.give("VERTICES");
    IncidenceMatrix<> maximalCones = fan.give("MAXIMAL_CONES");
    Matrix<Rational> facetNormals = fan.give("FACET_NORMALS");
    Matrix<Rational> facetNormalsInCones = fan.give("MAXIMAL_CONES_FACETS");
    Matrix<Rational> linearSpan = fan.give("LINEAR_SPAN_NORMALS");
    IncidenceMatrix<> linearSpanInCones = fan.give("MAXIMAL_CONES_LINEAR_SPAN_NORMALS");
    int fan_dim = fan.give("CMPLX_DIM");
    
    bool use_labels = !showWeights && clabels.size() > 0;
    
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
    
    //dbgtrace << "Computing bounding box..." << endl;
    
    //Compute facets of the bounding box
    Matrix<Rational> bbFacets(0,ambient_dim+1);
    
    //For each coordinate, contain minimum and maximum
    Vector<Rational> minCoord(ambient_dim);
    Vector<Rational> maxCoord(ambient_dim);
    //If bounding mode is relative, determine the maximal/minimal coordinates of the affine rays
    if(isRelative) {
      if(vertices.size() == 0) {
	minCoord = - bbDistance * ones_vector<Rational>(ambient_dim);
	maxCoord = bbDistance * ones_vector<Rational>(ambient_dim);
      }
      else {
	Matrix<Rational> bMatrix = boundingBox(rays.minor(vertices,~scalar2set(0)), bbDistance);
	minCoord = bMatrix.row(0);
	maxCoord = bMatrix.row(1);
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
    
    //dbgtrace << "Done." << endl;
    
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
      //dbgtrace << "Computing polytope of cone " << mc << endl;
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
      
      //dbgtrace << "Facets are " << facets << "Equalities are " << linspan << endl;
      
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
	if(showWeights || use_labels) {
	  Vector<Rational> barycenter = average(rows(polyRays));
	    //barycenter /= barycenter[0];
	  centermatrix = centermatrix / barycenter;
	  std::ostringstream wlabel;
	  wlabel << "# " << mc << ": ";
	  if(showWeights) wlabel << weights[mc];
	  else {
	    if(mc < clabels.size()) wlabel << clabels[mc];
	  }
	  centerlabels |= wlabel.str();
	}
      }
//       }
//       catch(...) { //An error should only occur if the polytope is empty. Then just omit it
// 	dbgtrace << "Cone " << mc << " not in bounding box. Omitting." << endl;
// 	pm::cout << "Test" << endl;
//       }
      
    }//END iterate rendered polyhedra
    
    if(showWeights || use_labels) {
      weightCenters.take("POINTS") << centermatrix;
      weightCenters.take("LABELS") << centerlabels;
    }
    result << weightCenters;
    
    return result;
  }
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
  
  Function4perl(&computeVisualPolyhedra, "computeVisualPolyhedra(WeightedComplex, Rational, $)");

  Function4perl(&computeBoundedVisual, "computeBoundedVisual(WeightedComplex, $, $, Rational, Matrix<Rational>, Array<String>)");
  
  Function4perl(&boundingBox, "compute_bounding_box(Matrix<Rational>,$)");
  
}}