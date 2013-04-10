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
 * Contains some helper functions for the lines-in-cubic method
 */

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/normalvector.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/atint/cdd_helper_functions.h"
#include "polymake/atint/lines_in_cubic_reachable.h"
#include "polymake/atint/lines_in_cubic.h"
#include "polymake/atint/WeightedComplexRules.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  //using namespace atintlog::dotrace;
  
  using polymake::polytope::cdd_interface::solver;
  
  /**
   @brief Takes a fan_intersection_result and cleans up the result so that no cone is contained in another and that cones are sorted according to their dimension. 
   @param fan_intersection_result fir
   @return DirectionIntersection
   */
  DirectionIntersection cleanUpIntersection(fan_intersection_result fir) {
    DirectionIntersection result;
      result.rays = fir.rays;
    IncidenceMatrix<> fir_cones = fir.cones;
    //First we sort all cells according to their dimension
    Set<int> cell_set;
    Set<int> edge_set;
    Set<int> point_set;
    for(int fc = 0; fc < fir_cones.rows(); fc++) {
	 if(fir_cones.row(fc).size() > 2) cell_set += fc;
	 if(fir_cones.row(fc).size() == 2) edge_set += fc;
	 if(fir_cones.row(fc).size() == 1) point_set += fc;
    }
    //Go through all edges, compare to cells, remove redundant ones
    Set<int> redundant_edges;
    for(Entire<Set<int> >::iterator e = entire(edge_set); !e.at_end(); e++) {
	bool found_container = false;
	for(Entire<Set<int> >::iterator c = entire(cell_set); !c.at_end(); c++) {
	    if ( (fir_cones.row(*c) * fir_cones.row(*e)).size() == fir_cones.row(*e).size()) {
		found_container = true; break;
	    }
	}
	if(found_container) redundant_edges += (*e);
    }
    edge_set -= redundant_edges;
    //Same for points
    Set<int> redundant_points;
    for(Entire<Set<int> >::iterator p = entire(point_set); !p.at_end(); p++) {
	bool found_container = false;
	Set<int> containers = cell_set + edge_set;
	for(Entire<Set<int> >::iterator c = entire(containers); !c.at_end(); c++) {
	    if ( (fir_cones.row(*c) * fir_cones.row(*p)).size() == fir_cones.row(*p).size()) {
		found_container = true; break;
	    }
	}
	if(found_container) redundant_points += (*p);
    }
    point_set -= redundant_points;
    
    result.cells = fir_cones.minor(cell_set,All);
    result.edges = fir_cones.minor(edge_set,All);
    result.points = fir_cones.minor(point_set,All);
    
    return result;
  }//END cleanUpIntersection
    
  // ------------------------------------------------------------------------------------------------
  
  /**
   @brief This takes a (two-dimensional) cone in R^3 in terms of a subset of rays and computes all codimension one faces. 
   @return A FacetData object (see above)
   */
  FacetData computeFacets(const Matrix<Rational> &rays,const Set<int> &cone) {
    //Compute facet equations and store them
    solver<Rational> sv;
    std::pair<Matrix<Rational>, Matrix<Rational> > ceq = sv.enumerate_facets(
	zero_vector<Rational>() | rays.minor(cone,All),Matrix<Rational>(0,rays.cols()+1),true,false);
    FacetData result;
      result.eq = ceq.second.row(0).slice(~scalar2set(0));
      result.ineqs = ceq.first.minor(All,~scalar2set(0));
    //Now go through all inequalities and find the corresponding vertices
    Vector<Set<int> > facets;
    for(int i = 0; i < result.ineqs.rows(); i++)  {
      Set<int> iset;
      for(Entire<Set<int> >::const_iterator c = entire(cone); !c.at_end(); c++) {
	if(result.ineqs.row(i) * rays.row(*c) == 0) iset += (*c);
      }
      facets |= iset;
    }
    //Each facet has to have at least one vertex
    for(int f = 0; f < facets.dim(); f++) {
      if(rays.col(0).slice(facets[f]) == zero_vector<Rational>(facets[f].size())) {
	facets = facets.slice(~scalar2set(f));
	result.ineqs = result.ineqs.minor(~scalar2set(f),All);
      }
    }
    
    result.facets = IncidenceMatrix<>(facets);
    
//     //Normlize inequalities such that normals lie in plane spanned by cone
//     
//     //Take the norm squared of the linear part of the equation
//     Rational eq_norm =
//       accumulate(attach_operation(result.eq.slice(~scalar2set(0)),operations::square()),operations::add());
//     //For each facet normal f, subtract the equation times g the scalar product of the linear parts f'*g'
//     for(int f = 0; f < result.ineqs.rows(); f++) {
//       Rational scal_prod = result.ineqs.row(f).slice(~scalar2set(0)) * 
// 			    result.eq.slice(~scalar2set(0));
//       result.ineqs.row(f) = result.ineqs.row(f) - ( scal_prod / eq_norm) * result.eq;
//     }
      
    return result;
  } 
    
  void testcf(Matrix<Rational> m, Set<int> c) {
    FacetData fd = computeFacets(m,c);
    pm::cout << fd.ineqs << endl;
    pm::cout << fd.eq << endl;
  }
    
  // ------------------------------------------------------------------------------------------------  
    
  /**
   @brief This takes a result of computeFacets and finds all the facets visible from "direction", i.e. all facets whose outer normal has strict positive scalar product with direction (it doesn't matter which normal we take: The direction must lie in the span of the two-dimensional cone, so the equation g of the cone is zero on direction. Any two representatives of a facet normal only differ by a multiple of g).
   */
  Vector<Set<int> > visibleFaces(FacetData fd, Vector<Rational> direction) {
    Vector<Set<int> > result;
    for(int f = 0; f < fd.ineqs.rows(); f++) {
      if(fd.ineqs.row(f) * direction < 0) { //<0, since ineqs has the inner normals
	result |= fd.facets.row(f);
      }
    }      
    return result;
  }
  
  void test(Matrix<Rational> m, Set<int> c, Vector<Rational> v) {
    FacetData fd = computeFacets(m,c);
      pm::cout << fd.facets << endl;
      pm::cout << fd.ineqs << endl;
      pm::cout << fd.eq << endl;
    pm::cout << visibleFaces(fd, v) << endl;
  }
  
  // ------------------------------------------------------------------------------------------------  
  
  /**
   @brief This takes a vertex in the cubic (whose function's domain is describes by frays and fcones) and a direction and computes the vertex w farthest away from vertex in this direction, such that the convex hull of vertex and w still lies in X. It returns the empty vertex, if the complete half-line lies in X.
   */
  Vector<Rational> maximalDistanceVector(const Vector<Rational> &vertex, const Vector<Rational> &direction, 
					  const Matrix<Rational> &frays, const IncidenceMatrix<> &fcones, const Matrix<Rational> &funmat) {
    //Create the one-dimensional half-line from vertex
    Matrix<Rational> hl_rays = vertex / direction;
    Vector<Set<int> > hl_cones; hl_cones |= sequence(0,2);
    Matrix<Rational> lin(0, hl_rays.cols());
    //Intersect with f-domain
    DirectionIntersection ref_line = cleanUpIntersection(
	cdd_fan_intersection(hl_rays, lin, hl_cones, frays,lin,fcones,true));
    //Find vertex
    int v_index = -1;
    for(int r = 0; r < ref_line.rays.rows(); r++) {
      if(ref_line.rays.row(r) == vertex) {v_index = r; break;}
    }
    IncidenceMatrix<> rays_in_cones = T(ref_line.edges);
    //Go through edges, starting at vertex and check if it is contained in X
    int current_edge = *(rays_in_cones.row(v_index).begin());
    int current_vertex = v_index;
    //When the current edge has only one vertex, we're done
    do {
      //Check if the edge lies in X
      Vector<Rational> interior_point = accumulate
	(rows(ref_line.rays.minor(ref_line.edges.row(current_edge),All)), operations::add()) /
	accumulate(ref_line.rays.minor(ref_line.edges.row(current_edge),All).col(0),operations::add());
      if(!maximumAttainedTwice(funmat * interior_point)) return ref_line.rays.row(current_vertex);
      else {
	current_vertex = *( (ref_line.edges.row(current_edge) - current_vertex).begin());
	if(ref_line.rays.row(current_vertex)[0] == 0) current_vertex = -1;
	else {
	    current_edge = *( (rays_in_cones.row(current_vertex) - current_edge).begin());
	}
      }
    }while(current_vertex >= 0);
  
    return Vector<Rational>();
  }
  
  // ------------------------------------------------------------------------------------------------ 
  
  /**
   @brief Takes two vectors v1 and v2 and a standard direction e_i + e_j for i,j in {0,..,3} and computes the rational number r, such that v1 + r*direction = v2. Returns also zero, if no such number exists. It also accepts a 0-dimensional vertex for v2 and will return 0 in that case.
   */
  Rational vertexDistance(const Vector<Rational> &v1, const Vector<Rational> &v2, const Vector<Rational> &direction) {
      if(v2.dim() == 0) return 0;
      Vector<Rational> diff = v2 - v1;
      Rational div = 0;
      for(int i = 1; i <= 3; i++) {
	if( (diff[i] == 0 && direction[i] != 0) || (diff[i] != 0 && direction[i] == 0)) return 0;
	if(diff[i] != 0) {
	    Rational d = diff[i] / direction[i];
	    if(div == 0) div = d;
	    else {
	      if(d != div) return 0;
	    }
	}
      }
      return div;
  }
  
  // ------------------------------------------------------------------------------------------------ 
  
  /**
   @brief Takes a vertex family and computes the index of the standard direction in 0,..,3 corresponding to its edge
   */
  int vertexFamilyDirection(VertexFamily f) {
    Vector<Rational> dir;
	if(f.edge(0,0) == 0) dir = f.edge.row(0);
	if(f.edge(1,0) == 0) dir = f.edge.row(1);
	if(dir.dim() == 0) dir = f.edge.row(0) - f.edge.row(1);
      if(dir[1] == 0 && dir[2] == 0) return 3;
      if(dir[1] == 0 && dir[3] == 0) return 2;
      if(dir[2] == 0 && dir[3] == 0) return 1;
      return 0;
  }
  
  // ------------------------------------------------------------------------------------------------ 
  
  /**
   @brief Computes all edge families lying in a 2-dimensional cone for a given direction
   @param DirectionIntersection cone A cone, refined along f
   @param Matrix<Rational> z_border The intersection of the cone with a cone in the 0-i-rechable locus
   @param Matrix<Rational> c_border The intersection of the cone with a cone in the j-k-reachable locus
   @param int leafAtZero The index of the leaf together with 0
   @param Matrix<Rational> funmat The function matrix of f, made compatible for vector multiplication
   @return std::pair< Vector<EdgeFamily>, Vector<EdgeLine> > A list of all edge families and edge lines lying in the cone	    
   */
  std::pair<Vector<EdgeFamily>, Vector<EdgeLine> > computeEdgeFamilies(DirectionIntersection cone, 
					 const Matrix<Rational> &z_border, 
					 const Matrix<Rational> &c_border,int leafAtZero, const Matrix<Rational> &funmat) {
					  
    Matrix<Rational> degree = (-1) *  unit_matrix<Rational>(3);
    degree = ones_vector<Rational>(3) / degree;
    degree = zero_vector<Rational>(4) | degree;
					 
    //First we project all vertices of the cone onto z_border
    Matrix<Rational> z_edge_rays = z_border;
    Vector<Set<int> > z_edges; z_edges |= sequence(0,2);
    Vector<Rational> direction = degree.row(0) + degree.row(leafAtZero);
    for(int dr = 0; dr < cone.rays.rows(); dr++) {
    if(cone.rays(dr,0) == 1) {
      //We go through all edges of z_edges until we find one that intersects (vertex + R_>=0 *direction)
      //This is computed as follows: Assume p is a vertex of an edge of z_edges and that w is the direction 
      //from p into that edge (either a ray or p2-p1, if the edge is bounded). Then we compute a linear representation of (vertex - p_1) in terms of w and 
      //direction. If a representation exists and the second edge generator is also a vertex, we have to check 
      //that the coefficient of w = p2-p1 is in between 0 and 1, otherwise it has to be > 0
      //
      //If we find an intersecting edge, we refine it (in z_edges) such that it contains the intersection point.
      for(int ze = 0; ze < z_edges.dim(); ze++) {
	Matrix<Rational> edge_generators = z_border.minor(z_edges[ze],All);
	Vector<Rational> p1 = edge_generators(0,0) == 0? edge_generators.row(1) : edge_generators.row(0);
	bool bounded = edge_generators(0,0) == edge_generators(1,0);
	Vector<Rational> w;
	  if(bounded) w = edge_generators.row(1) - edge_generators.row(0);
	  else w = (edge_generators(0,0) == 0? edge_generators.row(0) : edge_generators.row(1));
	Vector<Rational> lin_rep = linearRepresentation(cone.rays.row(dr) - p1, (w / direction));
	// 	Check that: 
	// 	- There is a representation
	// 	- The coefficient of w is > 0 (and < 1 if bounded)
	if(lin_rep.dim() > 0) {
	  if(lin_rep[0] > 0 && (!bounded || lin_rep[0] < 1)) {
	    //Then we add a vertex, stop searching for an edge and go to the next vertex
	    Vector<Rational> new_vertex = p1 + lin_rep[0]*w;
	    z_edge_rays /= new_vertex;
	    Vector<int> edge_index_list(z_edges[ze]);
	    z_edges = z_edges.slice(~scalar2set(ze));
	    Set<int> one_cone; one_cone += edge_index_list[0]; one_cone += (z_edge_rays.rows()-1);
	    Set<int> other_cone; other_cone += edge_index_list[1]; other_cone += (z_edge_rays.rows()-1);
	    z_edges |= one_cone;
	    z_edges |= other_cone;
	    break;
	  }
	}//END if linear rep exists
      }//END iterate z_edges
    }//END if vertex
    }//END project vertices
     
    //Then refine the cone along the new z_edges - direction
    //and compute codim one data
    Matrix<Rational> dummy_lineality(0, z_border.cols());
    DirectionIntersection refined_cone = cleanUpIntersection(
		  cdd_fan_intersection(cone.rays, dummy_lineality, cone.cells,
					z_edge_rays, dummy_lineality / direction, z_edges,true));
    CodimensionOneResult codimData =
	      calculateCodimOneData(refined_cone.rays, refined_cone.cells, true, dummy_lineality, IncidenceMatrix<>());
	      
    
    //Find the directionr ray
    int dir_index = -1;
    for(int r = 0; r < refined_cone.rays.rows(); r++) {
      if(refined_cone.rays.row(r) == direction) {
	dir_index = r; break;
      }
    }
    
    //Find all edges containing this ray
								   
     
  }//END computeEdgeFamilies
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
//   Function4perl(&maximalDistanceVector,"mdv(Vector<Rational>, Vector<Rational>, Matrix<Rational>, IncidenceMatrix, Matrix<Rational>)");
//   Function4perl(&vertexDistance,"vd(Vector<Rational>,Vector<Rational>, Vector<Rational>)");
  
}}