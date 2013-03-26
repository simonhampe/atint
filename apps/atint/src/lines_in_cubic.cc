#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/normalvector.h"
#include "polymake/atint/cdd_helper_functions.h"

namespace polymake { namespace atint{ 
  
//   using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
  using namespace atintlog::dotrace;

  
  /**
   * This contains the result of reachablePoints(...):
   * - The rays of the complex
   * - The maximal two-dimensional cells in terms of the rays
   * - The maximal one-dimensional cells in terms of the rays
   */
  struct ReachableResult {
    Matrix<Rational> rays;
    IncidenceMatrix<> cells;
    IncidenceMatrix<> edges;
  };
  
  /**
   * This contains the intersection of two ReachableResults:
   * - The rays of the complex
   * - The maximal two-dimensional cells
   * - The maximal one-dimensional cells
   * - The maximal zero-dimensional cells
   */
  struct DirectionIntersection {
    Matrix<Rational> rays;
    IncidenceMatrix<> cells;
    IncidenceMatrix<> edges;
    IncidenceMatrix<> points;    
  };
  
  /**
   @brief Takes a fan_intersection_result and cleans up the result so that no cone is contained in another and that cones are sorted according to their dimension
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
  
  /**
   @brief Computes whether in a list of values the maximum is attained at least twice.
   @param Vector<Rational> values A list of values
   @return True, if the maximum is attained at least twice, false otherwise
   */
  bool maximumAttainedTwice(Vector<Rational> values) {
    if(values.dim() <= 1) return false;
    Rational max = values[0];
    int count = 1;
    for(int j = 1; j < values.dim(); j++) {
      if(values[j] > max) {
	max = values[j]; count = 1; continue;
      }
      if(values[j] == max) count++;
    }
    return count >= 2;
  }
  
  /**
   @brief This takes a cubic surface defined by a tropical polynomial f and a direction index in 0,1,2,3 and computes the set of all points p such that the line from p in the direction of e_0,-e1,..,-e3 lies in X.
   @param MinMaxFunction f A tropical polynomial (with max) of degree 3
   @param WeightedComplex X The divisor of f (in R^3)
   @param int direction Lies in 0,1,2,3 and means we consider the direction e_0 = (1,1,1) or -e_i for i > 0
   @return ReachableResult
   */
  ReachableResult reachablePoints(perl::Object f, perl::Object X, int direction) {
    
    //Extract values    
        perl::Object lindom = f.give("DOMAIN");
    Matrix<Rational> funmat = f.give("FUNCTION_MATRIX");
      //Rearrange so that values of terms on points can be computed as product with matrix
      funmat = funmat.col(funmat.cols()-1) | funmat.minor(All,sequence(0,funmat.cols()-1));
    
    IncidenceMatrix<> cones = X.give("MAXIMAL_CONES");
    IncidenceMatrix<> codim = X.give("CODIM_1_FACES");
    IncidenceMatrix<> codim_by_max = X.give("CODIM_1_IN_MAXIMAL_CONES");
      codim_by_max = T(codim_by_max);
    Matrix<Rational> rays = X.give("RAYS");
    Set<int> vertices = X.give("VERTICES");
    
    Matrix<Rational> degree = (-1) *  unit_matrix<Rational>(3);
      degree = ones_vector<Rational>(3) / degree;
      degree = zero_vector<Rational>(4) | degree;
    if(direction < 0 || direction > 3) {
      throw std::runtime_error("Wrong direction index. Must lie in 0,1,2,3");
    }
    
    //Find the ray in X corresponding to the chosen directions
    int dir_index = -1;
    for(int r = 0; r < rays.rows(); r++) {
      if(rays.row(r) == degree.row(direction)) {
	  dir_index = r; break;
      }
    }
    if(dir_index == -1) {
      throw std::runtime_error("Cannot find direction ray in surface. Maybe not a cubic?");
    }
	

    dbgtrace << "Computing edges of cones with given direction..." << endl;
	
    //Find all cones that use the chosen direction and 
    //keep only the codimension one locus of theses
    Set<int> d_edges;
    for(int c = 0; c < cones.rows(); c++) {
	if(cones.row(c).contains(dir_index)) {
	    d_edges += codim_by_max.row(c);
	}
    }
    
    dbgtrace << "Intersecting with linearity domains" << endl;
    
    //For each of these edges, add +/-direction as lineality
    //Then refine the resulting complex along the linearity domains of f
    Set<int> used_rays_in_edges = accumulate(rows(codim.minor(d_edges,All)),operations::add());
    Matrix<Rational> d_rays = rays.minor(used_rays_in_edges,All);
    IncidenceMatrix<> d_cone_matrix = codim.minor(d_edges,used_rays_in_edges);
    Vector<Set<int> > d_cone_list;
       for(int dc = 0; dc < d_cone_matrix.rows(); dc++) {
	  d_cone_list |= d_cone_matrix.row(dc);
       }
    Matrix<Rational> d_lineality(0,degree.cols()); 
      d_lineality /= degree.row(direction);
    perl::Object d_complex("WeightedComplex");
      d_complex.take("RAYS") << d_rays;
      d_complex.take("MAXIMAL_CONES") << d_cone_list;
      d_complex.take("LINEALITY_SPACE") << d_lineality;
      d_complex.take("USES_HOMOGENEOUS_C") << true;
    d_complex = CallPolymakeFunction("intersect_container",d_complex,lindom);
    
    dbgtrace << "Projecting vertices " << endl;
    
    //We now go through all (new) vertices of this complex and "project" them onto the one-dimensional
    //complex of the edges defined by d_edges(i.e. we find an edge intersecting the ray (vertex + R*direction))
    Matrix<Rational> d_ref_rays = d_complex.give("RAYS");
    Set<int> d_ref_vertices = d_complex.give("VERTICES");
    Set<int> old_vertices = used_rays_in_edges * vertices;
    for(Entire<Set<int> >::iterator vert = entire(d_ref_vertices); !vert.at_end(); vert++) {
//       dbgtrace << "Vertex " << *vert << endl;
      //First we check if this vertex is an "old" vertex, i.e. already exists in the d_edges complex
      int old_index = -1;
      for(Entire<Set<int> >::iterator oldvert = entire(old_vertices); !oldvert.at_end(); oldvert++) {
	if(rays.row(*oldvert) == d_ref_rays.row(*vert)) {
	    old_index = *oldvert; break;
	}
      }
      if(old_index >= 0) continue;
      
      //Otherwise we go through all edges of d_edges until we find one that intersects (vertex + R*direction)
      //This is computed as follows: Assume p is a vertex of an edge of d_edges and that w is the direction 
      //from p into that edge (either a ray or p2-p1, if the edge is bounded). Then we compute a linear representation of (vertex - p_1) in terms of w and 
      //direction. If a representation exists and the second edge generator is also a vertex, we have to check 
      //that the coefficient of w = p2-p1 is in between 0 and 1, otherwise it has to be > 0
      //
      //If we find an intersecting edge, we refine it (in d_cone_list) such that it contains the intersection point.
      for(int ed = 0; ed < d_cone_list.dim(); ed++) {
// 	dbgtrace << "Edge " << ed << endl;
	Matrix<Rational> edge_generators = d_rays.minor(d_cone_list[ed],All);
	Vector<Rational> p1 = edge_generators(0,0) == 0? edge_generators.row(1) : edge_generators.row(0);
	bool bounded = edge_generators(0,0) == edge_generators(1,0);
	Vector<Rational> w;
	  if(bounded) w = edge_generators.row(1) - edge_generators.row(0);
	  else w = (edge_generators(0,0) == 0? edge_generators.row(0) : edge_generators.row(1));
	Vector<Rational> lin_rep = linearRepresentation(d_ref_rays.row(*vert) - p1, (w / degree.row(direction)));
	//Check that: 
	// - There is a representation
	// - The coefficient of w is > 0 (and < 1 if bounded)
	if(lin_rep.dim() > 0) {
	    if(lin_rep[0] > 0 && (!bounded || lin_rep[0] < 1)) {
	      //Then we add a vertex, stop searching for an edge and go to the next vertex
	      Vector<Rational> new_vertex = p1 + lin_rep[0]*w;
	      d_rays /= new_vertex;
	      Vector<int> edge_index_list(d_cone_list[ed]);
	      d_cone_list = d_cone_list.slice(~scalar2set(ed));
	      Set<int> one_cone; one_cone += edge_index_list[0]; one_cone += (d_rays.rows()-1);
	      Set<int> other_cone; other_cone += edge_index_list[1]; other_cone += (d_rays.rows()-1);
	      d_cone_list |= one_cone;
	      d_cone_list |= other_cone;
	      
	      break;
	    }
	}
	
	
      }//END intersect with edges.
    }//END project vertices
    
    dbgtrace << "Refine along projections" << endl;
    
    //Now take the refined edge complex, add lineality again and intersect it with d_complex to get
    //the final refined complex
    perl::Object d_vert_complex("WeightedComplex");
      d_vert_complex.take("RAYS") << d_rays;
      d_vert_complex.take("MAXIMAL_CONES") << d_cone_list;
      d_vert_complex.take("LINEALITY_SPACE") << d_lineality;
      d_vert_complex.take("USES_HOMOGENEOUS_C") << true;
    perl::Object final_refined = CallPolymakeFunction("intersect_container",d_complex, d_vert_complex);
      Matrix<Rational> final_rays = final_refined.give("RAYS");
      Set<int> final_vertices = final_refined.give("VERTICES");
      IncidenceMatrix<> final_cones = final_refined.give("MAXIMAL_CONES");
      IncidenceMatrix<> final_codim = final_refined.give("CODIM_1_FACES");
      IncidenceMatrix<> final_co_in_max = final_refined.give("CODIM_1_IN_MAXIMAL_CONES");
      IncidenceMatrix<> final_max_to_co = T(final_co_in_max);
      int final_direction_index = -1;
      for(int fr = 0; fr < final_rays.rows(); fr++) {					   
	if(final_rays.row(fr) == degree.row(direction)) {
	  final_direction_index = fr; break;
	}
      }
      
    //First find all codimension one cells (i.e. edges) whose span is generated by direction
    Set<int> direction_edges;
    for(int cc = 0; cc < final_codim.rows(); cc++) {
      Matrix<Rational> cc_rays = final_rays.minor(final_codim.row(cc),All);
      Vector<Rational> cc_span;
	if(cc_rays(0,0) == cc_rays(1,0)) cc_span = cc_rays.row(0) - cc_rays.row(1);
	else cc_span = (cc_rays(0,0) == 0? cc_rays.row(0) : cc_rays.row(1));
      if(direction == 0) {
	if(cc_span[1] == cc_span[2] && cc_span[2] == cc_span[3]) direction_edges += cc;
      }
      else {
	Vector<int> entries_to_check( (sequence(1,3) - direction) );
	if(cc_span[entries_to_check[0]] == 0 && cc_span[entries_to_check[1]] == 0) direction_edges += cc;
      }
    }
    Set<int> codimension_blacklist; //This will also later include edges we already used for iteration
      codimension_blacklist += direction_edges;
    
    dbgtrace << "Compute reachable 2-cells" << endl;
      
    //Find all maximal cells which have direction as ray and use them to iterate the refined complex:
    // In each step we find all maximal cells, that share an edge "in the right direction" with one 
    // of the currently selected cells and whose interior point is in the support of f. These cells
    // then serve as reference in the next iteration step. An edge in the right direction is an edge, whose 
    // span is not generated by direction and which has not been chosen before as edge
    //
    //Also, in each iteration we keep track of vertices lying in "edges in right direction", whose adjacent
    // cell is not in the support of f. Afterwards we apply a similar procedure to all such vertices
    // that are not in the reachable complex yet to find the one-dimensional part of this complex.
    Set<int> reachable_2_cells; //Indices of maximal cells that are "reachable"
    Set<int> problematic_vertices; //Vertices we might have to check afterwards to compute 1-dim. part
    Set<int> reference_cells = T(final_cones).row(final_direction_index);
    while(reference_cells.size() > 0) {
      reachable_2_cells += reference_cells;
      //First compute all edges in right direction
      Set<int> edges_in_direction = 
	accumulate(rows(final_max_to_co.minor(reference_cells,All)),operations::add()) - 
	codimension_blacklist;
      //... and add them to edges we don't need in the future
      codimension_blacklist += edges_in_direction;
      //For each of these edges, compute the adjacent cell in the correct direction and check if its 
      //in the support of X
      Set<int> next_reference_cells;
      for(Entire<Set<int> >::iterator eid = entire(edges_in_direction); !eid.at_end(); eid++) {
	int adjacent_cell =  *( (final_co_in_max.row(*eid) - reference_cells).begin());
	//Compute interior vector of adjacent cell
	Vector<Rational> ip = accumulate(rows(final_rays.minor(final_cones.row(adjacent_cell),All)),
						      operations::add()) /
					   (final_cones.row(adjacent_cell) * final_vertices).size();
	//Compute values of f on ip
	Vector<Rational> ipv = funmat * ip;
	//If the cell is in the support of X, add it to the set of reference cells for the next iteration
	//Otherwise keep the vertices as potentially problematic
	if(maximumAttainedTwice(ipv)) {
	  next_reference_cells += adjacent_cell;
	}
	else {
	  problematic_vertices += final_codim.row(*eid) * final_vertices;
	}
      }//END iterate edges in direction
      reference_cells = next_reference_cells;
    }//END iterate 2-cells
    
    //Now we have already computed the 2-dimensional part of the reachable set.
//     Set<int> rays_in_two_set = accumulate(rows(final_cones.minor(reachable_2_cells,All)),operations::add());
//     perl::Object reach_two("WeightedComplex");
//       reach_two.take("RAYS") << final_rays.minor(rays_in_two_set,All);
//       reach_two.take("MAXIMAL_CONES") << final_cones.minor(reachable_2_cells,rays_in_two_set);
//       reach_two.take("USES_HOMOGENEOUS_C") << true;
   
    dbgtrace << "Find problematic vertices " << endl;
      
    //Remove from the "problematic" vertices all those that were not problematic after all:
    //A vertex is only problematic, if it is only contained in one-valent codimension one faces of reach_two
    // (counting only those not in the span of the direction)
    Set<int> not_so_problematic_vertices;
    Set<int> reach_two_codim_indices =
      accumulate(rows(final_max_to_co.minor(reachable_2_cells,All)),operations::add());
    Set<int> reach_two_codim_nospan_indices = reach_two_codim_indices - direction_edges;
    IncidenceMatrix<> reach_two_codim = final_codim.minor(reach_two_codim_nospan_indices,All);
    IncidenceMatrix<> reach_two_co_in_max =
			final_co_in_max.minor(reach_two_codim_nospan_indices,reachable_2_cells);
    for(Entire<Set<int> >::iterator pv = entire(problematic_vertices); !pv.at_end(); pv++) {
      Set<int> containing_codim = T(reach_two_codim).row(*pv);
      for(Entire<Set<int> >::iterator cocodim = entire(containing_codim); !cocodim.at_end(); cocodim++) {
	if( reach_two_co_in_max.row(*cocodim).size() > 1) {
	    not_so_problematic_vertices += *pv; continue;
	}
      }
    }//END sort out non-problematic vertices
    problematic_vertices -= not_so_problematic_vertices;
    
    //Now we apply the same procedure as before. This time we have edges as reference cells
    //We start with the edges in direction span containing our problematic vertices that "point away" 
    //from the reach-2-locus (i.e. are not a codimension one face of it)
    Set<int> interesting_dir_edges = direction_edges - reach_two_codim_indices; //Set of all dir edges not in reach_two
    Set<int> reference_edges = accumulate(rows(T(final_codim).minor(problematic_vertices,All)),operations::add());
      reference_edges *= interesting_dir_edges;
    
    Set<int> blacklisted_edges = reach_two_codim_indices;
    Set<int> reachable_1_cells;
    
    dbgtrace << "Compute reachable 1-cells" << endl;
    
    //Go through all reference edges. Test on each one, if it lies in the support. If so, we take
    // the "next" edge as a new reference edge.
    while(reference_edges.size() > 0) {
      blacklisted_edges += reference_edges;
      Set<int> next_reference_edges;
      for(Entire<Set<int> >::iterator re = entire(reference_edges); !re.at_end(); re++) {
	//Compute an interior point and compute values on it
	Vector<Rational> ip = accumulate(rows(final_rays.minor(final_codim.row(*re),All)),
						      operations::add()) /
					   (final_codim.row(*re) * final_vertices).size();
	Vector<Rational> ipv = funmat * ip;
	if(maximumAttainedTwice(ipv)) {
	    reachable_1_cells += *re;
	    //Find the "next" edge, i.e. a direction edge, sharing a vertex with this edge
	    //that is NOT a codimension one edge of reach_two and that is NOT blacklisted
	    next_reference_edges += (  (accumulate(rows(T(final_codim).minor(final_codim.row(*re)*final_vertices,All)),operations::add())*direction_edges) - blacklisted_edges);
	}
      }//END iterate reference edges
      reference_edges = next_reference_edges;
    }//END iterate 1-cells
    
    dbgtrace << "Prepare result " << endl;
    
    //Now we also have the one-dim. part
//     Set<int> rays_in_one_set = accumulate(rows(final_codim.minor(reachable_1_cells,All)),operations::add());
//     perl::Object reach_one("WeightedComplex");
//       reach_one.take("RAYS") << final_rays.minor(rays_in_one_set,All);
//       reach_one.take("MAXIMAL_CONES") << final_codim.minor(reachable_1_cells,rays_in_one_set);
//       reach_one.take("USES_HOMOGENEOUS_C") << true;
    Set<int> all_rays_used =
      accumulate(rows(final_cones.minor(reachable_2_cells,All)),operations::add()) + 
      accumulate(rows(final_codim.minor(reachable_1_cells,All)),operations::add()); 
    Matrix<Rational> reachable_rays = final_rays.minor(all_rays_used,All);
    IncidenceMatrix<> reach_cells = final_cones.minor(reachable_2_cells, all_rays_used);
    IncidenceMatrix<> reach_edges = final_codim.minor(reachable_1_cells, all_rays_used);
      
    ReachableResult result;
      result.rays = reachable_rays;
      result.cells = reach_cells;
      result.edges = reach_edges;
    return result;
    
  }//END reachablePoints
  
  //Documentation see perl wrapper
  perl::ListReturn linesInCubic(perl::Object f) {
    //First, we compute the divisor of f
    perl::Object r3 = CallPolymakeFunction("linear_nspace",3);
    perl::Object X = CallPolymakeFunction("divisor",r3,f);
    
    //Then we compute all reachable points for each direction
    Vector<ReachableResult> reachable_points;
    for(int i = 0; i < 4; i++) {
	dbglog << "Computing reachable points from direction " << i << endl;
	reachable_points |= reachablePoints(f,X,i);
    }
    
    //Then we compute the six intersection sets of the reachable loci
    //This contains the intersection of the 0-reachable points with the i-reachable points 
    // at i-1
    Vector<DirectionIntersection> zero_reachable;
    //If zero_reachable[i] contains (0 cap i) and the other two directions are j and k, then
    // at position i we have (j cap k)
    Vector<DirectionIntersection> complement_reachable;
    Matrix<Rational> dummy_lineality(0,4);
    for(int i = 1; i <= 3; i++) {
	dbglog << "Intersecting reachable loci from direction 0 and " << i << endl;
	fan_intersection_result z_inter = cdd_fan_intersection(
	  reachable_points[0].rays, dummy_lineality, reachable_points[0].cells / reachable_points[0].edges,
	  reachable_points[i].rays, dummy_lineality, reachable_points[i].cells / reachable_points[i].edges,true
	);
	Vector<int> remaining(sequence(1,3) - i);
	fan_intersection_result c_inter = cdd_fan_intersection(
	  reachable_points[remaining[0]].rays, dummy_lineality, reachable_points[remaining[0]].cells / reachable_points[remaining[0]].edges,
	  reachable_points[remaining[1]].rays, dummy_lineality, reachable_points[remaining[1]].cells / reachable_points[remaining[1]].edges,true
	);
	//Clean up the intersection
	DirectionIntersection z_inter_clean = cleanUpIntersection(z_inter);
	DirectionIntersection c_inter_clean = cleanUpIntersection(c_inter);
	
	 
    }
    
    return perl::ListReturn();
    
  }//END linesInCubic
  
  // PERL WRAPPERS /////////////////////////////////////////
  
  UserFunction4perl("# @category Enumerative geometry"
		    "# This takes a tropical polynomial (using max) of degree 3 and computes all"
		    "# lines in the corresponding cubic"
		    "# @param MinMaxFunction f A tropical polynomial of degree 3"
		    "# @return WeightedComplex An array, containing all isolated lines in the cubic"
		    "# and each family of lines as a two-dimensional complex"
		    ,&linesInCubic,"lines_in_cubic(MinMaxFunction)");
  
//   Function4perl(&reachablePoints,"rp(MinMaxFunction,WeightedComplex,$)");
  
}}