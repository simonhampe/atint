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

namespace polymake { namespace atint{ 
  
//   using namespace atintlog::donotlog;
//   using namespace atintlog::dolog;
  using namespace atintlog::dotrace;

  
  using polymake::polytope::cdd_interface::solver;
  
  /**
   * This contains the intersection of two ReachableResults:
   * - The rays of the complex
   * - The maximal three-dimensional cells
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
   * This contains the facet data of a two-dimensional cone in R^3:
   * - The codimension one faces in terms of ray indices (in a ray matrix that was given when creating this object)
   * - In the same order as the faces, a list of the corresponding inequalities in std polymake format, given as a matrix
   * - The single equation of the cone, given as a vector
   */
  struct FacetData {
    IncidenceMatrix<> facets;
    Matrix<Rational> ineqs;
    Vector<Rational> eq;
  };
  
  /**
   * Describes a line with a single vertex or a family starting in such a line:
   * - The vertex of the line
   * - A list of indices in 0,..,5 indicating which rays span a 2-dimensional cell in the family:
   * 	0 = 0,1
   * 	1 = 0,2
   * 	2 = 0,3
   * 	3 = 1,2
   * 	4 = 1,3
   * 	5 = 2,3
   */
  struct VertexLine {
    Vector<Rational> vertex;
    Set<int> cells;
  };
  
  /**
   * Describes a one-dimensional family of a line with a single vertex:
   * - A matrix (with 2 rows) describing the edge of the family 
   */
  struct VertexFamily {
    Matrix<Rational> edge;
  };
  
  /**
   * Describes a single line with a bounded edge or a family starting at such a line:
   * - The vertex at 0
   * - The vertex away from 0
   * - The index of the other leaf at 0
   * - Whether the leafs at 0 span a cell
   * - Whether the leafs away from 0 span a cell
   */
  struct EdgeLine {
    Vector<Rational> vertexAtZero;
    Vector<Rational> vertexAwayZero;
    int leafAtZero;
    bool spanAtZero;
    bool spanAwayZero;
  };
  
  /**
   * Describes a family of lines with bounded edge:
   * - A matrix (with 2 rows) describing the edge at 0
   * - A matrix (with 2 rows) describing the edge away from 0
   * - The index of the other leaf at 0
   */
  struct EdgeFamily {
    Matrix<Rational> edgeAtZero;
    Matrix<Rational> edgeAwayZero;
    int leafAtZero;
  };
  
  // ------------------------------------------------------------------------------------------------
  
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
    
  //Documentation see perl wrapper
  perl::ListReturn linesInCubic(perl::Object f) {
    //First, we compute the divisor of f
    perl::Object r3 = CallPolymakeFunction("linear_nspace",3);
    perl::Object X = CallPolymakeFunction("divisor",r3,f);
    perl::Object lindom = f.give("DOMAIN");
      Matrix<Rational> lindom_rays = lindom.give("RAYS");
      IncidenceMatrix<> lindom_cones = lindom.give("MAXIMAL_CONES");
    Matrix<Rational> funmat = f.give("FUNCTION_MATRIX");
      //Rearrange so that values of terms on points can be computed as product with matrix
      funmat = funmat.col(funmat.cols()-1) | funmat.minor(All,sequence(0,funmat.cols()-1));
    
    Matrix<Rational> degree = (-1) *  unit_matrix<Rational>(3);
    degree = ones_vector<Rational>(3) / degree;
    degree = zero_vector<Rational>(4) | degree;
    
    //Then we compute all reachable points for each direction
    Vector<ReachableResult> reachable_points;
    for(int i = 0; i < 4; i++) {
	dbglog << "Computing reachable points from direction " << i << endl;
	reachable_points |= reachablePoints(f,X,i);
    }
    
    //This maps pairs of leafs to increasing numbers and vice versa
    Map<int, Map<int,int> > edge_pair_map;
    Map<int,Set<int> > index_to_pair_map;
    int count = 0;
    for(int a = 0; a <= 3; a++) {
      edge_pair_map[a] = Map<int,int>();
      for(int b = a+1; b <= 3; b++) {
	(edge_pair_map[a])[b] = count;
	index_to_pair_map[count] = (scalar2set(a) + scalar2set(b));
	count++;
      }
    }
    
    //Result variables
    
    //These variables contain the lines in the cubic
    Vector<VertexLine> vertex_line;
    Vector<VertexFamily> vertex_family;
    Vector<EdgeLine> edge_line;
    Vector<EdgeFamily> edge_family;
    
    
    //Then we compute the six intersection sets of the reachable loci, paired in complements
    //In each step we go through all pairs of maximal cones of these two complexes. We add to each
    // cone the direction of the potential bounded edge of the line and then intersect the two
    // resulting polyhedra. Afterwards we check if the result lies in X.
    
    Matrix<Rational> dummy_lineality(0,4);
    for(int i = 1; i <= 3; i++) {
      Vector<int> remaining(sequence(1,3) - i);
	dbglog << "Intersecting reachable loci from direction 0 and " << i << " and from direction " << remaining[0] << " and " << remaining[1] << endl;
	//All points reachable from direction 0 and i
	fan_intersection_result z_inter = cdd_fan_intersection(
	  reachable_points[0].rays, dummy_lineality, reachable_points[0].cells / reachable_points[0].edges,
	  reachable_points[i].rays, dummy_lineality, reachable_points[i].cells / reachable_points[i].edges,true
	);
	//All points reachable from the other two directions
	fan_intersection_result c_inter = cdd_fan_intersection(
	  reachable_points[remaining[0]].rays, dummy_lineality, reachable_points[remaining[0]].cells / reachable_points[remaining[0]].edges,
	  reachable_points[remaining[1]].rays, dummy_lineality, reachable_points[remaining[1]].cells / reachable_points[remaining[1]].edges,true
	);
	//Clean up the intersection
	DirectionIntersection z_inter_clean = cleanUpIntersection(z_inter);
	DirectionIntersection c_inter_clean = cleanUpIntersection(c_inter);
	
	IncidenceMatrix<> all_z_cones = z_inter_clean.cells / z_inter_clean.edges / z_inter_clean.points;
	IncidenceMatrix<> all_c_cones = c_inter_clean.cells / c_inter_clean.edges / c_inter_clean.points;

// 	perl::ListReturn test;
// 	  perl::Object t("WeightedComplex");
// 	  t.take("RAYS") << z_inter_clean.rays;
// 	  t.take("MAXIMAL_CONES") << all_z_cones;
// 	  t.take("USES_HOMOGENEOUS_C") << true;
// 	  test << t; return test;
	
	dbglog << "Computing potential families in between..." << endl;
	//Go through all pairs of cells, add the appropriate direction vector and intersect the two polyhedra
	Vector<Rational> dir_z = - degree.row(0) - degree.row(i);
	Vector<Rational> dir_c = - dir_z;
	for(int zc = 0; zc < all_z_cones.rows(); zc++) {
	    //We compute the number of relevant cells in zc. If the dimension of zc is <= 1, this is just zc
	    //If the dimension is 2, we take all the codimension one faces of zc that are visible 
	    //"from the direction of dir_z", i.e. the scalar product of the outer facet normal with dir_z is
	    //strictly positive. We will then do computations separately for each such codim one face
	    Vector<Set<int> > z_cones;
	    FacetData z_facets;
	    int z_dim = (all_z_cones.row(zc).size() <= 2? all_z_cones.row(zc).size()-1 : 2);
	    if(z_dim < 2) z_cones |= all_z_cones.row(zc);
	    else {
	      z_facets = computeFacets(z_inter_clean.rays,all_z_cones.row(zc));
	      z_cones |= visibleFaces(z_facets, dir_z);
	    }
	    
	    //Now we go through all cells of c and apply the same procedure
	    for(int cc = 0; cc < all_c_cones.rows(); cc++) {
	      Vector<Set<int> > c_cones;
	      FacetData c_facets;
	      int c_dim = (all_c_cones.row(cc).size() <= 2? all_c_cones.row(cc).size()-1 : 2);
	      if(c_dim < 2) c_cones |= all_c_cones.row(cc);
	      else {
		c_facets = computeFacets(c_inter_clean.rays,all_c_cones.row(cc));
		c_cones |= visibleFaces(c_facets,dir_c);
	      }
	      
	      //Now we go trough all pairs of elements of z_cones and c_cones and find families in between them
	      for(int zlist = 0; zlist < z_cones.dim(); zlist++) {
		for(int clist = 0; clist < c_cones.dim(); clist++) {
		  //First we compute the intersection of the two polyhedra we obtain if we add to the
		  //cells in zlist/clist a ray in edge direction. The result cannot have lineality space,
		  //as neither of the two factors has one (not even implicit).
		  Matrix<Rational> center = cdd_cone_intersection(
		      z_inter_clean.rays.minor(z_cones[zlist],All) / dir_z, dummy_lineality, 
		      c_inter_clean.rays.minor(c_cones[clist],All) / dir_c, dummy_lineality, true).first ;
		  if(center.rows() == 0) continue;
		  if(center.col(0) == zero_vector<Rational>(center.rows())) continue;
		  
		  //If the intersection is 0-dimensional, we just check if the maximum is attained twice,
		  //then add it as solution
		  if(center.rows() == 1) {
		    if(maximumAttainedTwice(funmat * center.row(0))) {
		      VertexLine vl;
			vl.vertex = center.row(0);
			Set<int> spans;
			if(z_dim == 2) spans += (edge_pair_map[0])[i];
			Vector<int> rem(sequence(1,3) - i);
			if(c_dim == 2) spans += (edge_pair_map[rem[0]])[rem[1]];
			vl.cells = spans;
		      vertex_line |= vl;
		    }
		    continue;
		  }//END if 0-dimensional
		  
		  //First we compute the "border" of center, i.e. the intersection with the two cones 
		  Matrix<Rational> z_border = cdd_cone_intersection(
		      z_inter_clean.rays.minor(z_cones[zlist],All), dummy_lineality,
		      center, dummy_lineality,true).first;
		  Matrix<Rational> c_border = cdd_cone_intersection(
		      c_inter_clean.rays.minor(c_cones[clist],All), dummy_lineality,
		      center, dummy_lineality, true).first;
		  
		  //Then we refine this polyhedron along the domains of f
		  Vector<Set<int> > center_cone; center_cone |= sequence(0,center.rows());
		  DirectionIntersection center_ref = cleanUpIntersection(cdd_fan_intersection(
		      center, dummy_lineality, center_cone, lindom_rays, dummy_lineality, lindom_cones,true));
		  
		  //If the intersection is 1-dimensional , we have two possibilities:
		  // (1) The center cone is equal to the intersection of the two border cones. Then the 
		  //     whole cone is a one-dimensional family of lines without bounded edge.
		  // (2) The center cone intersects the border cones each in a vertex. Then we check if all 
		  //     cells of center_ref are contained in X. If so, we have an isolated line. 
		  if(center.rows() == 2 ) {
		    if(z_border.rows() == 2) {
		      VertexFamily vf;
			vf.edge = z_border;
		      vertex_family |= vf;
		    }//END case (1)
		    else {
		      bool found_bad_cell = false;
		      for(int cencone = 0; cencone < center_ref.edges.rows(); cencone++) {
			Matrix<Rational> cencone_rays = center_ref.rays.minor(center_ref.edges.row(cencone),All);
			Vector<Rational> interior_point = 
			    accumulate(rows(cencone_rays),
				      operations::add()) / accumulate(cencone_rays.col(0),operations::add());
			if(!maximumAttainedTwice(funmat*interior_point)) {
			  found_bad_cell = true; break;
			}
		      }//END iterate edges of refined center cell
		      //If all edges lie in X, we add the line as isolated line
		      if(!found_bad_cell) {
			EdgeLine el;
			  el.vertexAtZero = z_border.row(0);
			  el.vertexAwayZero = c_border.row(0);
			  el.leafAtZero = i;
			  el.spanAtZero = (z_dim == 2);
			  el.spanAwayZero = (c_dim == 2);			
			edge_line |= el;
		      }
		    }//END case (2)
		    continue;
		  }//END if 1-dimensional
		  
		  
		  //TODO: 2-dimensional case
		  
		}//END iterate relevant cones from c
	      }//END iterate relevant cones from z
	      
	      
	    }//END iterate j-k-reachable cones
	    
	}//END iterate 0-i-reachable cones
	
	
	
    }//END compute intersections of reachable loci
    
    dbglog << "Cleaning up result" << endl;
    
    // Step I: Clean up each type individually by checking if one contains the other
        
    Set<int> double_vertices;
    for(int vl = 0; vl < vertex_line.dim(); vl++) {
      if(!double_vertices.contains(vl)) {
	for(int ovl = 0; ovl < vertex_line.dim(); ovl++) {
	  if(vl != ovl) {
	    if(vertex_line[vl].vertex == vertex_line[ovl].vertex && 
		(vertex_line[vl].cells * vertex_line[ovl].cells).size() == vertex_line[ovl].cells.size() ) {
		double_vertices += ovl;
	    }	  
	  }
	}
      }
    }//END clean up vertex_line
    vertex_line = vertex_line.slice(~double_vertices);
    
    //TODO: For now we only check if both vertices are equal
//     Set<int> double_lines;
//     for(int el = 0; el < edge_line.dim(); el++) {
//       if(!double_lines.contains(el)) {
// 	for(int oel = el+1; oel < edge_line.dim(); oel++) {
// 	    if( edge_line[el].vertexAtZero == edge_line[oel].vertexAtZero &&
// 	        edge_line[el].vertexAwayZero == edge_line[oel].vertexAwayZero)
// 	      double_lines += oel;
// 	}
//       }
//     }
//     edge_line = edge_line.slice(~double_lines);
    
    
    //Create corresponding line objects ...............................................................
    perl::ListReturn result;
    
    //Create vertex_line objects:
    // If two rays in such an object span a 2-dim-cell, this is computed as follows:
    // We check, how far the line in direction of the sum of the two corr. rays lies in X
    // If all of it lies in X, the 2-dim cell is just the vertex + the two rays. If not, let
    // w be the end vertex of the line. Then we have two 2-dim. cells: conv(vertex,w) + each of the rays
    for(int ivert = 0; ivert < vertex_line.dim(); ivert++) {
      perl::Object var("WeightedComplex");
	Matrix<Rational> var_rays = degree / vertex_line[ivert].vertex ;
	Vector<Set<int> > var_cones;
	//Find all rays that are NOT involved in a 2-dim cell
	Set<int> rays_in_cells;
	  for(Entire<Set<int> >::iterator vcells = entire(vertex_line[ivert].cells); !vcells.at_end(); vcells++) {
	    rays_in_cells += index_to_pair_map[*vcells];
	  }
	Set<int> rays_not_in_cells = sequence(0,4) - rays_in_cells;
	//Add those rays not in 2-dim. cells
	for(Entire<Set<int> >::iterator i = entire(rays_not_in_cells); !i.at_end(); i++) {
	  var_cones|= (scalar2set(4) + scalar2set(*i));
	}
	//Now we compute the 2-dimensional cells as described above
	for(Entire<Set<int> >::iterator span = entire(vertex_line[ivert].cells); !span.at_end(); span++) {
	    Vector<Rational> md_vector = maximalDistanceVector(
		    vertex_line[ivert].vertex, 
		    accumulate(rows(degree.minor(index_to_pair_map[*span],All)),operations::add()),
		    lindom_rays, lindom_cones, funmat);
	    if(md_vector.dim() == 0) {
	      var_cones |= (scalar2set(4) + index_to_pair_map[*span]);
	    }
	    else {
	      var_rays /= md_vector;
	      Vector<int> dirs(index_to_pair_map[*span]);
	      var_cones |= (scalar2set(4) + scalar2set(var_rays.rows()-1) + dirs[0]);
	      var_cones |= (scalar2set(4) + scalar2set(var_rays.rows()-1) + dirs[1]);
	    }
	}
	var.take("RAYS") << var_rays;
	var.take("MAXIMAL_CONES") << var_cones;
	var.take("USES_HOMOGENEOUS_C") << true;
      result << var;  
    }
    
    //Create vertex_family objects: Find the direction spanned by the family and only add the remaining three
    //as rays
    for(int fvert = 0; fvert < vertex_family.dim(); fvert++) {
      Vector<Rational> dir;
	if(vertex_family[fvert].edge(0,0) == 0) dir = vertex_family[fvert].edge.row(0);
	if(vertex_family[fvert].edge(1,0) == 0) dir = vertex_family[fvert].edge.row(1);
	if(dir.dim() == 0) dir = vertex_family[fvert].edge.row(0) - vertex_family[fvert].edge.row(1);
      int missing_dir = 0;
      if(dir[1] == 0 && dir[2] == 0) missing_dir = 3;
      if(dir[1] == 0 && dir[3] == 0) missing_dir = 2;
      if(dir[2] == 0 && dir[3] == 0) missing_dir = 1;
      Matrix<Rational> var_rays = vertex_family[fvert].edge / degree.minor(~scalar2set(missing_dir),All);
      Vector<Set<int> > var_cones;
      for(int r = 2; r < 5; r++) {
	var_cones |= (sequence(0,2) + r);
      }
      perl::Object var("WeightedComplex");
	var.take("RAYS") << var_rays;
	var.take("MAXIMAL_CONES") << var_cones;
	var.take("USES_HOMOGENEOUS_C") << true;
      result << var;
      
    }
    
    
    
    //Create edge_lines
    // Two-dimensional cells at each end are computed as for vertex_line
    for(int el = 0; el < edge_line.dim(); el++) {
      perl::Object var("WeightedComplex");
	Matrix<Rational> var_rays = edge_line[el].vertexAtZero / edge_line[el].vertexAwayZero / degree;
	Vector<Set<int> > var_cones;
	var_cones |= sequence(0,2);
	if(edge_line[el].spanAtZero) {
	  Vector<Rational> z_mdvector = maximalDistanceVector(
	      var_rays.row(0), degree.row(0) + degree.row(edge_line[el].leafAtZero),
		lindom_rays, lindom_cones, funmat);
	  if(z_mdvector.dim() == 0) {
	    var_cones |= (scalar2set(0) + scalar2set(2) + scalar2set(edge_line[el].leafAtZero+2));
	  }
	  else {
	    var_rays /= z_mdvector;
	    var_cones |= (scalar2set(0) + scalar2set(var_rays.rows()-1) + scalar2set(2));
	    var_cones |= 
	      (scalar2set(0) + scalar2set(var_rays.rows()-1) + scalar2set(edge_line[el].leafAtZero+2));
	  }
	}
	else {
	  var_cones |= (scalar2set(0) + scalar2set(2));
	  var_cones |= (scalar2set(0) + scalar2set(edge_line[el].leafAtZero+2));
	}
	
	Vector<int> rem(sequence(1,3) - edge_line[el].leafAtZero);	  
	if(edge_line[el].spanAwayZero) {
	  Vector<Rational> c_mdvector = maximalDistanceVector(
	      var_rays.row(1), degree.row(rem[0]) + degree.row(rem[1]),
	      lindom_rays, lindom_cones,funmat);
	  if(c_mdvector.dim() == 0) {
	    var_cones |= (scalar2set(1) + scalar2set(rem[0]+2) + scalar2set(rem[1]+2));
	  }
	  else {
	    var_rays /= c_mdvector;
	    var_cones |= (scalar2set(1) + scalar2set(var_rays.rows()-1) + scalar2set(rem[0]+2));
	    var_cones |= (scalar2set(1) + scalar2set(var_rays.rows()-1) + scalar2set(rem[1]+2));
	  }
	}
	else {	    
	  var_cones |= (scalar2set(1) + (rem[0]+2));
	  var_cones |= (scalar2set(1) + (rem[1]+2));
	}
	
	
	
	var.take("RAYS") << var_rays;
	var.take("MAXIMAL_CONES") << var_cones;
	var.take("USES_HOMOGENEOUS_C") << true;
      result << var;
      
    }

    
    
    return result;
  }//END linesInCubic
  
  // PERL WRAPPERS ///////////////////////////////////////////////////////////////////////////////////////
  
  UserFunction4perl("# @category Enumerative geometry"
		    "# This takes a tropical polynomial (using max) of degree 3 and computes all"
		    "# lines in the corresponding cubic"
		    "# @param MinMaxFunction f A tropical polynomial of degree 3 such that the corresponding "
		    "# cubic surface does not have an explicit lineality space"
		    "# @return WeightedComplex An array, containing all isolated lines in the cubic"
		    "# and each family of lines as a two-dimensional complex"
		    ,&linesInCubic,"lines_in_cubic(MinMaxFunction)");

  Function4perl(&maximalDistanceVector,"mdv(Vector<Rational>, Vector<Rational>, Matrix<Rational>, IncidenceMatrix, Matrix<Rational>)");

  
  
}}