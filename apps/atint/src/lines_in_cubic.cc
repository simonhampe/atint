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
#include "polymake/atint/lines_in_cubic.h"
#include "polymake/atint/lines_in_cubic_helper.h"
#include "polymake/atint/lines_in_cubic_reachable.h"

namespace polymake { namespace atint{ 
  
//   using namespace atintlog::donotlog;
//   using namespace atintlog::dolog;
  using namespace atintlog::dotrace;

  
  using polymake::polytope::cdd_interface::solver;
    
  // ------------------------------------------------------------------------------------------------  
    
  //Documentation see perl wrapper
  perl::Object linesInCubic(perl::Object f) {
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
			  Vector<int> rem(sequence(1,3) - i);
			  el.maxDistAtZero = 
			    maximalDistanceVector(el.vertexAtZero, 
						  degree.row(0) + degree.row(el.leafAtZero),
						  lindom_rays, lindom_cones, funmat);
			  el.maxDistAwayZero =
			    maximalDistanceVector(el.vertexAwayZero,
						  degree.row(rem[0]) + degree.row(rem[1]),
						  lindom_rays, lindom_cones, funmat);
			edge_line |= el;
		      }
		    }//END case (2)
		    continue;
		  }//END if 1-dimensional
		  
		  //If the cone is 2-dimensional, we have to refine, project and refine again,
		  // similarly to the procedure used when computing the reachable points
		  if(center.rows() > 2) {
		    
		    //TODO: 2-dimensional case
// 		    edge_family |= computedEdgeFamilies(center_ref, z_border, c_border, i, lindom_rays, lindom_cones, funmat);
		  }
		  
		  
		  
		}//END iterate relevant cones from c
	      }//END iterate relevant cones from z
	      
	      
	    }//END iterate j-k-reachable cones
	    
	}//END iterate 0-i-reachable cones
	
	
	
    }//END compute intersections of reachable loci
    
    // .............................................................................................
    dbglog << "Cleaning up result" << endl;
    
    // Step I: Clean up each type individually by checking if one contains the other or has to be glued
    // together somehow
        
    //Vertex line: Only check if one is contained in the other
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
    
    //For vertex families we have to check if two families can be glued together
    std::list<VertexFamily> vf_queue;
    for(int vf = 0; vf < vertex_family.dim(); vf++) {
      vf_queue.push_back(vertex_family[vf]);
    }
    Vector<VertexFamily> approved_vfamilies;
    while(vf_queue.size() > 0) {
      VertexFamily fam = vf_queue.front(); vf_queue.pop_front();
      bool is_approved = true;
      for(int vf = 0; vf < approved_vfamilies.dim(); vf++) {
	//Case I: Both are identical: Just ignore the new one and add nothing new
	if( (approved_vfamilies[vf].edge.row(0) == fam.edge.row(0) &&
	      approved_vfamilies[vf].edge.row(1) == fam.edge.row(1)) ||
	    (approved_vfamilies[vf].edge.row(0) == fam.edge.row(1) &&
	      approved_vfamilies[vf].edge.row(1) == fam.edge.row(0))) {
	  is_approved = false;
	  break;
	}
	
	//Case II: At least one consists of two vertices, one of which is equal to
	// a vertex of the other. Then we glue both together, add the new element to the vf_queue and
	// go to the next vf_queue element
	if( (approved_vfamilies[vf].edge(0,0) + approved_vfamilies[vf].edge(1,0) == 2) ||
	    (fam.edge(0,0) + fam.edge(1,0) == 2) ) {
	  int avf_equal = -1;
	  int fam_equal = -1;
	  if(approved_vfamilies[vf].edge.row(0) == fam.edge.row(0)) avf_equal = 0; fam_equal = 0;
	  if(approved_vfamilies[vf].edge.row(0) == fam.edge.row(1)) avf_equal = 0; fam_equal = 1;
	  if(approved_vfamilies[vf].edge.row(1) == fam.edge.row(0)) avf_equal = 1; fam_equal = 0;
	  if(approved_vfamilies[vf].edge.row(1) == fam.edge.row(1)) avf_equal = 1; fam_equal = 1;
	  
	  if(avf_equal != -1) {
	    int other_avf = avf_equal == 0? 1 : 0;
	    int other_fam = fam_equal == 0? 1 : 0;
	    VertexFamily replacement;
	      replacement.edge = approved_vfamilies[vf].edge.row(other_avf) / fam.edge.row(other_fam);
	    vf_queue.push_back(replacement);
	    approved_vfamilies = approved_vfamilies.slice(~scalar2set(vf));
	    is_approved = false;
	    break;
	  } //END create new vertex family
	}//END Case II	
      }//END iterate approved families
      if(is_approved) approved_vfamilies |= fam;
    }//END clean up vertex families
    vertex_family = approved_vfamilies;
    
    //For edge lines, we check if two are equal or can be glued together
    std::list<EdgeLine> el_queue;
    for(int el = 0; el < edge_line.dim(); el++) {
      el_queue.push_back(edge_line[el]);
    }
    Vector<EdgeLine> approved_elfamilies;
    while(el_queue.size() > 0) {
      EdgeLine fam = el_queue.front(); el_queue.pop_front();
      bool is_approved = true;
      for(int el = 0; el < approved_elfamilies.dim(); el++) {
	//Interesting things only happen if both have the same direction
	if(approved_elfamilies[el].leafAtZero == fam.leafAtZero) {
	  //Compatibility = Both families can be glued, i.e. in each direction either:
	  // - the vertices at that end agree -> keep vertex with largest span
	  // - one vertex has unbounded span and the other vertex differs by a positive multiple of 
	  //   the edge direction -> keep vertex with unbounded span
	  // - one vertex has bounded span, which contains the other vertex. -> keep vertex with bounded span
	  bool isCompatibleAtZero = false;
	  Vector<Rational> vertexAtZero;
	  Vector<Rational> mdistAtZero;
	  bool spanAtZero;
	  if(approved_elfamilies[el].vertexAtZero == fam.vertexAtZero) {
	    isCompatibleAtZero = true;
	    vertexAtZero = fam.vertexAtZero;
	    spanAtZero = approved_elfamilies[el].spanAtZero || fam.spanAtZero;
	  }
	  else {
	    //Check if difference of vertices is multiple of edge direction
	    Vector<Rational> edir = degree.row(0) + degree.row(fam.leafAtZero);
	    Rational vdiff = vertexDistance(approved_elfamilies[el].vertexAtZero,
					     fam.vertexAtZero, edir);
	    //If so, check if span is compatible
	    if(vdiff != 0) {
	      if(vdiff > 0 && approved_elfamilies[el].spanAtZero) {
		Rational sdiff = vertexDistance(approved_elfamilies[el].vertexAtZero,
					       approved_elfamilies[el].maxDistAtZero, edir);
		isCompatibleAtZero = sdiff == 0 || sdiff >= vdiff;
		vertexAtZero = approved_elfamilies[el].vertexAtZero;
		mdistAtZero = approved_elfamilies[el].maxDistAtZero;
		spanAtZero = true;
	      }		
	      if(vdiff < 0 && fam.spanAtZero) {
		Rational sdiff = vertexDistance(fam.vertexAtZero, fam.maxDistAtZero,edir);
		isCompatibleAtZero = sdiff == 0 || sdiff >= -vdiff;
		vertexAtZero = fam.vertexAtZero;
		mdistAtZero = fam.maxDistAtZero;
		spanAtZero = true;
	      }
	    }//END if is multiple of edge direction
	  }//END check compat at zero
	  if(!isCompatibleAtZero) continue;
	  
	  bool isCompatibleAwayZero = false;
	  Vector<Rational> vertexAwayZero;
	  Vector<Rational> mdistAwayZero;
	  bool spanAwayZero;
	  if(approved_elfamilies[el].vertexAwayZero == fam.vertexAwayZero) {
	    isCompatibleAwayZero = true;
	    vertexAwayZero = fam.vertexAwayZero;
	    spanAwayZero = approved_elfamilies[el].spanAwayZero || fam.spanAwayZero;
	  }
	  else {
	    //Check if difference of vertices is multiple of edge direction
	    Vector<int> rem(sequence(1,3) - fam.leafAtZero);
	    Vector<Rational> edir = degree.row(rem[0]) + degree.row(rem[1]);
	    Rational vdiff = vertexDistance(approved_elfamilies[el].vertexAwayZero,
					     fam.vertexAwayZero, edir);
	    //If so, check if span is compatible
	    if(vdiff != 0) {
	      if(vdiff > 0 && approved_elfamilies[el].spanAwayZero) {
		Rational sdiff = vertexDistance(approved_elfamilies[el].vertexAwayZero,
					       approved_elfamilies[el].maxDistAwayZero, edir);
		isCompatibleAwayZero = sdiff == 0 || sdiff >= vdiff;
		vertexAwayZero = approved_elfamilies[el].vertexAwayZero;
		mdistAwayZero = approved_elfamilies[el].maxDistAwayZero;
		spanAwayZero = true;
	      }		
	      if(vdiff < 0 && fam.spanAwayZero) {
		Rational sdiff = vertexDistance(fam.vertexAwayZero, fam.maxDistAwayZero,edir);
		isCompatibleAwayZero = sdiff == 0 || sdiff >= -vdiff;
		vertexAwayZero = fam.vertexAwayZero;
		mdistAwayZero = fam.maxDistAwayZero;
		spanAwayZero = true;
	      }
	    }//END if is multiple of edge direction
	  }//END check compat away zero
	  if(!isCompatibleAwayZero) continue;
	  
	  //If we arrive here, we have to glue
	  EdgeLine newel;
	    newel.vertexAtZero = vertexAtZero;
	    newel.vertexAwayZero = vertexAwayZero;
	    newel.leafAtZero = fam.leafAtZero;
	    newel.spanAtZero = spanAtZero;
	    newel.spanAwayZero = spanAwayZero;
	    newel.maxDistAtZero = mdistAtZero;
	    newel.maxDistAwayZero = mdistAwayZero;
	  approved_elfamilies = approved_elfamilies.slice(~scalar2set(el));
	  el_queue.push_back(newel);
	  is_approved = false;
	  break;	  
	}//END if same direction
      }//END iterate approved families
      if(is_approved) approved_elfamilies |= fam;
    }//END clean up edge lines
    edge_line = approved_elfamilies;
    
    // STEP II: A family of one type may be contained in a family of a different type
    
    //Check if any edge line is contained in a spanning vertex line
    Set<int> redundant_el;
    for(int el = 0; el < edge_line.dim(); el++) {
      for(int vl = 0; vl < vertex_line.dim(); vl++) {
	//One vertex has to be equal and the other vertex has to be in the span of the vertex line
	if(edge_line[el].vertexAtZero == vertex_line[vl].vertex ||
	   edge_line[el].vertexAwayZero == vertex_line[vl].vertex) {
	  Vector<Rational> other_vertex;
	  Set<int> dir_indices;
	  if(edge_line[el].vertexAtZero == vertex_line[vl].vertex) {
	      other_vertex = edge_line[el].vertexAwayZero;
	      dir_indices = sequence(1,3) - scalar2set(edge_line[el].leafAtZero);
	  }
	  else {
	      other_vertex = edge_line[el].vertexAtZero;
	      dir_indices = scalar2set(0) + scalar2set(edge_line[el].leafAtZero);
	  }
	  bool span_has_dir = false;
	  for(Entire<Set<int> >::iterator c = entire(vertex_line[vl].cells); !c.at_end(); c++) {
	    if(  (dir_indices * index_to_pair_map[*c]).size() == 2) {
	      span_has_dir = true; break;
	    }
	  }
	  bool is_compatible = false;
	  if(span_has_dir) {
	    Vector<Rational> edir = accumulate(rows(degree.minor(dir_indices,All)), operations::add());
	    Rational vdiff = vertexDistance(vertex_line[vl].vertex, other_vertex, edir);
	    if(vdiff > 0) {
	      Vector<Rational> mdist = maximalDistanceVector(vertex_line[vl].vertex, edir, lindom_rays, lindom_cones, funmat);
	      Rational mdiff = vertexDistance(vertex_line[vl].vertex, mdist, edir);
	      is_compatible = mdiff == 0 || mdiff >= vdiff;
	    }
	  }
	  if(is_compatible) redundant_el += el;
	}
      }
    }
    edge_line = edge_line.slice(~redundant_el);
    
    //Check if any vertex_line is contained in a vertex_family
    Set<int> redundant_vl;
    for(int vl = 0; vl < vertex_line.dim(); vl++) {
      for(int vf = 0; vf < vertex_family.dim(); vf++) {
	//Check if the vertex is a vertex of the family and the spans contain the direction
	int fam_dir = vertexFamilyDirection(vertex_family[vf]);
	if( vertex_line[vl].vertex == vertex_family[vf].edge.row(0) ||
	     vertex_line[vl].vertex == vertex_family[vf].edge.row(1) ){
	  Set<int> cells = vertex_line[vl].cells;
	  bool is_compatible = true;
	  for(Entire<Set<int> >::iterator c = entire(cells); !c.at_end(); c++) {
	    if(!index_to_pair_map[*c].contains(fam_dir)) {
	      is_compatible = false; break;
	    }
	  }
	  if(is_compatible) redundant_vl += vl;
	}
      }//END iterate vertex_family
    }//END compare vertex_line to vertex_family
    vertex_line = vertex_line.slice(~redundant_vl);
    
   
    
    //Create corresponding line objects ...............................................................    
    dbglog << "Creating lines as complexes..." << endl;
    
    perl::Object result("LinesInCubic");

    
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
      if(vertex_line[ivert].cells.size() == 0) 
	result.add("LIST_ISOLATED_NO_EDGE",var);
      else 
	result.add("LIST_FAMILY_FIXED_VERTEX", var);
    }
    
    //Create vertex_family objects: Find the direction spanned by the family and only add the remaining three
    //as rays
    for(int fvert = 0; fvert < vertex_family.dim(); fvert++) {
      int missing_dir = vertexFamilyDirection(vertex_family[fvert]);
      Matrix<Rational> var_rays = vertex_family[fvert].edge / degree.minor(~scalar2set(missing_dir),All);
      Vector<Set<int> > var_cones;
      for(int r = 2; r < 5; r++) {
	var_cones |= (sequence(0,2) + r);
      }
      perl::Object var("WeightedComplex");
	var.take("RAYS") << var_rays;
	var.take("MAXIMAL_CONES") << var_cones;
	var.take("USES_HOMOGENEOUS_C") << true;
      result.add("LIST_FAMILY_MOVING_VERTEX", var);
      
    }
    
    
    
    //Create edge_lines
    // Two-dimensional cells at each end are computed as for vertex_line
    for(int el = 0; el < edge_line.dim(); el++) {
      perl::Object var("WeightedComplex");
	Matrix<Rational> var_rays = edge_line[el].vertexAtZero / edge_line[el].vertexAwayZero / degree;
	Vector<Set<int> > var_cones;
	var_cones |= sequence(0,2);
	if(edge_line[el].spanAtZero) {
	  Vector<Rational> z_mdvector = edge_line[el].maxDistAtZero;
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
	  Vector<Rational> c_mdvector = edge_line[el].maxDistAwayZero;
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
      if(edge_line[el].spanAtZero || edge_line[el].spanAwayZero) 
	result.add("LIST_FAMILY_FIXED_EDGE",var);
      else result.add("LIST_ISOLATED_EDGE",var);
      
    }

    dbglog << "Done." << endl;
    
    return result;
  }//END linesInCubic
  
  // PERL WRAPPERS ///////////////////////////////////////////////////////////////////////////////////////
  
  UserFunction4perl("# @category Enumerative geometry"
		    "# This takes a tropical polynomial (using max) of degree 3 and computes all"
		    "# lines in the corresponding cubic"
		    "# @param MinMaxFunction f A tropical polynomial of degree 3 such that the corresponding "
		    "# cubic surface does not have an explicit lineality space"
		    "# @return LinesInCubic An object containing the complete result of the computation"
		    ,&linesInCubic,"lines_in_cubic(MinMaxFunction)");

  
  
  
}}