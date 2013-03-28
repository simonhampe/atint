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
#include "polymake/atint/lines_in_cubic_reachable.h"

namespace polymake { namespace atint{ 
  
//   using namespace atintlog::donotlog;
  using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;

  
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
   @brief This takes a two-dimensional cone in R^3 in terms of its rays and vertices in homogeneous coordinates and computes all codimension one faces "visible" from a given direction. I.e. it computes the codimension one locus and keeps all faces whose normal has a positive scalar product with direction
   */
  Vector<Set<int> > visibleFaces(Matrix<Rational> rays, Vector<Rational> direction) {
    //Compute codimension one locus
        
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
    
    Matrix<Rational> degree = (-1) *  unit_matrix<Rational>(3);
    degree = ones_vector<Rational>(3) / degree;
    degree = zero_vector<Rational>(4) | degree;
    
    //Then we compute all reachable points for each direction
    Vector<ReachableResult> reachable_points;
    for(int i = 0; i < 4; i++) {
	dbglog << "Computing reachable points from direction " << i << endl;
	reachable_points |= reachablePoints(f,X,i);
    }
    
    //Result variables
    
    //All isolated fourvalent lines
    Matrix<Rational> isolated_vertices(0,4);
    
    //Then we compute the six intersection sets of the reachable loci, paired in complements
    //In each step we go through all pairs of maximal cones of these two complexes. We add to each
    // cone the direction of the potential bounded edge of the line and then intersect the two
    // resulting polyhedra. Afterwards we check if the result lies in X.
    
    //This contains the intersection of the 0-reachable points with the i-reachable points 
    // at i-1
    Vector<DirectionIntersection> zero_reachable;
    //If zero_reachable[i] contains (0 cap i) and the other two directions are j and k, then
    // at position i we have (j cap k)
    Vector<DirectionIntersection> complement_reachable;
    
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

	dbglog << "Computing potential families in between..." << endl;
	//Go through all pairs of cells, add the appropriate direction vector and intersect the two polyhedra
	Vector<Rational> dir_z = - degree.row(0) - degree.row(i);
	Vector<Rational> dir_c = - dir_z;
	for(int zc = 0; zc < all_z_cones.rows(); zc++) {
	    //We take only the "visible" cells from the edge direction. In the case of 0/1 dimenional
	    //cells, this is always the complete cell, in the case of 2-dimensional cells, this is a 
	    //part of the codimension one-locus of that cell
	    Matrix<Rational> zc_rays = z_inter_clean.rays.minor(all_z_cones.row(zc),All);
	    Vector<Set<int> > zc_cones;
	      if(zc_rays.rows() == 1) zc_cones |= scalar2set(0);
	      if(zc_rays.rows() == 2) zc_cones |= sequence(0,2);
	      if(zc_rays.rows() > 2) zc_cones |= visibleFaces(zc_rays, dir_z) ;
	  
	    for(int cc = 0; cc < all_c_cones.rows(); cc++) {
	      std::pair< Matrix<Rational>, Matrix<Rational> > connection = cdd_cone_intersection(
		z_inter_clean.rays.minor(all_z_cones.row(zc),All) / dir_z,dummy_lineality,
		c_inter_clean.rays.minor(all_c_cones.row(cc),All) / dir_c,dummy_lineality,true);
	      Matrix<Rational> conn_rays = connection.first;
	      Matrix<Rational> conn_lin = connection.second;
	      
	      //If the intersection is not empty, we extract the corresponding (potential) family and
	      //check if it is contained in X
		
	      // CASE I: The intersection is 0-dimensional: This is always an isolated line
	      if(conn_rays.rows() == 1 && conn_lin.rows() == 0) {
		insert_rays(isolated_vertices, conn_rays,true,true);		
	      }
	    }//END iterate complement cones
	}//END iterate reachable cones
	
	
	
    }//END compute intersections of reachable loci
    
    //TODO: Clean up result
    
    //Create corresponding line objects
    perl::ListReturn result;
    
    for(int l = 0; l < isolated_vertices.rows(); l++) {
      perl::Object var("WeightedComplex");
	var.take("RAYS") << isolated_vertices.minor(scalar2set(l),All) / degree;
	Vector<Set<int> > var_rays;
	for(int i = 1; i <= 4; i++) {
	  var_rays |= (scalar2set(0) + scalar2set(i));
	}
	var.take("MAXIMAL_CONES") << var_rays;
	var.take("USES_HOMOGENEOUS_C") << true;
      result << var;
    }
    
    
    return result;
  }//END linesInCubic
  
  // PERL WRAPPERS ///////////////////////////////////////////////////////////////////////////////////////
  
  UserFunction4perl("# @category Enumerative geometry"
		    "# This takes a tropical polynomial (using max) of degree 3 and computes all"
		    "# lines in the corresponding cubic"
		    "# @param MinMaxFunction f A tropical polynomial of degree 3"
		    "# @return WeightedComplex An array, containing all isolated lines in the cubic"
		    "# and each family of lines as a two-dimensional complex"
		    ,&linesInCubic,"lines_in_cubic(MinMaxFunction)");
  
//   Function4perl(&reachablePoints,"rp(MinMaxFunction,WeightedComplex,$)");
  
}}