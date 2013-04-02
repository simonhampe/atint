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
  using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;

  
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
    
    //Result variables
    
    //All isolated fourvalent lines, represented by the position of the apex
    Matrix<Rational> isolated_vertices(0,4);
    //All isolated lines with nontrivial bounded edge, represented by the positions of the vertices
    Vector< std::pair<Vector<Rational>, Vector<Rational> > > isolated_lines;
    //All families. Each family is represented by two cones, representing the possible placements of the
    //two vertices. Each cone is represented as a matrix of rays
    Vector< std::pair<Matrix<Rational>, Matrix<Rational> > > families;
    
    
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
		  //then add it as isolated solution
		  if(center.rows() == 1) {
		    if(maximumAttainedTwice(funmat * center.row(0))) {
		      isolated_vertices /= center.row(0); 
		    }
		    continue;
		  }//END if 0-dimensional
		  
		  //Then we refine this polyhedron along the domains of f
		  Vector<Set<int> > center_cone; center_cone |= sequence(0,center.rows());
		  DirectionIntersection center_ref = cleanUpIntersection(cdd_fan_intersection(
		      center, dummy_lineality, center_cone, lindom_rays, dummy_lineality, lindom_cones,true));
		  
		  //If the intersection is 1-dimensional, we have to check, if ALL cells of center_ref lie in X.
		  //If so, we add the result as an isolated line.
		  if(center.rows() == 2) {
		    
		  }
		  
		  
		}//END iterate relevant cones from c
	      }//END iterate relevant cones from z
	      
	      
	    }//END iterate j-k-reachable cones
	    
	}//END iterate 0-i-reachable cones
	
	
	
    }//END compute intersections of reachable loci
    
    dbglog << "Cleaning up result" << endl;
    
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
		    "# @param MinMaxFunction f A tropical polynomial of degree 3 such that the corresponding "
		    "# cubic surface does not have an explicit lineality space"
		    "# @return WeightedComplex An array, containing all isolated lines in the cubic"
		    "# and each family of lines as a two-dimensional complex"
		    ,&linesInCubic,"lines_in_cubic(MinMaxFunction)");

  Function4perl(&test,"test(Matrix<Rational>, Set<Int>, Vector<Rational>)");

}}