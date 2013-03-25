#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/atint/divisor.h"
#include "polymake/atint/normalvector.h"

namespace polymake { namespace atint{ 
  
  using namespace atintlog::donotlog;
  //using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;

  /**
   @brief This takes a cubic surface defined by a tropical polynomial f and a direction index in 0,1,2,3 and computes the set of all points p such that the line from p in the direction of e_0,-e1,..,-e3 lies in X.
   @param MinMaxFunction f
   @param int direction Lies in 0,1,2,3 and means we consider the direction e_0 = (1,1,1) or -e_i for i > 0
   @return WeightedComplex A polyhedral complex. This may be non-pure, so visualization is not immediately possible.
   */
  perl::Object reachablePoints(perl::Object f, int direction) {
   
    perl::Object r3 = CallPolymakeFunction("linear_nspace",3);
    perl::Object X = CallPolymakeFunction("divisor",r3,f);
    perl::Object lindom = f.give("DOMAIN");
    
    IncidenceMatrix<> cones = X.give("MAXIMAL_CONES");
    IncidenceMatrix<> codim = X.give("CODIM_1_FACES");
    IncidenceMatrix<> codim_by_max = X.give("CODIM_1_IN_MAXIMAL_CONES");
      codim_by_max = T(codim_by_max);
    Matrix<Rational> rays = X.give("RAYS");
    
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
	
    
    //Find all cones that use the chosen direction and 
    //keep only the codimension one locus of theses
    Set<int> d_cones;
    for(int c = 0; c < cones.rows(); c++) {
	if(cones.row(c).contains(dir_index)) {
	    d_cones += codim_by_max.row(c);
	}
    }
    
    Matrix<Rational> reach_rays = rays;
    Vector<Set<int> > reach_cones;
    
    //For each of these edges, add the direction as lineality space and
    //intersect the resulting complex with the linearity domain of f
    for(Entire<Set<int> >::iterator edge = entire(d_cones); !edge.at_end(); edge++) {
	perl::Object stripe("WeightedComplex");
	  Matrix<Rational> stripe_rays = rays.minor(codim.row(*edge),All);
	    stripe_rays /= (-degree.row(direction));
	    stripe_rays /= -degree.row(direction);
	  Vector<Set<int> > stripe_cones;
	    stripe_cones |= (sequence(0,stripe_rays.rows()-1));
	    stripe_cones |= (sequence(0,stripe_rays.rows()-2) + (stripe_rays.rows()-1));
	  stripe.take("RAYS") << stripe_rays;
	  stripe.take("MAXIMAL_CONES") << stripe_cones;
	  stripe.take("USES_HOMOGENEOUS_C") << true;
	stripe = CallPolymakeFunction("intersect_container",stripe,lindom);
	
	//Extract the vertices of the refined stripe and project them onto the edge
	Matrix<Rational> sref_rays = stripe.give("RAYS");
	Set<int> sref_verts = stripe.give("VERTICES");
	
	//We project each vertex along the direction-axis onto the space spanned by the
	//generators of the edge. We then sort them according to their second coordinate 
	//(we make sure the first edge generator is a vertex)
	Matrix<Rational> generators = stripe_rays.minor(sequence(0,stripe_rays.rows()-2),All);
	if(generators(0,0) == 0) {
	    for(int c = 0; c < generators.cols(); c++) generators(0,c).swap(generators(1,c));
	}
	
	Matrix<Rational> ordered_projected_vertices(0, sref_rays.cols());
	Vector<Rational> second_projection_coordinate;
	for(Entire<Set<int> >::iterator v = entire(sref_verts); !v.at_end(); v++) {
	    Vector<Rational> proj_vert = 
	      (generators.row(0) * sref_rays.row(*v)) * generators.row(0) +
	      (generators.row(1) * sref_rays.row(*v)) * generators.row(1);
	    int vert_spc = generators.row(1) * sref_rays.row(*v);
	    int pos = -1;
	    for(int spc = 0; spc < second_projection_coordinate.dim(); spc++) {
		if(second_projection_coordinate[spc] > vert_spc) {
		    pos = spc; break;
		}
	    }
	    
	}
	
	
	
    }
    
    return perl::Object("WeightedComplex");
    
    
    
  }
  
  Function4perl(&reachablePoints,"rp(MinMaxFunction,$)");
  
}}