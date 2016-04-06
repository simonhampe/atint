/*
	This program is free software; you can redistribute it and/or
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
	Copyright (C) 2011 - 2015, Simon Hampe <simon.hampe@googlemail.com>

	Functions to compute the weight space/cone/polytope of a cycle and check
	for irreducibility.
	*/

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Array.h"
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"
#include "polymake/tropical/lattice.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/LoggingPrinter.h"

namespace polymake { namespace tropical {
	using namespace atintlog::donotlog;
	//using namespace atintlog::dolog;
	//   using namespace atintlog::dotrace;


	///////////////////////////////////////////////////////////////////////////////////////

	Matrix<Rational> cycle_weight_space(perl::Object C) {

		//Extract values
		Matrix<Rational> rays = C.give("VERTICES");
                if (rays.rows()==0)
                  return Matrix<Rational>();
                rays = tdehomog(rays);
		Matrix<Rational> linspace = C.give("LINEALITY_SPACE");
                linspace = tdehomog(linspace);

		int ambient_dim = C.give("PROJECTIVE_AMBIENT_DIM");
		int dim = C.give("PROJECTIVE_DIM");
		IncidenceMatrix<> maximal_cones = C.give("MAXIMAL_POLYTOPES");

		//If there is only one maximal polytope, this is (possibly locally) only a linear space
		if(maximal_cones.rows() == 1) {
			Matrix<Rational> linw(1,1); linw(0,0) = 1;
			return linw;
		}

		IncidenceMatrix<> codim_1_faces = C.give("CODIMENSION_ONE_POLYTOPES");
		IncidenceMatrix<> codim_in_maximal = C.give("MAXIMAL_AT_CODIM_ONE");
		IncidenceMatrix<> maximal_to_codim = T(codim_in_maximal);
		Map<std::pair<int,int>, Vector<Integer> > lattice_normals = C.give("LATTICE_NORMALS");
		//Dehomogenize lattice normals as well
		for (auto& ln : lattice_normals) {
                  lattice_normals[ln.first] = tdehomog_vec(ln.second);
		}

		if(dim == 0) {
			return unit_matrix<Rational>(maximal_cones.rows());
		}

		//dbgtrace << "Computing equivalence classes" << endl;

		//Compute equivalence classes of maximal cones
		//and while doing so identify those classes which can not be balanced (i.e. which have to have
		// weight 0).
		Vector<Set<int> > subdivision_classes;
		Set<int> zero_classes; //List of classes that must have weight 0.
		Set<int> hasBeenAdded; //List of cones that know their class
		Map<int,int> class_index; //Maps each cone to the index of its class in subdivision_classes
		for(int mc = 0; mc < maximal_cones.rows(); mc++) {
			if(!hasBeenAdded.contains(mc)) {
				Set<int> mc_class; mc_class += mc;
				class_index[mc] = subdivision_classes.dim();
				bool is_zero_class = false;
				//Elements in this queue are already in the class but their neighbours might not:
				std::list<int> queue;
				queue.push_back(mc);
				while(queue.size() > 0) {
					int node = queue.front(); queue.pop_front();
					Set<int> node_codim = maximal_to_codim.row(node);
					for(Entire<Set<int> >::iterator nc = entire(node_codim); !nc.at_end(); nc++) {
						if(codim_in_maximal.row(*nc).size() == 2) {
							int other_cone = *((codim_in_maximal.row(*nc) - node).begin());
							if(!hasBeenAdded.contains(other_cone)) {
								hasBeenAdded += other_cone;
								mc_class += other_cone;
								class_index[other_cone] = subdivision_classes.dim();
								//If it is not a zero class yet, check if it is now
								if(!is_zero_class) {
									//The lattice normals must add up to 0 mod V_{nc}
									Vector<Rational> lsum (lattice_normals[std::make_pair(*nc,node)] + lattice_normals[std::make_pair(*nc,other_cone)]);
									Matrix<Rational> vtau = rays.minor(codim_1_faces.row(*nc),All) / linspace;
									if(rank(vtau/lsum) > rank(vtau)) is_zero_class = true;
								}//END check if zero class
							}
						}//END if two-valent
					}//END iterate all codim-1-faces
				}//END iterate queue
				if(is_zero_class) zero_classes += subdivision_classes.dim();
				subdivision_classes |= mc_class;	
			}//END if not in a class yet
		}//END iterate maximal cones to find equiv. classes

		//dbgtrace << "Zero classes: " << zero_classes << endl;
		//dbgtrace << "Clas indices: " << class_index << endl;

		//Linear system matrix
		Matrix<Rational> system_matrix(0,subdivision_classes.dim());

		//dbgtrace << "Computing system matrices" << endl;

		//Iterate codimension one faces to compute local signature systems
		for(int tau = 0; tau < codim_1_faces.rows(); tau++) {

			//dbgtrace << "Computing signature neighbours" << endl;
			//Create local system matrix for tau
			Matrix<Rational> Mtau(ambient_dim+1,0);
			//Find signature neighbours
			Set<int> nonsig_neighbours;
			Set<int> all_neighbours = codim_in_maximal.row(tau);
			for(Entire<Set<int> >::iterator mc = entire(all_neighbours); !mc.at_end(); mc++) {
				//See if there is another one with the same class
				for(Entire<Set<int> >::iterator other_mc = entire(all_neighbours); !other_mc.at_end(); other_mc++) {
					if(*other_mc != *mc && class_index[*mc] == class_index[*other_mc]) {
						nonsig_neighbours += *mc; break;
					}
				}//END iterate rest of neighbours
			}//END iterate all neighbours to find signature ones
			Vector<int> sig_neighbours(all_neighbours - nonsig_neighbours);

			//dbgtrace << "Signature neighbours " << sig_neighbours << endl;

			for(int smc = 0; smc < sig_neighbours.dim(); smc++) {
				Mtau |= lattice_normals[std::make_pair(tau, sig_neighbours[smc])];
			}//END iterate signature neighbours

			//dbgtrace << "Computing lattice basis" << endl;

			if(dim > 1) {
				Matrix<Integer> lbasis = 	
					lattice_basis_of_cone(rays.minor(codim_1_faces.row(tau),All),linspace,dim-1,true);
				Mtau |= T( zero_vector<Integer>() | lbasis);
			}//END compute lattice basis


			//dbgtrace << "Appending result " << endl;

			//Compute kernel
			Matrix<Rational> Ktau = null_space(Mtau);
			//dbgtrace << "Kernel: " << Ktau << endl;
			//Compute the conversion to the total weight space:
			Set<int> remaining_classes = sequence(0,subdivision_classes.dim());
			Matrix<Rational> Ptau(Ktau.rows(), subdivision_classes.dim());
			for(int sig = 0; sig < sig_neighbours.dim(); sig++) {
				Ptau.col(class_index[sig_neighbours[sig]]) = Ktau.col(sig);
				remaining_classes -= (class_index[sig_neighbours[sig]]);
			}//END copy results

			//The remaining classes can have arbitrary weights here, so we append the appropriate unit vectors
			for(Entire<Set<int> >::iterator rc = entire(remaining_classes); !rc.at_end(); rc++) {
				Ptau /= unit_vector<Rational>(subdivision_classes.dim(),*rc);
			}//END add unit vectors for remaining classes

			//dbgtrace << "Transform: " << Ptau << endl;

			//Compute equations, attach to total system and reduce
			system_matrix /= null_space(Ptau);
			system_matrix = system_matrix.minor(basis_rows(system_matrix),All);


		}//END iterate codimension one faces

		//dbgtrace << "Transforming..." << endl;

		//To compute the final subdivision weight space, we add the equation that all zero classes
		// must have weight 0
		for(Entire<Set<int> >::iterator zc = entire(zero_classes); !zc.at_end(); zc++) {
			system_matrix /= unit_vector<Rational>(system_matrix.cols(), *zc);
		}//END add zero class equations

		//Compute the solution
		Matrix<Rational> subdiv_space = null_space(system_matrix);

		//Transform to weight space of actual complex
		Matrix<Rational> result(subdiv_space.rows(), maximal_cones.rows());
		for(int mc = 0; mc < maximal_cones.rows(); mc++) {
			result.col(mc) = subdiv_space.col(class_index[mc]);
		}//END transform result

		return result;

	}//END weight_properties

	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see perl wrapper
	bool is_irreducible(perl::Object C) {
		//First we compute the gcd of the weights of C
		Vector<Integer> weights = C.give("WEIGHTS");
		if(weights.dim() == 0) return true;
		Integer g = weights[0];
		for(int w = 1; w < weights.dim(); w++) {
			g = gcd(g,weights[w]);
		}
		if(g != 1) {
			return false;
		}

		//     Map<int,int> dummy;
		Matrix<Integer> F = C.give("WEIGHT_SPACE");
		return F.rows() == 1;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see perl wrapper
	perl::Object weight_cone(perl::Object fan, Set<int> negative_directions) {
		//Extract weight system
		Matrix<Rational> wsystem = fan.give("WEIGHT_SYSTEM");
		int N = fan.give("N_MAXIMAL_POLYTOPES");

		//Take facets of orthant and invert chosen rows
		Matrix<Rational> orthant = unit_matrix<Rational>(N);
		for(Entire<Set<int> >::iterator coord = entire(negative_directions); !coord.at_end(); coord++) {
			orthant.row(*coord) *= -1;
		}

		perl::Object cone("polytope::Cone");
		if(wsystem.rows() > 0)
			cone.take("EQUATIONS") << wsystem;
		cone.take("INEQUALITIES") << orthant;
		return cone;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Computes the decomposition polytope of a tropical variety: 
	 * It computes the irreducible varieties in the 
	 * weight cone (where the corresponding orthant is taken such that the weight vector of X 
	 * lies in that orthant). It then
	 * computes the polytope of all positive linear combinations of those irreducible varieties that 
	 * produce the original weight vector.
	 * @param Cycle fan
	 * @return polytope::Polytope
	 */
	perl::Object decomposition_polytope(perl::Object fan) {
		//Extract values
		Vector<Integer> weights = fan.give("WEIGHTS");
		Set<int> negative_entries;
		for(int i = 0; i < weights.dim(); i++) { if(weights[i] < 0) negative_entries += i;}

		perl::Object cone = weight_cone(fan,negative_entries);

		//Compute equations: The linear combination of the rays of the cone must be the weight vector of the variety
		Vector<Rational> rweights(weights);
		Matrix<Rational> crays = cone.give("RAYS");
		crays = Matrix<Rational>(makePrimitiveInteger(crays));
		crays = (- rweights) | T(crays);

		//Facets = positive orthant
		Matrix<Rational> facets = unit_matrix<Rational>(crays.cols()-1);
		facets = zero_vector<Rational>() | facets;

		perl::Object p("polytope::Polytope");
		p.take("INEQUALITIES") << facets;
		p.take("EQUATIONS") << crays;
		return p;


	}

	// ------------------------- PERL WRAPPERS ---------------------------------------------------

	Function4perl(&is_irreducible,"is_irreducible(Cycle)");
	//     Function4perl(&cycle_weight_space,"cycle_weight_space(Cycle)");
	//     Function4perl(&cycle_weight_system,"cycle_weight_system(Cycle)");
	//     Function4perl(&cycle_weight_cone,"cycle_weight_cone(Cycle)");
	Function4perl(&cycle_weight_space,"cycle_weight_space(Cycle)");
	
	UserFunction4perl("# @category Weight space"
			"# Computes the possible positive decompositions into irreducible subvarieties of the same "
			"# weight positivity signature (i.e. the weight on a cone has to have the same sign as in the cycle)"
			"# To be precise, it computes the irreducible varieties as rays of the weight cone"
			"# (where the corresponding orthant is taken such that the weight vector of X "
			"# lies in that orthant). It then computes the polytope of all positive linear combinations "
			"# of those irreducible varieties that produce the original weight vector."
			"# @param Cycle A weighted complex"
			"# @return polytope::Polytope",			
			&decomposition_polytope,"decomposition_polytope(Cycle)");

	UserFunction4perl("# @category Weight space"
			"# Takes a polyhedral complex and computes a weight cone, i.e. "
			"# intersects the [[WEIGHT_SPACE]] with a chosen orthant (by default the positive orthant)"
			"# @param Cycle X A polyhedral complex"
			"# @param Set<int> negative A subset of the coordinates {0,..,N-1}, where N is "
			"# the number of maximal cells of X. Determines the orthant to"
			"# intersect the weight space with: All coordinates in the set are negative, the others positive"
			"# If the set is not given, it is empty by default (i.e. we take the positive orthant)",
			&weight_cone,
			"weight_cone(Cycle, Set<Int>)");


}}

