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

	Computes all [[LATTICE...]- related properties
	*/

#include "polymake/client.h"
#include "polymake/Set.h"
#include "polymake/Array.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Map.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/linear_algebra_tools.h"
#include "polymake/tropical/LoggingPrinter.h"


namespace polymake { namespace tropical {

	typedef Map< std::pair<int,int>, Vector<Integer> > LatticeMap;

	/*
	 * @brief Computes [[LATTICE_NORMAL_SUM]]
	 */
	void computeLatticeNormalSum(perl::Object cycle) {
		LatticeMap latticeNormals = cycle.give("LATTICE_NORMALS");
		int ambient_dim = cycle.give("FAN_AMBIENT_DIM");
		IncidenceMatrix<> codimOneCones = cycle.give("CODIMENSION_ONE_POLYTOPES");
		Vector<Integer> weights = cycle.give("WEIGHTS");
		IncidenceMatrix<> codimInc = cycle.give("MAXIMAL_AT_CODIM_ONE");

		//This will contain the result
		Matrix<Integer> summatrix(0,ambient_dim);

		//Iterate over all codim one faces
		for(int facet = 0; facet < codimOneCones.rows(); facet++) {
			//This will contain the weighted sum of the lattice normals
			Vector<Integer> result = zero_vector<Integer>(ambient_dim);
			Set<int> adjacentCones = codimInc.row(facet);
			//Go through all adjacent cones
			for(Entire<Set<int> >::iterator e=entire(adjacentCones); !e.at_end(); ++e) {
				result = result +latticeNormals[std::make_pair(facet,*e)] * weights[*e];
			}
			summatrix = summatrix / result;
		}

		cycle.take("LATTICE_NORMAL_SUM") << summatrix;

	}//END computeLatticeNormalSum

	/*
	 * @brief Computes properties [[LATTICE_NORMAL_FCT_VECTOR]], [[LATTICE_NORMAL_SUM_FCT_VECTOR]]
	 */
	void computeLatticeFunctionData(perl::Object cycle) {
		//Extract properties from the cycle
		Matrix<Rational> linealitySpace = cycle.give("LINEALITY_SPACE");
			linealitySpace = tdehomog(linealitySpace);	
			int lineality_dim = linealitySpace.rows();

		Matrix<Rational> rays = cycle.give("SEPARATED_VERTICES");
			rays = tdehomog(rays);
		LatticeMap latticeNormals = cycle.give("LATTICE_NORMALS");
		Matrix<Rational> normalsums = cycle.give("LATTICE_NORMAL_SUM");
			normalsums = tdehomog(normalsums);
		IncidenceMatrix<> codimOneCones = cycle.give("SEPARATED_CODIMENSION_ONE_POLYTOPES");
		IncidenceMatrix<> maximalCones = cycle.give("SEPARATED_MAXIMAL_POLYTOPES");
		IncidenceMatrix<> coneIncidences = cycle.give("MAXIMAL_POLYTOPES");


		//Result variables
		Map<std::pair<int,int>, Vector<Rational> > summap;
		Matrix<Rational> summatrix;
		Vector<bool> balancedFaces(codimOneCones.rows());

		//Iterate over all codim 1 faces
		for(int fct = 0; fct < codimOneCones.rows(); fct++) {
			//dbgtrace << "Facet: " << fct << endl;

			Set<int> adjacentCones = coneIncidences.row(fct);
			for(Entire<Set<int> >::iterator mc = entire(adjacentCones); !mc.at_end(); ++mc) {
				//dbgtrace << "Maxcone " << *mc << endl;
				Vector<Rational> normalvector(latticeNormals[std::make_pair(fct,*mc)]);
				//Dehomogenize this by hand
				normalvector.slice(~scalar2set(0)) -= normalvector[1] * ones_vector<Rational>(normalvector.dim()-1);
				normalvector = normalvector.slice(~scalar2set(1));
				//Compute the representation of the normal vector
				summap[std::make_pair(fct,*mc)]= functionRepresentationVector(
						maximalCones.row(*mc),
						normalvector,
						rays, linealitySpace);
			}

			//Now compute the representation of the sum of the normals
			try {
				summatrix = summatrix / functionRepresentationVector(
						codimOneCones.row(fct),
						normalsums.row(fct),
						rays, linealitySpace);
			}
			catch(std::runtime_error &e) { //This goes wrong, if X is not balanced at a given codim 1 face
				summatrix /= zero_vector<Rational>(rays.rows() + lineality_dim);
			}
		}

		//Set fan properties
		cycle.take("LATTICE_NORMAL_FCT_VECTOR") << summap;
		cycle.take("LATTICE_NORMAL_SUM_FCT_VECTOR") << summatrix; 
	}//END computeLatticeFunctionData

	Function4perl(&computeLatticeNormalSum,"computeLatticeNormalSum(Cycle)");
	Function4perl(&computeLatticeFunctionData,"computeLatticeFunctionData(Cycle)");

}}
