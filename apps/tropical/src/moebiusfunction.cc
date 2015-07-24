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

	Contains functions to compute the affine transform of a cycle 
	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {

	/*
	 * @brief Returns a Vector<int>, whose i-th entry is the value of the moebius function 
	 * of the i-th entry in LATTICE_OF_FLATS->FACES.
	 * If inverse=true, the inclusion relation is inverted.
	 */
	Vector<int> moebius(perl::Object matroid, bool inverse) {
		perl::Object lattice = matroid.give("LATTICE_OF_FLATS");
		IncidenceMatrix<> flats = lattice.give("FACES");
		
		bool is_dual = flats.row(0).size() != 0;
		int start_index = (inverse != is_dual)? flats.rows()-1 : 0;
		int end_index = (inverse != is_dual)? -1 : flats.rows(); 
		int increment = (inverse != is_dual)? -1 : 1;

		

		Vector<int> result(flats.rows());
		result[start_index] = 1;
		for(int f = start_index+increment; f != end_index; f+= increment) {
			if(f == start_index) {
				result[f] = 1; continue;
			}
			int sum = 0;
			for(int g = start_index; g != f; g+= increment) {
				Set<int> inter = flats.row(f) * flats.row(g);
				if(inter.size() == (inverse? flats.row(f).size() : flats.row(g).size())) {
					sum += result[g];
				}
			}
			result[f] = - sum;
		}
		return result;
	}


	UserFunction4perl("",&moebius, "moebius(matroid::Matroid;$=0)");

}}
