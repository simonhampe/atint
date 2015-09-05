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

	*/

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/Polynomial.h"
#include "polymake/TropicalNumber.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"
#include "polymake/tropical/thomog.h"
#include "polymake/tropical/specialcycles.h"

namespace polymake { namespace tropical {

	template <typename Addition>
		perl::Object refined_hypersurface(perl::Object skeleton, Polynomial<TropicalNumber<Addition> > p) {
			Matrix<Rational> skel_rays = skeleton.give("VERTICES");
			IncidenceMatrix<> skel_cones = skeleton.give("MAXIMAL_POLYTOPES");
			Set<int> hypercones;

			Matrix<int> monomials = p.monomials_as_matrix();
			Vector<Rational> coefficients(p.coefficients_as_vector());

			for(int c = 0; c < skel_cones.rows(); c++) {
				Vector<Rational> interior = accumulate( rows(skel_rays.minor(skel_cones.row(c),~scalar2set(0))),
															operations::add());
				Vector<Rational> values = monomials * interior + coefficients;
				Rational maximum = Addition::orientation() == 1 ?
					accumulate(values, operations::min()) :
					accumulate(values, operations::max());

				int count = 0;
				for(int entry = 0; entry < values.dim(); entry++) {
					if(values[entry] == maximum) count ++;
					if(count >= 2) { hypercones += c; break;}
				}
			}

			perl::Object result(perl::ObjectType::construct<Addition>("Cycle"));
				result.take("VERTICES") << skel_rays;
				result.take("MAXIMAL_POLYTOPES") << skel_cones.minor(hypercones, All);
				result.take("WEIGHTS") << ones_vector<Integer>(hypercones.size());
			return result;
		}


	IncidenceMatrix<> translated_cones(const Matrix<Rational> &rays, const Matrix<Rational> &ref_rays,
													const IncidenceMatrix<> &cones) {
		Map<int,int> raymap;
		for(int r = 0; r < rays.rows(); r++) {
			for(int s = 0; s < ref_rays.rows(); s++) {
				if(ref_rays.row(s) == rays.row(r)) {
					raymap[r] = s; break; 
				}
			}
		}

		Vector<Set<int> > result;
		for(int c = 0; c < cones.rows(); c++) {
			result |= attach_operation(cones.row(c), pm::operations::associative_access<Map<int,int>, int>(&raymap));
		}

		return IncidenceMatrix<>(result);
	}//END translated_cones

	int index_of(const Set<int> &set, const IncidenceMatrix<> &cones) {
		for(int c = 0; c < cones.rows(); c++) {
			if(cones.row(c) == set) {
				return c;
			}
		}
		throw std::runtime_error("Cone not found");
	}


	template <typename Addition>
		Matrix<Rational> fan_system(perl::Object container, Array<perl::Object> fans) {
			Matrix<Rational> ref_rays = container.give("VERTICES");
			IncidenceMatrix<> ref_cones = container.give("MAXIMAL_POLYTOPES");

			Matrix<Rational> result(fans.size(), ref_cones.rows());

			for(int f = 0; f < fans.size(); f++) {
				cout << "Checking fan " << (f+1) << " of " << fans.size() << endl;
				Matrix<Rational> frays = fans[f].give("VERTICES");
				IncidenceMatrix<> fcones = fans[f].give("MAXIMAL_POLYTOPES");
				Vector<Integer> fweight = fans[f].give("WEIGHTS");
				IncidenceMatrix<> tcones = translated_cones(frays, ref_rays, fcones);
				for(int tc = 0; tc < tcones.rows(); tc++) {
					int t_index = index_of(tcones.row(tc), ref_cones);
					result(f,t_index) = fweight[tc];
				}
			}

			return result;

		}//END add_fans

	UserFunctionTemplate4perl("","refined_hypersurface<Addition>(Cycle<Addition>, Polynomial<TropicalNumber<Addition> >)");
	UserFunctionTemplate4perl("","fan_system<Addition>(Cycle<Addition>, Cycle<Addition>+)");

}}