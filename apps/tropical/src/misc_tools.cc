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

	Implementations of miscelleaneous tools
	*/

#include "polymake/client.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/tropical/misc_tools.h"


namespace polymake { namespace tropical { 

	std::pair<Set<int>, Set<int> > far_and_nonfar_vertices(const Matrix<Rational> &m) {
		Set<int> nonfar;
		Set<int> far;
		for(int r = 0; r < m.rows(); r++) {
			(m(r,0) == 0 ? far : nonfar) += r;
		}
		return std::pair<Set<int>,Set<int> >(far,nonfar);
	}

}}
