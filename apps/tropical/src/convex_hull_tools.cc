
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

	Implements convex_hull_tools.h
	*/

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/tropical/LoggingPrinter.h"
#include "polymake/polytope/cdd_interface.h"
#include "polymake/tropical/convex_hull_tools.h"
#include "polymake/tropical/solver_def.h"

namespace polymake { namespace tropical {

	using namespace atintlog::donotlog;
	//using namespace atintlog::dolog;
	//using namespace atintlog::dotrace;


	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see header
	Vector<int> insert_rays(Matrix<Rational> &rays, Matrix<Rational> nrays, bool is_normalized) {
		//Normalize new rays, if necessary
		if(!is_normalized) {
			cdd_normalize_rays(nrays,true);
		}

		//Insert rays
		Vector<int> new_ray_indices;
		for(int nr = 0; nr < nrays.rows(); nr++) {
			int new_rayindex = -1;
			for(int oray = 0; oray < rays.rows(); oray++) {
				if(rays.row(oray) == nrays.row(nr)) {
					new_rayindex = oray; break;
				}
			}
			if(new_rayindex == -1) {
				rays /= nrays.row(nr);
				new_rayindex = rays.rows()-1;
			}
			new_ray_indices |= new_rayindex;
		}

		return new_ray_indices;
	}//END insert_rays

	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see header
	void cdd_normalize_rays(Matrix<Rational> &rays, bool has_leading_coordinate) {
		for(int r = 0; r < rays.rows(); r++) {
			//Vertices are normalized to have first coordinate 1
			if(rays(r,0) != 0 && has_leading_coordinate) {
				rays.row(r) *= (1/rays(r,0));
			}
			//Rays are normalized to have first nonzero coordinate +-1
			else {
				for(int c = 0; c < rays.cols(); c++) {
					if(rays(r,c) != 0) {
						rays.row(r) *= abs(1/rays(r,c));  
						break;
					}
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see header
	std::pair<Matrix<Rational>, Matrix<Rational> > cdd_cone_intersection(
			const Matrix<Rational> &xrays, const Matrix<Rational> &xlin, 
			const Matrix<Rational> &yrays, const Matrix<Rational> &ylin) {

		solver<Rational> sv;


		//Compute facets
		std::pair<Matrix<Rational>, Matrix<Rational> > x_eq =
			sv.enumerate_facets(xrays, xlin,false,false);
		std::pair<Matrix<Rational>, Matrix<Rational> > y_eq =
			sv.enumerate_facets(yrays, ylin,false,false);

		//Compute intersection rays
		std::pair<Matrix<Rational>, Matrix<Rational> > inter;
		try {
			inter = sv.enumerate_vertices(x_eq.first / y_eq.first, x_eq.second / y_eq.second, false,true);
		}
		catch(...) {
			inter.first = Matrix<Rational>(0,std::max(xrays.cols(), xlin.cols()));
			inter.second = Matrix<Rational>(0,std::max(xrays.cols(), xlin.cols()));
		}

		//normalize
		cdd_normalize_rays(inter.first, true);

		return inter;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	//Documentation see header
	fan_intersection_result cdd_fan_intersection(	
			Matrix<Rational> xrays, Matrix<Rational> xlin, IncidenceMatrix<> xcones,
			Matrix<Rational> yrays, Matrix<Rational> ylin, IncidenceMatrix<> ycones) {

		solver<Rational> sv;

		//Precompute h-representations of the x-cones and y-cones
		Vector<std::pair<Matrix<Rational>, Matrix<Rational> > > xequations;
		for(int xc = 0; xc < xcones.rows(); xc++) {
			xequations |= sv.enumerate_facets(xrays.minor(xcones.row(xc),All),
					xlin,false,false);
		}
		Vector<std::pair<Matrix<Rational>, Matrix<Rational> > > yequations;
		for(int yc = 0; yc < ycones.rows(); yc++) {
			yequations |= sv.enumerate_facets(yrays.minor(ycones.row(yc),All),
					ylin,false,false);
		}

		//Now compute intersections
		Matrix<Rational> interrays;
		Matrix<Rational> interlineality;
		bool lineality_computed = false;
		Vector<Set<int> > intercones;

		Vector<Set<int> > xcontainers;
		Vector<Set<int> > ycontainers;

		for(int xc = 0; xc < xcones.rows(); xc++) {
			for(int yc = 0; yc < ycones.rows(); yc++) {	
				//Compute intersection
				std::pair<Matrix<Rational>, Matrix<Rational> > inter;
				try {
					inter = sv.enumerate_vertices(xequations[xc].first / yequations[yc].first,
							xequations[xc].second / yequations[yc].second,
							false,true);
				}
				catch(...) {
					inter.first = Matrix<Rational>(0,std::max(xrays.cols(), xlin.cols()));
					inter.second = Matrix<Rational>(0,std::max(xrays.cols(), xlin.cols()));
				}

				if(!lineality_computed) {
					interlineality = inter.second.rows() > 0 ? 
						inter.second :
						Matrix<Rational>();
					lineality_computed = true;
				}

				//The empty cone will not be included 
				if(inter.first.rows() == 0) {
					continue;
				}

				//If cone contains no vertices (i.e. the intersection is actually 
				//empty), we leave it out
				if(inter.first.col(0) == zero_vector<Rational>(inter.first.rows())) continue;

				cdd_normalize_rays(inter.first, true);

				//Insert rays into ray list and create cone
				Set<int> new_cone_set;
				// 	bool new_ray_added = false;
				for(int r = 0; r < inter.first.rows(); r++) {
					int ray_index = -1;
					for(int oray = 0; oray < interrays.rows(); oray++) {
						if(inter.first.row(r) == interrays.row(oray)) {
							ray_index = oray; 
							new_cone_set += ray_index;
							break;
						}
					}
					if(ray_index == -1) {
						interrays /= inter.first.row(r);
						new_cone_set += (interrays.rows()-1);
						// 	      new_ray_added = true;
					}
				}

				//Make sure we don't add a cone twice
				//Also: Remember intersections that are contained in this one or contain this one
				Set<int> containedCones;
				Set<int> containerCones;
				int new_cone_index = -1;
				for(int j = 0; j < intercones.dim(); j++) {
					if(intercones[j] == new_cone_set) {
						new_cone_index = j; 
					}
					else {
						int sz = (intercones[j] * new_cone_set).size();
						if(sz == intercones[j].size()) containedCones += j;
						if(sz == new_cone_set.size()) containerCones += j;
					}
				}
				if(new_cone_index == -1) {
					intercones |= new_cone_set;
					new_cone_index = intercones.dim()-1;
					xcontainers |= Set<int>();
					ycontainers |= Set<int>();
				}

				//First add all containers from the containing cones
				for(Entire<Set<int> >::iterator sup = entire(containerCones); !sup.at_end(); sup++) {
					xcontainers[new_cone_index] += xcontainers[*sup];
					ycontainers[new_cone_index] += ycontainers[*sup];
				}
				//Add xc and yc as containers
				xcontainers[new_cone_index] += xc;
				ycontainers[new_cone_index] += yc;
				//Add all current containers to the contained cones
				for(Entire<Set<int> >::iterator sub = entire(containedCones); !sub.at_end(); sub++) {
					xcontainers[*sub] += xcontainers[new_cone_index];
					ycontainers[*sub] += ycontainers[new_cone_index];
				}



			}//END iterate ycones
		}//END iterate xcones

		//Create result:
		fan_intersection_result f;
		f.rays = interrays;
		if(interlineality.rows() == 0) interlineality = Matrix<Rational>(0,interrays.cols());
		f.lineality_space = interlineality;
		f.cones = IncidenceMatrix<>(intercones);
		f.xcontainers = IncidenceMatrix<>(xcontainers);
		f.ycontainers = IncidenceMatrix<>(ycontainers);
		return f;

	}

	// ------------------------- PERL WRAPPERS ---------------------------------------------------

	Function4perl(&cdd_cone_intersection, "cdd_cone_intersection(Matrix<Rational>,Matrix<Rational>,Matrix<Rational>,Matrix<Rational>,$)");

	Function4perl(&cdd_normalize_rays, "cdd_normalize_rays(Matrix<Rational>,$)");


}}
