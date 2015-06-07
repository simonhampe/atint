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

	This file contains the algorithms to compute bergman fans as described
	by Felipe Rincón in his paper "Computing tropical linear spaces".
	See also http://math.berkeley.edu/~felipe/bergman
	*/


#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Vector.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Map.h"
#include "polymake/tropical/specialcycles.h"
#include "polymake/tropical/LoggingPrinter.h"


namespace polymake { namespace tropical { 

	using namespace atintlog::donotlog;
	//using namespace atintlog::dolog;
	//using namespace atintlog::dotrace;

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Takes a matrix and computes all sets of column indices, s.t. the corresponding column vectors form a basis of the column space of the matrix.
	  @param Matrix<Rational> m 
	  @return IncidenceMatrix An incidence matrix whose rows correspond to the bases and where each row contains the colum indices contained in that basis.
	  */
	IncidenceMatrix<> computeMatrixBases(Matrix<Rational> m) {
		//First we determine the rank of the matrix
		int r = rank(m);

		Vector<Set<int> > result;
		Array<Set<int> > rsets = all_subsets_of_k(sequence(0,m.cols()),r);
		//pm::Subsets_of_k<Set<int> > ( sequence(0,m.cols()),r );
		for(int i = 0; i < rsets.size(); i++) {
			if(rank(m.minor(All,rsets[i])) == r) {
				result |= rsets[i];
			}
		}
		return IncidenceMatrix<>(result);
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Takes a matrix and finds all columns that are coloops in the matrix matroid, i.e. all column vectors that are contained in every basis of the column space. In other words, the rank of the matrix decreases, if I remove such a column
	  @param Matrix<Rational> m
	  @return Set<int> A set of column indices of columns that are coloops
	  */
	Set<int> computeMatrixColoops(Matrix<Rational> m) {
		//First we determine the rank of the matrix
		int r = rank(m);

		Set<int> coloops;
		for(int c = 0; c < m.cols(); c++) {
			if(rank(m.minor(All,~scalar2set(c))) < r) {
				coloops += c;
			}
		}
		return coloops;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief For a given base B of a matroid, this method computes all the sets C(k,B)-{k}, where C(k,B) is the fundamental circuit of k over B
	  @param IncidenceMatrix<> bases The set of all bases of the matroid
	  @param int mybase The row index of the relevant base in the matrix bases
	  @param Vector<int> complement The complement of the basis as an ordered list
	  @return Vector<Set<int> > The set C(k,B) - {k}, in the order of complement
	  */
	Vector<Set<int> > computeFk(const IncidenceMatrix<> &bases, int mybase, Vector<int> complement) {
		Vector<Set<int> > result;
		for(int k = 0; k < complement.dim(); k++) {
			Set<int> Fk;
			//We go through all elements i of mybase and check if B - i + k is independent, i.e. equal to a basis
			Set<int> B = bases.row(mybase);
			for(Entire<Set<int> >::iterator i = entire(B); !i.at_end(); i++) {
				Set<int> C = B; C -= (*i); C += complement[k];
				bool independent = false;
				//Check if the set C is contained in any basis
				for(int row = 0; row < bases.rows(); row++) {
					if(row != mybase) {
						if(C == bases.row(row)) {
							independent = true;
							break;
						}
					}
				}
				if(independent) {
					Fk += (*i);
				}
			}
			result |= Fk;
		}
		return result;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief For a given base B of a linear matroid, this method computes all the sets C(k,B)-{k}, where C(k,B) is the fundamental circuit of k over B. This is much faster then computeFk for general matroids
	  @param IncidenceMatrix<> bases The set of all bases of the matroid
	  @param int mybase The row index of the relevant base in the matrix bases
	  @param Matrix<Rational> m The matrix representing the matroid
	  @param Vector<int> complement The complement of the basis as an ordered list
	  @return Vector<Set<int> > The set C(k,B) - {k}, in the order of complement
	  */
	Vector<Set<int> > computeFkLinear(const IncidenceMatrix<> &bases, int mybase, const Matrix<Rational> &m, const Vector<int> &complement) {
		Vector<Set<int> > result;
		Vector<int> B(bases.row(mybase));
		//Compute the row reduced matrix s.t. the minor corresponding to mybase is the identity
		Matrix<Rational> I = inv(m.minor(All,bases.row(mybase)));
		Matrix<Rational> A = I * m;
		for(int k = 0; k < complement.dim(); k++) {
			Set<int> Fk;
			for(int i = 0; i < A.rows(); i++) {
				if(A(i,complement[k]) != 0) {
					Fk += B[i];
				}
			}
			result |= Fk;
		}
		return result;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Given a regressive compatible pair (p,L) on a basis B, computes the corresponding cone
	  @param int n The number of elements of the matroid ground set
	  @param Matrix<Rational> rays A reference to the ray matrix. New rays will be appended
	  @param Set<int> B The basis on which the pair is defined
	  @param Vector<Set<int> > Fksets The sets F_k for all k not in B
	  @param Vector<int> p The map p: [n] - B -> B, p[i] = image of complement[i]
	  @param Vector<int> L The ordering on Im(p) (first = smallest)
	  @param Vector<Set<int> > Qb The i-th entry (where i in B) contains the preimage of i under p
	  @return The ray indices of the cone in the new ray matrix
	  */
	Set<int> computeCone(int n, Matrix<Rational> &rays, const Set<int> &B, const Vector<Set<int> > &Fksets, const Vector<int> &p, const Vector<int> &L, const Vector<Set<int> > &Qb) {
		//First we compute the ray w_b for each b in B
		Matrix<Rational> wbs(0,n);
		//Compute the complement of L (i.e. the singleton blocks of the defining tree)
		Set<int> Lcomp = B - Set<int>(L);
		//dbgtrace << "Computing leaves" << endl;
		//For each element in L find the attached leafs
		Vector<Set<int> > leafsAttached(L.dim());
		for(int l = 0; l < leafsAttached.dim(); l++) { leafsAttached[l] = Set<int>();}
		for(Entire<Set<int> >::const_iterator c = entire(Lcomp); !c.at_end(); c++) {
			//Start from the right end of L to find the correct node
			bool found = false;
			for(int l = L.dim() -1; l >= 0 && !found; l--) {
				//Find a complement element k that maps to L[l] and s.t. F_k contains c
				for(int k = 0; k < Fksets.dim(); k++) {
					if(p[k] == L[l]) {
						if(Fksets[k].contains(*c)) {
							found = true;
							leafsAttached[l] += *c;
							break;
						}
					}
				}
			}
			//Also, we can already compute the ray w_c
			wbs /= unit_vector<Rational>(n,*c);
		} //End attach leaves

		//dbgtrace << "Computing block vectors" << endl;

		//Now compute the rays w_b for b in L (except for the smallest one, which gives (1,..1))
		for(int b = 1; b < L.dim(); b++) {
			Vector<Rational> wb(n);
			for(int c = b; c < L.dim(); c++) {
				wb[L[c]] = 1;
				//Add all elements from Q_c and the leaves
				Set<int> X = Qb[L[c]] + leafsAttached[c];
				for(Entire<Set<int> >::iterator x = entire(X); !x.at_end(); x++) {
					wb[*x] = 1;
				}
			}
			wbs /= wb;
		}

		//dbgtrace << "Cone rays: " << wbs << endl;
		//dbgtrace << "Checking rays" << endl;

		//Now we check for existence of the rays and compute indices
		Set<int> result;
		for(int r = 0; r < wbs.rows(); r++) {
			int index = -1;
			for(int oray = 0; oray < rays.rows(); oray++) {
				if(rays.row(oray) == wbs.row(r)) {
					index = oray; break;
				}
			}
			if(index >= 0) {
				result += index;
			}
			else {
				rays /= wbs.row(r);
				int newindex = rays.rows()-1;
				result += newindex;
			}
		}
		//dbgtrace << "Indices: " << result << endl;
		return result;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Computes the bergman fan of a matroid that is supposed to be loop- and coloop-free
	  @param int n The number of elements of the ground set of the matroid
	  @param IncidenceMatrix<> bases The bases of the matroid
	  @param bool is_linear Whether the matroid is represented by a Q-matrix
	  @param Matrix<Rational> m if is_linear is true, this contains the matrix representing the matroid
	  @return A pair consisting of the rays (in Min and not yet homogeneous) and the maximal cones.
	  */
	std::pair<Matrix<Rational>, IncidenceMatrix<> > bergman_fan(int n, const IncidenceMatrix<> &bases, bool is_linear, const Matrix<Rational> &m) {
		//Prepare result variables
		Matrix<Rational> rays(0,n);
		Vector<Set<int> > cones;
		Vector<Integer> weights;

		//Now compute cones for each basis
		Set<int> complete = sequence(0,n);
		for(int B = 0; B < bases.rows(); B++) {
			//Initialize fundamental circuits
			Vector<int> complement(complete - bases.row(B));
			Vector<Set<int> > Fksets = 
				is_linear? 	computeFkLinear(bases,B,m,complement) : 
				computeFk(bases,B,complement);

			Vector<Vector<int> >  Fklists(Fksets.dim());
			for(int k = 0; k < Fksets.dim(); k++) {
				Fklists[k] = Vector<int>(Fksets[k]);
			}

			//i-th element will contain all the k not in B mapped to i by p
			Vector<Set<int> > Qb(n); 
			for(int l = 0; l < Qb.dim(); l++) { Qb[l] = Set<int>();}

			//Now we go through all the regressive compatible pairs (p,L)

			// Contains the images p(k), more precisely: p[i] is the image of complement[i]
			Vector<int> p(complement.dim()); 
			Vector<int> L; //the ordering on Im(p) (first element = smallest, etc)
			int k = 0; //The k we currently define an image for (actually: index in complement)
			//For given k, iterations[0]..iterations[k] describe the iterations for all k' <= k
			// Semantics: (-1,-1): We just started with this k, check F_k cap L
			// (0,-1): Start adding elements from F_k \ L
			// (b,x): In the last iteration we took the element at index b in F_k and inserted it
			// before position x. Find the next valid position
			Vector<std::pair<int,int> > iterations(complement.dim());
			for(int t = 0; t < iterations.dim(); t++) {
				iterations[t] = std::make_pair(-1,-1);
			}	
			while(k != -1) {
				//If we're through all k's, make a cone -------------------------
				if(k >= complement.dim()) {
					//dbgtrace << "Compatible pair: \n: Order: " << L <<"\np: " << p << endl;
					//dbgtrace << "Qb: " << Qb << endl;
					cones |= computeCone(n,rays,bases.row(B),Fksets,p,L,Qb);
					k--;
					continue;
				}
				//Remove k from the preimage of its image
				Qb[p[k]] -= complement[k];

				//Check F_k \cap L for elements --------------------------------
				int smallestInL = -1; //contains the smallest element in Fk cap L (or -1 if none)
				int smallestx = -1; // contains the position of the above element in L
				if(iterations[k].first == -1) {
					//dbgtrace << "Iteration for " << complement[k] << " at " << iterations[k] << endl;
					iterations[k] = std::make_pair(0,-1); 
					for(int l = 0; l < L.dim(); l++) {
						if(Fksets[k].contains(L[l])) {
							smallestInL = L[l]; break;
							smallestx = l;
						}
					}
					if(smallestInL >= 0) {
						p[k] = smallestInL;
						Qb[smallestInL] += complement[k];
						k++;
						continue;
					}
				}
				//Go through the basis elements -------------------------------
				//Remove the last b we added
				if(iterations[k].second != -1) {
					L = L.slice(~scalar2set(iterations[k].second));  
				}
				//dbgtrace << "Iteration for " << complement[k] << " at " << iterations[k] << endl;
				//dbgtrace << "L is " << L << endl;
				//dbgtrace << "Smallest x is " << smallestx << endl; 
				//Have to recompute minimal element from Fk cap L (any other b from Fk must
				//be placed before that element)
				smallestx = -1;	      
				for(int l = 0; l < L.dim(); l++) {
					if(Fksets[k].contains(L[l])) {
						smallestx = l; break;
					}
				}
				//Find the next valid (b,x)
				int b = iterations[k].first;
				int x = iterations[k].second;
				bool found = false;
				if(smallestx < 0) {
					smallestx = L.dim();
				}
				do {
					x++;
					//If we arrive after the last valid position, take the next b
					if(x > smallestx) {
						x = 0; b++;
					}
					//Check that b is not in L
					bool binL = false;
					for(int l = 0; l < L.dim(); l++) {
						if(Fklists[k][b] == L[l]) {
							binL = true; break;
						}
					}
					if(binL) {
						x = 0; b++;
					}
					//Check that b is smaller than k (and the list size)
					if(b >= Fklists[k].dim()) {
						break;
					}
					if(Fklists[k][b] >= complement[k]) {
						break;
					}
					//Check validity of b
					bool isvalid = true;
					for(int l = 0; isvalid && l < complement.dim() && l < k; l++) {
						if(Fksets[l].contains(Fklists[k][b])) {
							//check if (p(complement[l]) > b), i.e. p(c[l]) occurs at position x or later
							for(int y = x; y < L.dim(); y++) {
								if(L[y] == p[l]) {
									isvalid = false; break;
								}
							}
						}
					}
					if(isvalid) {
						//If we arrive here, we found a valid position
						found = true; break;
					}
				} while(true);
				//If we found a position, continue the iteration
				if(found) {
					p[k] = Fklists[k][b];
					Qb[Fklists[k][b]] += complement[k];
					iterations[k] = std::make_pair(b,x);
					Vector<int> newL(L.dim()+1);
					for(int y = 0; y < x; y++) newL[y] = L[y];
					newL[x] = Fklists[k][b];
					for(int y = x+1; y < newL.dim(); y++) newL[y] = L[y-1];
					L = newL;
					k++;
					continue;
				}
				//Otherwise step down
				else {
					iterations[k] = std::make_pair(-1,-1);
					k--;
					continue;
				}

			} //End compute cones

		} //End for all bases

		return std::make_pair(rays,cones);
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Takes a bergman fan of a matrix or general matroid and computes the final result, 
	  i.e. adds lineality spaces for missing coloops, leading coordinates and a vertex at the origin.
	  @param int n The total number of elements of the ground set
	  @param Matrix<Rational> bergman_rays The rays of the reduce matroid (without coloops)
	  @param IncidenceMatrix<> cones The cones of the matroid fan.
	  @param Set<int> coloops The element indices that were removed from the ground set, since they were coloops
	  */
	template <typename Addition>
	perl::Object modify_fan(int n, Matrix<Rational> bergman_rays, IncidenceMatrix<> pre_cones, Set<int> coloops) {
		Matrix<Rational> bergman_lineality = Matrix<Rational>(0,n);
		Vector<Integer> weights = ones_vector<Integer>(pre_cones.rows()); 

		//Next, we have to add appropriate zero columns for each coloop and a lineality space
		Matrix<Rational> new_bergman_rays(bergman_rays.rows(),n);
		Matrix<Rational> new_bergman_lin(bergman_lineality.rows() + coloops.size(), n);
		new_bergman_rays.minor(All,~coloops) = bergman_rays;
		Set<int> first_rows = sequence(0,bergman_lineality.rows());
		new_bergman_lin.minor(first_rows,~coloops) = bergman_lineality;
		new_bergman_lin.minor(~first_rows,coloops) = unit_matrix<Rational>(coloops.size());
		bergman_rays = new_bergman_rays;
		bergman_lineality = new_bergman_lin;

		//We have to reduce the lineality space mod (1,..,1) - since we only get
		//lineality space from the coloops, we only have to do this, if everything is a coloop
		if(coloops.size() == n)
			bergman_lineality = bergman_lineality.minor(~scalar2set(0),All);
		

		//Make max or min and add vertex
		bergman_rays *= Addition::orientation();
		bergman_rays = zero_vector<Rational>(bergman_rays.rows()) | bergman_rays;
		bergman_rays /= unit_vector<Rational>(n+1,0);
		bergman_lineality = zero_vector<Rational>(bergman_lineality.rows()) | bergman_lineality;
		IncidenceMatrix<> cones(0,bergman_rays.rows());
		for(int mc = 0; mc < pre_cones.rows(); mc++) {
			cones /= (pre_cones.row(mc) + (bergman_rays.rows()-1));
		}

		perl::Object result(perl::ObjectType::construct<Addition>("Cycle"));
		result.take("PROJECTIVE_VERTICES") << bergman_rays;
		if(bergman_lineality.rows() > 0)
			result.take("LINEALITY_SPACE") << bergman_lineality;
		result.take("MAXIMAL_POLYTOPES") << cones;
		result.take("WEIGHTS") << weights;
		return result;
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Takes a matrix and computes its vector matroid's bergman fan
	  @param Matrix<Rational> m Any rational matrix
	  @return The bergman fan of the given matroid in non-homog. coordinates
	  */
	template <typename Addition>
	perl::Object prepareBergmanMatrix(Matrix<Rational> m) {
		int n = m.cols();
		//First we check for loops - if there is one, the fan is empty - and coloops
		Set<int> coloops;
		int mrank = rank(m);
		for(int c = 0; c < m.cols(); c++) {
			if(m.col(c) == zero_vector<Rational>(m.rows())) {
				return empty_cycle<Addition>(m.cols()-1);
			}
			if(rank(m.minor(All,~scalar2set(c))) < mrank) {
				coloops += c;
			}
		}
		m = m.minor(All,~coloops);
		//Now we make sure that m.rows = rank(m)
		Set<int> rbasis = basis_rows(m);
		m = m.minor(rbasis,All);
		//Now we compute the bergman fan and extract values
		IncidenceMatrix<> I = computeMatrixBases(m);
		std::pair<Matrix<Rational>, IncidenceMatrix<> > computation = bergman_fan(m.cols(),I,1,m);
		return modify_fan<Addition>(n,computation.first, computation.second ,coloops);
	}

	///////////////////////////////////////////////////////////////////////////////////////

	/**
	  @brief Takes a matroid and computes its bergman fan
	  @param int n The number of elements of the ground set
	  @param IncidenceMatrix<> bases The bases of the matroid
	  @return The bergman fan of the given matroid in non-homog. coordinates
	  */
	template <typename Addition>
	perl::Object prepareBergmanMatroid(perl::Object matroid) {
		int n = matroid.give("N_ELEMENTS");
		Array<Set<int> > bases_array = matroid.give("BASES");
		IncidenceMatrix<>	bases = IncidenceMatrix<>(bases_array);
		//First we check for loops and coloops
		Set<int> coloops = matroid.give("COLOOPS");
		Set<int> loops = matroid.give("LOOPS");
		//If it has any loops, the fan is empty
		if(loops.size() > 0) {
			return empty_cycle<Addition>(n-1);
		}
		bases = bases.minor(All,~coloops);
		std::pair<Matrix<Rational>, IncidenceMatrix<> > computation = bergman_fan(n-coloops.size(),bases,0,Matrix<Rational>());
		return modify_fan<Addition>(n,computation.first, computation.second,coloops);
	}
	/*
		void measure1(IncidenceMatrix<> bases, Matrix<Rational> m) {
		for(int b = 0; b < bases.rows(); b++) {
		Vector<int> complement(sequence(0,m.cols()) - bases.row(b));
		computeFkLinear(bases,b,m, complement);	      
		}
		}*/

	// ------------------------- PERL WRAPPERS ---------------------------------------------------

	FunctionTemplate4perl("prepareBergmanMatrix<Addition>(Matrix<Rational>)");
	FunctionTemplate4perl("prepareBergmanMatroid<Addition>(matroid::Matroid)");

}}

