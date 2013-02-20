/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 * 
 * ---
 * Copyright (C) 2012, Simon Hampe <hampe@mathematik.uni-kl.de>
 * 
 * 
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Graph.h"
#include "polymake/linalg.h"
#include "polymake/atint/LoggingPrinter.h"

namespace polymake { namespace atint { 
  
  using namespace atintlog::donotlog;
//   using namespace atintlog::dolog;
//   using namespace atintlog::dotrace;
  
//   int insert_ray(Matrix<Rational> &rays, Vector<Rational> nray) {
//     //Normalize ray
//     for(int c = 0; c < nray.dim(); c++) {
//       if(nray[c] != 0) {
// 	nray /= abs(nray[c]); break;
//       }
//     }
//     
//     //Insert
//     int new_rayindex = -1;
//     for(int oray = 0; oray < rays.rows(); oray++) {
//       if(rays.row(oray) == nray) {
// 	new_rayindex = oray; break;
//       }
//     }
//     if(new_rayindex == -1) {
//       rays /= nray;
//       new_rayindex = rays.rows()-1;
//     }
//     
//     return new_rayindex;
//   }
//   
//   void insert_cone_in_list(Vector<Set<int> > &cones, Vector<Integer> &weights, Set<int> ncone, Integer nweight) {
     //dbgtrace << "Inserting a cone" << endl;
     //dbgtrace << "Cones are " << cones << ", weights are " << weights << " new cone: " << ncone << endl;
//     int ncone_index = -1;
//     for(int c = 0; c < cones.dim(); c++) {
//       Set<int> inter = ncone * cones[c];
//       if( inter.size() == ncone.size() && inter.size() == cones[c].size()) {
// 	ncone_index = c; break;
//       }
//     }
//     if(ncone_index == -1) {
//       cones |= ncone;
//       weights |= nweight;
//     }
//     else {
//       weights[ncone_index] += nweight;
//     }
//     
     //dbgtrace << "New weights " << weights << endl;
//   }
//   
//   ///////////////////////////////////////////////////////////////////////////////////////
//   
//   perl::Object hurwitz6_recession(perl::Object pre_cycle) {
//     
//     //Extract values
//     Matrix<Rational> rays = pre_cycle.give("RAYS");
//     IncidenceMatrix<> cones = pre_cycle.give("MAXIMAL_CONES");
//     Set<int> dir = pre_cycle.give("DIRECTIONAL_RAYS");
//     Vector<Integer> weights = pre_cycle.give("TROPICAL_WEIGHTS");
//     
//     //Group cones according to their combinatorial type
//     Vector<Set<int> > cone_groups;
//     Matrix<Rational> cone_group_representative(0,9);
//     
//     Set<int> leaves_to_remove = sequence(7,2);
//       
//     //SORT CONES
//     for(int mc = 0; mc < cones.rows(); mc++) {
//       //If the recession fan has not the right dimension, ignore the cone
//       if( (cones.row(mc) * dir).size() < 2) continue;
//       
//       perl::Object mccurve = CallPolymakeFunction("rational_curve_from_cone",pre_cycle,8,mc);
//       Vector<Set<int> > mcsets = mccurve.give("SETS");
//       for(int mcs = 0; mcs < mcsets.dim(); mcs++) {
// 	mcsets[mcs] = mcsets[mcs] - leaves_to_remove;
//       }
//       perl::Object curve_to_search("RationalCurve");
// 	curve_to_search.take("N_LEAVES") << 6;
// 	curve_to_search.take("SETS") << mcsets;
// 	curve_to_search.take("COEFFS") << ones_vector<Rational>(mcsets.dim());
//       Vector<Rational> curve_vector = curve_to_search.CallPolymakeMethod("matroid_vector");
// 	
//       bool match_found = 0;
//       for(int cg = 0; cg < cone_groups.dim(); cg++) {
// 	if(curve_vector == cone_group_representative.row(cg)) {
// 	    cone_groups[cg] += mc;
// 	    match_found = 1;
// 	    break;
// 	}
//       }
//       if(!match_found) {
// 	Set<int> newcg;
// 	  newcg += mc;
// 	cone_groups |= newcg;
// 	cone_group_representative /= curve_vector;
//       }
//     }//END sort maximal cones
//     
//     //DEBUG:
// //     cone_groups = cone_groups.slice(scalar2set(0));
      //dbgtrace << "Cones: " << cone_groups << endl;
//     
//     perl::Object ffm = CallPolymakeFunction("forgetful_map", 8,leaves_to_remove);
//     Matrix<Rational> ff_matrix = ffm.give("MATRIX");
//     
     //dbgtrace << "Computing recession fan " << endl;
//     
//     //COMPUTE RECESSION FAN
//     Matrix<Rational> rfan_rays(0, 9);
//     Vector<Set<int> > rfan_cones;
//     Vector<Integer> rfan_weights;
//     
//     for(int cg = 0; cg < cone_groups.dim(); cg++) {
       //dbgtrace << "Subdividing cone group " << cg << endl;
//       Set<int> cones_in_group = cone_groups[cg];
//       Set<int> group_ray_set = 
// 	  accumulate(rows(cones.minor(cones_in_group,All)), operations::add()) * dir;
//       Vector<int> group_ray_list(group_ray_set);
//       
       //dbgtrace << "Computing new rays\n";
//       Matrix<Rational> group_rays = rays.minor( group_ray_set, ~scalar2set(0));
//       group_rays = group_rays * T(ff_matrix);
//       Map<int,int> new_index;
//       for(int gr = 0; gr < group_rays.rows(); gr++) {
// 	new_index[group_ray_list[gr]] = insert_ray(rfan_rays, group_rays.row(gr));
//       }
//       
       //dbgtrace << "Computing new cones\n";
//       Vector<Set<int> > queue;
//       Vector<Integer> wqueue = weights.slice(cones_in_group);
//       for(Entire<Set<int> >::iterator cig = entire(cones_in_group); !cig.at_end(); cig++) {
// 	queue |= attach_operation( cones.row(*cig) * dir, pm::operations::associative_access<Map<int,int>, int>(&new_index));
//       }
//       
       //dbgtrace << "Iterating queue...\n";
//       
//       while(queue.dim() > 0) {
 	//dbgtrace << "Queue: " << queue << endl;
// 	//Extract first element
// 	Set<int> front = queue[0]; 
 	//dbgtrace << "Checking cone " << front << endl;
// 	queue = queue.slice(~scalar2set(0));
// 	Integer front_weight = wqueue[0];
// 	wqueue = wqueue.slice(~scalar2set(0));
// 	//Find badly intersecting cone
// 	bool found_bad_cone = false;
// 	for(int q = 0; q < queue.dim(); q++) {
// 	  if( (front * queue[q]).size() == 0) {
// 	    //Compute intersection
 	    //dbgtrace << "Intersecting with cone " << queue[q] << endl;
// 	    perl::Object c1("polytope::Cone");
// 	      c1.take("RAYS") << rfan_rays.minor(front,All);
// 	    perl::Object c2("polytope::Cone");
// 	      c2.take("RAYS") << rfan_rays.minor(queue[q],All);
// 	    perl::Object intercone = CallPolymakeFunction("polytope::intersection",c1,c2);
// 	    Matrix<Rational> interrays = intercone.give("RAYS");
// 	    if(interrays.rows() > 0) {
// 		int nindex = insert_ray(rfan_rays, interrays.row(0));
 		//dbgtrace << "Intersection in ray  " << nindex << endl;
// 		//If a cone already contains the new ray, just add it to the queue
// 		//otherwise subdivide it
// 		Vector<Set<int> > cones_sub;
// 		  cones_sub |= front; cones_sub |= queue[q];
// 		for(int i = 0; i < 2; i++) {
// 		    Set<int> ssub = cones_sub[i];
// 		    if(ssub.contains(nindex)) {
// 		      queue |= ssub;
// 		      wqueue |= (i == 0? front_weight : wqueue[q]);
// 		    }
// 		    else {
// 		      Vector<int> slist(ssub);
// 		      for(int v = 0; v < 2; v++) {
// 			Set<int> ssubset;
// 			ssubset += nindex; ssubset += slist[v];
// 			queue |= ssubset;
// 			wqueue |= (i == 0? front_weight : wqueue[q]);
// 		      }
// 		    }
// 		}//END subdivide cones and add to queue
// 		//Also remove badly intersecting cone
// 		queue = queue.slice(~scalar2set(q));
// 		wqueue = wqueue.slice(~scalar2set(q));
// 		found_bad_cone = true;
// 		break;
// 	    }//END if intersection is bad
// 	  }//END if no intersection a priori
// 	}//END iterate queue
// 	
// 	//If we arrive here, there was none, we can add the current cone
 	//dbgtrace << "Inserting cone " << front << endl;
// 	if(!found_bad_cone) insert_cone_in_list(rfan_cones, rfan_weights, front, front_weight);
//       }
//       
//       
//     }//END iterate cone groups
//     
//     perl::Object result("WeightedComplex");
//       result.take("RAYS") << rfan_rays;
//       result.take("MAXIMAL_CONES") << rfan_cones;
//       result.take("TROPICAL_WEIGHTS") << rfan_weights;
//     
//     return result;
//     
//   }
//   
//   perl::Object refinem6(perl::Object recfan) {
//     
//     Matrix<Rational> rays = recfan.give("RAYS");
//     IncidenceMatrix<> cones = recfan.give("MAXIMAL_CONES");
//     
//     //Compute M_0,6
//     perl::Object m6 = CallPolymakeFunction("tropical_m0n",6);
//     Matrix<Rational> mrays = m6.give("RAYS");
//     IncidenceMatrix<> mcones = m6.give("MAXIMAL_CONES");
//     
//     //Result variables
//     Matrix<Rational> rrays = rays;
//     Map<int,int> m6tores;
//     for(int mr = 0; mr < mrays.rows(); mr++) {
//       m6tores[mr] = insert_ray(rrays, mrays.row(mr));
//     }
//     Vector<Set<int> > rcones;
//     
//     //Iterate cones of m6
//     for(int mc = 0; mc < mcones.rows(); mc++) {
//       //Find cones contained in this one 
//       perl::Object mccurve = CallPolymakeFunction("rational_curve_from_cone",m6,6,mc);
//       Vector<Set<int> > mcsetlist = mccurve.give("SETS");
       //dbgtrace << "Refining cone " << mc << ": " << mcsetlist <<  endl;
// 	for(int v = 0; v < mcsetlist.dim(); v++) {
// 	    if(!mcsetlist[v].contains(1)) mcsetlist[v] = sequence(1,6) - mcsetlist[v];
// 	}
//       Set<Set<int> > mcsets(mcsetlist);
//       
//       Vector<int> relevant_cones;
//       Vector<int> expected_cones;
//       bool found_interior = false;
//       for(int rc = 0; rc < cones.rows(); rc++) {
// 	perl::Object rccurve = CallPolymakeFunction("rational_curve_from_cone",recfan,6,rc);
// 	Vector<Set<int> > rcsetlist = rccurve.give("SETS");
// 	  for(int w = 0; w < rcsetlist.dim(); w++) {
// 	      if(!rcsetlist[w].contains(1)) rcsetlist[w] = sequence(1,6) - rcsetlist[w];
// 	  }
// 	if( (mcsets * Set<Set<int> >(rcsetlist)).size() == rcsetlist.dim()) {
// 	  relevant_cones |= rc;
// 	  if(rcsetlist.dim() == 3) found_interior = true;
// 	  expected_cones |= (rcsetlist.dim() == 3 ? 2 : 1);
// 	}
// 	
//       }//END iterate recession fan cones
//       
       //dbgtrace << "Relevant cones: " << relevant_cones << endl;
       //dbgtrace << "Expected cones: " << expected_cones << endl;
//       
//       //If there were no interior cones, just add the cone of m6
//       if(!found_interior) {
// 	Set<int> translated_cone = attach_operation(mcones.row(mc), pm::operations::associative_access<Map<int,int>, int>(&m6tores));
// 	rcones |= translated_cone;
//       }
//       else {
// 	Vector<int> found_cones(relevant_cones.dim());
// 	for(int c = 0; c < relevant_cones.dim(); c++) {
// 	    if(found_cones[c] < expected_cones[c]) {
 	      //dbgtrace << "Starting a cone with " << c << endl;
// 	      found_cones[c]++;
// 	      Set<int> newcone = cones.row(relevant_cones[c]);
// 	      //Find first cone to add
// 	      for(int d1 = 0; d1 < relevant_cones.dim(); d1++) {
// 		if( d1 != c && found_cones[d1] < expected_cones[d1] &&
// 		  (newcone * cones.row(relevant_cones[d1])).size()  == newcone.size() -1) {
// 		  //Find second cone to add
// 		  Set<int> coned1 = newcone + cones.row(relevant_cones[d1]);
// 		  bool found_one = false;
// 		  for(int d2 = 0; d2 < relevant_cones.dim(); d2++) {
// 		    if(d1 != d2 && d2 != c && found_cones[d2] < expected_cones[d2] &&
// 		      (coned1 * cones.row(relevant_cones[d2])).size()  == coned1.size() -1) {
 			//dbgtrace << "Adding " << d1 << " and " << d2 << endl;
// 			coned1 += cones.row(relevant_cones[d2]);
// 			found_cones[d1]++;found_cones[d2]++;
// 			rcones |= coned1;
// 			found_one = true;
// 			break;
// 		    }
// 		  }
// 		  if(found_one)  break;
// 		}
// 		
// 	      }
// 	      
// 	    }
// 	}
// //  	if(mc == 9) return perl::Object("WeightedComplex");
//       }
       //dbgtrace << "Refined cones: " << rcones << endl;
//       
//     }//END iterate cones of m6
//     
//     perl::Object result("WeightedComplex");
//       result.take("RAYS") << rrays;
//       result.take("MAXIMAL_CONES") << rcones;
//       result.take("TROPICAL_WEIGHTS") << ones_vector<Integer>(rcones.dim());
//       
//     return result;  
//     
//   }
  
  Vector<Integer> weight_aim(perl::Object mref, perl::Object hurwitz) {
    
    Matrix<Rational> mref_rays = mref.give("RAYS"); 
    IncidenceMatrix<> mref_codim = mref.give("CODIM_1_FACES");
    
    Matrix<Rational> hurwitz_rays = hurwitz.give("RAYS");
    IncidenceMatrix<> hurwitz_cones = hurwitz.give("MAXIMAL_CONES");
    Vector<Integer> hurwitz_weights = hurwitz.give("TROPICAL_WEIGHTS");
    
    //Identify rays
    Map<int,int> htom;
    for(int h = 0; h < hurwitz_rays.rows(); h++) {
      for(int m = 0; m < mref_rays.rows(); m++) {
	if(hurwitz_rays.row(h) == mref_rays.row(m)) {
	    htom[h] = m; break;
	}
      }
    }
    
    //Identify cones
    Map<int,int> cone_htom;
    Vector<Integer> aim(mref_codim.rows());
    for(int hc = 0; hc < hurwitz_cones.rows(); hc++) {
      Set<int> hcone_conv = attach_operation(hurwitz_cones.row(hc), pm::operations::associative_access<Map<int,int>,int>(&htom));
      for(int mc = 0; mc < mref_codim.rows(); mc++) {
	if( (hcone_conv * mref_codim.row(mc)).size() == hcone_conv.size()) {
	  aim[mc] = hurwitz_weights[hc];
	  break;
	}
      }
    }
    
    return aim;
    
    
  }
  
  /*
  Vector<int> findShelling(perl::Object cycle) {
    
    //Extract values
    IncidenceMatrix<> cones = cycle.give("MAXIMAL_CONES");
    IncidenceMatrix<> codim = cycle.give("CODIM_1_FACES");
    IncidenceMatrix<> maxInCodim = cycle.give("CODIM_1_IN_MAXIMAL_CONES");
      maxInCodim = T(maxInCodim);
    
    //Initialize backtrack variables
    Vector<int> currentSequence;
    Vector<Set<int> > triedCones(cones.rows());
    
    while(currentSequence.dim() < cones.rows()) {
      //dbgtrace << "Current sequence: " << currentSequence << endl;
      //dbgtrace << "Tried cones: " << triedCones[currentSequence.dim()] << endl;
      //Find next untried cone
      int tryIndex = 0;
      while( (Set<int>(currentSequence).contains(tryIndex) || triedCones[currentSequence.dim()].contains(tryIndex)) && tryIndex < triedCones.dim()) {
	tryIndex++;
      }
      
      //If we find none, go back a step
      if(tryIndex >= triedCones.dim()) {
	if(currentSequence.dim() > 0) {
	  triedCones[currentSequence.dim()] = Set<int>();
	  currentSequence = currentSequence.slice(sequence(0,currentSequence.dim()-1));
	}
	continue;
      }
      
      //dbgtrace << "Trying to add " << tryIndex << endl;
      
      //Otherwise, check if this works
      triedCones[currentSequence.dim()] += tryIndex;
      
      //Compute union of all previous cones
      Set<int> prevUnion;
      for(int i = 0; i < currentSequence.dim(); i++) {
	prevUnion += cones.row(currentSequence[i]);
      }
      prevUnion = prevUnion * cones.row(tryIndex);
      
      //dbgtrace << "Union of previous cones: " << prevUnion << endl;
      
      //Then see if it is a union of codim 1 faces of tryIndex - cone
      Set<int> tryCodims = maxInCodim.row(tryIndex);
      Set<int> codimUnion;
      for(Entire<Set<int> >::iterator cd = entire(tryCodims); !cd.at_end(); cd++) {
	if( (codim.row(*cd) * prevUnion).size() == codim.row(*cd).size()) {
	    codimUnion += codim.row(*cd);
	}
      }
      
      //dbgtrace << "Union of codim 1 faces: " << codimUnion << endl;
      
      if(codimUnion.size() == prevUnion.size()) {
	currentSequence |= tryIndex;
      }
    }
    
    return currentSequence;
    
    
  }*/
  Map<int,int> cf_loc(perl::Object mref, perl::Object hurwitz) {
    Matrix<Rational> mref_rays = mref.give("RAYS"); 
    IncidenceMatrix<> mref_max = mref.give("MAXIMAL_CONES");
    IncidenceMatrix<> mref_codim = mref.give("CODIM_1_FACES");
    IncidenceMatrix<> mref_codim_in_max = mref.give("CODIM_1_IN_MAXIMAL_CONES");
    
    Matrix<Rational> hurwitz_rays = hurwitz.give("RAYS");
    IncidenceMatrix<> hurwitz_cones = hurwitz.give("MAXIMAL_CONES");
    Vector<Integer> hurwitz_weights = hurwitz.give("TROPICAL_WEIGHTS");
    
    //Identify rays
    Map<int,int> htom;
    for(int h = 0; h < hurwitz_rays.rows(); h++) {
      for(int m = 0; m < mref_rays.rows(); m++) {
	if(hurwitz_rays.row(h) == mref_rays.row(m)) {
	    htom[h] = m; break;
	}
      }
    }
    
    //Identify cones
    Map<int,int> cone_htom;
    Vector<Integer> aim(mref_codim.rows());
    for(int hc = 0; hc < hurwitz_cones.rows(); hc++) {
      Set<int> hcone_conv = attach_operation(hurwitz_cones.row(hc), pm::operations::associative_access<Map<int,int>,int>(&htom));
      for(int mc = 0; mc < mref_codim.rows(); mc++) {
	if( (hcone_conv * mref_codim.row(mc)).size() == hcone_conv.size()) {
	  aim[mc] = hurwitz_weights[hc];
	  cone_htom[hc] = mc;
	  break;
	}
      }
    }
    
    return cone_htom;
    
    
    
//     pm::cout << "Weight aim: " << aim << endl;
//     
//     //Now do a local computation for each codimension one cone
//     Matrix<Rational> global_cutting_eq(0, mref_rays.rows()+1);
//     for(int mc = 0; mc < mref_codim.rows(); mc++) {
//       pm::cout << "Codim 1 face " << (mc+1) << " of " << mref_codim.rows() << endl;
//       Set<int> surrounding_rays = accumulate(rows(mref_max.minor(mref_codim_in_max.row(mc),All)),operations::add());
//       perl::Object mloc = CallPolymakeFunction("local_codim_1",mref,mc);
//       Vector<Integer> loc_aim; loc_aim |= aim[mc];
//       Matrix<Rational> local_cutting = 
// 	CallPolymakeFunction("cutting_functions",mloc, loc_aim);
// 	local_cutting = null_space(local_cutting);
//       Matrix<Rational> local_transform(local_cutting.rows(),mref_rays.rows()+1);
// 	local_transform.minor(All, surrounding_rays) = local_cutting.minor(All,sequence(0,local_cutting.cols()-1));
// 	local_transform.col(local_transform.cols()-1) = local_cutting.col(local_cutting.cols()-1);
//       global_cutting_eq /= local_transform;
//     }
//     
//     return null_space(global_cutting_eq);
  }
  
  
  
  Matrix<Rational> graph_type_solutions(perl::Object subdivision, perl::Object cycle, Vector<int> degree, int min, int max) {
    
    Vector<Integer> aim = weight_aim(subdivision, cycle);
    pm::cout << "Computing cutting functions" << endl;
    Matrix<Rational> cutting_functions = CallPolymakeFunction("cutting_functions", subdivision, aim);
      cutting_functions = null_space(cutting_functions);
    //Sort rays of subdivision by equivalence class
    Vector<Set<int> > eq_classes;
    Matrix<Rational> rays = subdivision.give("RAYS");
    perl::ListResult types = ListCallPolymakeFunction("rational_curve_list_from_moduli",rays);
    Vector<Rational> linear_comb;
    Map<int,int> type_class;
    for(int t = 0; t < types.size(); t++) {
      perl::Object this_type = types[t];
      int cl = -1;
      for(int c = 0; c < eq_classes.dim(); c++) {
	perl::Object comp_type = types[ *(eq_classes[c].begin())];
	Graph<> g1 = this_type.give("GRAPH.ADJACENCY");
	Graph<> g2 = comp_type.give("GRAPH.ADJACENCY");
	if( CallPolymakeFunction("graph::isomorphic",g1,g2)) {
	    cl = c; break;
	}
      }
      if(cl == -1) {
	Set<int> nset; nset += t;
	eq_classes |= nset;
	type_class[t] = eq_classes.dim()-1;
      }
      else {
	eq_classes[cl] += t;
	type_class[t] = cl;
      }
      //Also compute linear combination of absolute values of sums of splits
      Vector<Set<int> > tsets = this_type.give("SETS");
      Vector<Rational> tcoeffs = this_type.give("COEFFS");
      Rational s(0);
      for(int tc = 0; tc < tcoeffs.dim(); tc++) {
	int csum = 0;
	Set<int> cs = tsets[tc];
	for(Entire<Set<int> >::iterator i = entire(cs); !i.at_end(); i++) {
	  csum += degree[*i - 1];
	}
	s += tcoeffs[tc] * abs(csum);
      }
      linear_comb |= s;
    }//END sort rays by type
    
    //Transform cutting functions into polyhedron equation
    cutting_functions = 
      cutting_functions.col(cutting_functions.cols()-1) | cutting_functions.minor(All,sequence(0,cutting_functions.cols()-1)) | Matrix<Rational>(cutting_functions.rows(), eq_classes.dim());

    //Create equations for equivalence classes: phi(ray_i) (= i-th entry) = linear_comb[i] * c_corr_class
    for(int r = 0; r < types.size(); r++) {
	Vector<Rational> r_eq = zero_vector<Rational>(cutting_functions.cols());
	r_eq[r+1] = 1;
	r_eq[types.size() + type_class[r] + 1] = -linear_comb[r];
	cutting_functions /= r_eq;
    }
    
    //Now create min/max inequalities
    Matrix<Rational> ineq_min = unit_matrix<Rational>(types.size());
    ineq_min = ( - min * ones_vector<Rational>(types.size())) | ineq_min;
    ineq_min = ineq_min | Matrix<Rational>(ineq_min.rows(),eq_classes.dim());
    
    Matrix<Rational> ineq_max = - unit_matrix<Rational>(types.size());
    ineq_max = ( max * ones_vector<Rational>(types.size())) | ineq_max;
    ineq_max = ineq_max | Matrix<Rational>(ineq_max.rows(), eq_classes.dim());
  
    //Print out classes 
    pm::cout << "Number of classes: " << eq_classes.dim() << endl;
    for(int c = 0; c < eq_classes.dim(); c++) {
      perl::Object t = types[ *(eq_classes[c].begin())];
      std::string ts = t.CallPolymakeMethod("to_string");
      pm::cout << "Class " << (c+1) << ": " << ts << endl;
    }
    
    pm::cout << "Computing null space" << endl;
    Matrix<Rational> result = null_space(cutting_functions);
    for(int rr = 0; rr < result.rows(); rr++) {
      if(result(rr,0) != 0) {
	result.row(rr) *= (1 / result(rr,0));
      }
    }
    
    return result;
    
//     perl::Object poly("polytope::Polytope");
//       poly.take("EQUATIONS") << cutting_functions;
//       poly.take("INEQUALITIES") << //(ineq_max / ineq_min);
//     return poly;
    
  }
  
  
  // ------------------------- PERL WRAPPERS ---------------------------------------------------
  
//   Function4perl(&findShelling,"fshell(WeightedComplex)");
//     Function4perl(&hurwitz6_recession, "h6r(WeightedComplex)");
//     Function4perl(&refinem6,"rm6(WeightedComplex)");
    Function4perl(&weight_aim,"wa6(WeightedComplex,WeightedComplex)");
    Function4perl(&graph_type_solutions, "gts(WeightedComplex, WeightedComplex, Vector<Int>, $,$)");
    Function4perl(&cf_loc,"cfl(WeightedComplex,WeightedComplex)");
}}