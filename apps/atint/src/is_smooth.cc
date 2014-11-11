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
 * Copyright (C) 2014, Dennis Diefenbach <diefenba@mathematik.uni-kl.de>
 * 
 * This file contains a function to test smoothness of tropical fans (at least for dimension 0,1,2
 * and codimension 0,1).
 */

#include "polymake/client.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/atint/LoggingPrinter.h"
#include "polymake/list"
#include "polymake/Vector.h"
#include "polymake/PowerSet.h"
#include "polymake/atint/lll.h"
#include "polymake/atint/converters.h"

namespace polymake {
    namespace atint {
      
        using namespace atintlog::donotlog;
	//using namespace atintlog::dolog;
      //   using namespace atintlog::dotrace;
      
        //This is the type returned by many functions, it is a tripel consisting of:
            //-An int which is 0 if the curve is not smooth, 1 if it is smooth, 2 if it is not clear
            //-if s=1 the "matroid" contains a matroid such that the corresponding matroid fan is isomrphic to the input fan
            //-"map" is a Z-isomorphisem between the input fan and the matroid fan
        struct result{
            int is_smooth;
            perl::Object matroid;
            Matrix<int> map;
        };
        
        //This functions takes the lineality space of the fan, searches a Z-isomorphisem mapping it to the last coordinates, projects the last coordinates down and returnes the projection of the fan and by reference the Z-isomorphisem
        perl::Object cut_out_lineality_space(perl::Object f, Matrix <int> &map_lineality_space);
        
        //If the linear spann U of the rays of the input fan is not the hole space, the function can find a Z-isomorphisem h that maps U to the first coordinates and gives back the image of the fan under h and by reference h
        perl::Object full_dimensional(perl::Object f, Matrix <int> &map_full_dimensional);
        
        //Takes as an input a fan of dimension 1 and returnes 0 if they are NOT smooth and 1 if they are smooth
        result smooth_dim1(perl::Object f);
        
        //Takes as an input a fan of dimension 2 and returnes 0 if they are NOT smooth and 1 if they are smooth
        result smooth_dim2(perl::Object f);
        
        //Checks if the graph is the combinatorial graph of a matroid
        result is_combiantorial_right(perl::Object f, Vector< Vector<int> > combinatorics, int N_Nodes, int N_Edges, int Fan_Ambient_Dimension);
        
        //Checks if the Flats satisfies the axiom of a matroid
        int check_matroid(Vector< Vector<int> > &combinatorics,  Vector < Set <int> > &Flats, int Fan_Ambient_Dimension);
        
        //Takes as an input a fan of codimension 1 and returnes 0 if they are NOT smooth and 1 if they are smooth
        result smooth_codim1(perl::Object f);
        
        
        //This functions takes a fan and checks if it is smooth or not, i.e. it checks if there exists a matroid fan such that the corresponing matroid fan is Z-isomorphic to the given fan
        perl::ListReturn is_smooth(perl::Object f){
	    f = f.CallPolymakeMethod("dehomogenize");	  
	  
            result R;
            R.is_smooth=2;
            //Check that all tropical weights are 1
            Vector<int> Tropical_Weights=f.give("TROPICAL_WEIGHTS");
            for (int i=0; i<Tropical_Weights.size(); i++) {
                if (Tropical_Weights[i]!=1) {
                    R.is_smooth=0;
                }
            }
            
            //Find the coarsest subdivision of the fan f
            perl::Object g;
            if (R.is_smooth!=0) {
            
                try {
                    g=CallPolymakeFunction("coarsen",f,1);
                }
                catch (std::runtime_error) {
                    //The coarsest subdivision does not exist
                    R.is_smooth=0;
                }
            }
            
            if (R.is_smooth!=0) {

                //Restrict to the case where the fan has no lineality space
                int Lineality_Dim=g.give("LINEALITY_DIM");
                Matrix <int> map_lineality_space;
                if (Lineality_Dim>0) {
                    g=cut_out_lineality_space(g, map_lineality_space);
                }
                
                //Determine some elementary properties of the fan f
                int Fan_Dim= g.give("FAN_DIM");
                int N_Rays=g.give("N_RAYS");
                Matrix<Rational> Rays=g.give("RAYS");
                int Linear_Spann_Dim=rank(Rays);
                
                R.is_smooth=2; //s=0 the curve is NOT smooth, s=1 the curve is smooth, s=2 we cannot say if it is smooth or not
                
                //Check if the number of rays is not to big, i.e. upper bounds for the number of rays in the coarsets subdivision of the fan are known, if they are not satisfied the fan is NOT smooth
                if (Fan_Dim==1){
                    if (N_Rays>Linear_Spann_Dim+1){
                        R.is_smooth=0;
                    }
                }
                else if (Fan_Dim==2){
                    if ((Linear_Spann_Dim+1)<=9) { //Upper bound for dimension 1 cells found experimentally
                        int bound[7]={0,4,6,10,14,16,21};
                        if (N_Rays>bound[(Linear_Spann_Dim+1)-3]) {
                             R.is_smooth=0;
                        }
                    }
                    else if (N_Rays>((Linear_Spann_Dim+1)*(Linear_Spann_Dim+1))){ //Theoretical upper bound
                        R.is_smooth=0;
                    }
                }
                else if (Fan_Dim==Linear_Spann_Dim-1){
                    if (N_Rays>Linear_Spann_Dim+1){
                        R.is_smooth=0;
                    }
                }
                
                
                
                //If the above criterions did not suffice to determine if the curve is smooth or not the algorithem enters here
                Matrix <int> map_full_dimensional;
                if (R.is_smooth==2){
                    //Restrict to the case where the fan is full dimensional
                    int Fan_Ambient_Dimension= g.give("FAN_AMBIENT_DIM");
                    int codim_linear_spann=Fan_Ambient_Dimension-Linear_Spann_Dim;
                    if (codim_linear_spann!=0) {
                        g=full_dimensional(g, map_full_dimensional);
                        Fan_Ambient_Dimension= g.give("FAN_AMBIENT_DIM");
                    }
                    if (Fan_Dim==0){
                        Array< Set<int> > Bases(1);
                        Bases[0]=Bases[0]+0;
                        perl::Object M("matroid::Matroid");
                        M.take("BASES") << Bases;
                        M.take("N_ELEMENTS") << 1;
                        R.is_smooth=1;
                        R.matroid=M;
                    }
                    else if (Fan_Dim==1){
                        R=smooth_dim1(g);
                    }
                    else if (Fan_Dim==Fan_Ambient_Dimension-1){
                        R=smooth_codim1(g);
                    }
                    else if (Fan_Dim==2){
                        R=smooth_dim2(g);
                    }
                    else if (Fan_Dim==Fan_Ambient_Dimension){
                        R.is_smooth=0;
                    }
                    else{
                        R.is_smooth=2;
                    }
                    if (R.is_smooth==1) {
                        Array< Set<int> > Bases=(R.matroid).give("BASES");
                        int N_Elements=(R.matroid).give("N_ELEMENTS");
                        if (codim_linear_spann!=0) {
                            //Compute the bases
                            for (int i=0; i<codim_linear_spann; i++) {
                                for (int j=0; j<Bases.size(); j++) {
                                    if (Bases[j].contains(N_Elements-1)) {
                                        Set<int> NewBase=Bases[j]-(N_Elements-1);
                                        NewBase=NewBase+(N_Elements);
                                        Bases.resize(Bases.size()+1,NewBase);
                                    }
                                }
                                N_Elements=N_Elements+1;
                            }
                            //Compute the map
                            Matrix <int> temp((R.map).rows()+codim_linear_spann,(R.map).cols()+codim_linear_spann);
                            for (int i=0; i<(R.map).rows(); i++) {
                                for (int j=0; j<(R.map).rows(); j++) {
                                    temp[i][j]=R.map[i][j];
                                }
                            }
                            for (int i=0; i<codim_linear_spann; i++) {
                                for (int j=0; j<codim_linear_spann; j++) {
                                    if (i==j) {
                                        temp[i+(R.map).rows()][j+(R.map).rows()]=1;
                                    }
                                    else{
                                        temp[i+(R.map).rows()][j+(R.map).rows()]=0;
                                    }
                                }
                            }
                            R.map=map_full_dimensional*temp;
                        }
                        if (Lineality_Dim>0) {
                            //Compute the matroid
                            for (int i=0; i<Lineality_Dim; i++) {
                                for (int j=0; j<Bases.size(); j++) {
                                    if (Bases[j].contains(N_Elements-1)){
                                        Bases[j]=Bases[j]+(N_Elements);
                                    }
                                    else{
                                        Bases[j]=Bases[j]+(N_Elements-1);
                                    }
                                }
                                N_Elements=N_Elements+1;
                            }
                            //Compute the map
                            Matrix <int> temp((R.map).rows()+Lineality_Dim,(R.map).cols()+Lineality_Dim);
                            for (int i=0; i<(R.map).rows(); i++) {
                                for (int j=0; j<(R.map).rows(); j++) {
                                    temp[i][j]=R.map[i][j];
                                }
                            }
                            for (int i=0; i<Lineality_Dim; i++) {
                                for (int j=0; j<Lineality_Dim; j++) {
                                    if (i==j) {
                                        temp[i+(R.map).rows()][j+(R.map).rows()]=1;
                                    }
                                    else{
                                        temp[i+(R.map).rows()][j+(R.map).rows()]=0;
                                    }
                                }
                            }
                            R.map=map_lineality_space*temp;
                        }
                        
                        (R.matroid).take("BASES") << Bases;
                        (R.matroid).take("N_ELEMENTS") << N_Elements;
                    }
                }
            }
//             if (R.is_smooth==0) {
//                 pm::cout<<"The fan is not smooth :)"<<endl;
//             }
//             else if(R.is_smooth==1){
//                 pm::cout<<"The fan is smooth :)"<<endl;
//             }
//             else{
//                 pm::cout<<"I'm sorry, but I do not know if the fan is smooth or not :("<<endl;
//             }
            perl::ListReturn S;
            S<<R.is_smooth;
            if (R.is_smooth==1) {
                S<<R.matroid;
                S<<R.map;
            }
            return S;
        }
        
        
        perl::Object cut_out_lineality_space(perl::Object f, Matrix <int> &map_lineality_space){
            //A basis of the lineality space is computed
            Matrix<Rational> Lineality_Space=f.give("LINEALITY_SPACE");
            int Fan_Ambient_Dimension=f.give("CMPLX_AMBIENT_DIM");
            int Lineality_Dim=f.give("LINEALITY_DIM");
            Matrix<Rational> Matrix_Description=null_space (Lineality_Space);
            
            //If the lineality space is the hole space
            if (Matrix_Description.rows()==0) {
                Matrix_Description.resize(1,Fan_Ambient_Dimension);
                for (int i=0; i<Fan_Ambient_Dimension; i++) {
                    Matrix_Description[0][1]=0;
                }
            }
            
            //Make the matrix integer without changing the kernel
            Matrix<Integer> Matrix_Description_Integer(Matrix_Description.rows(), Matrix_Description.cols());
            for (int i=0; i<Matrix_Description.rows(); i++) {
                Integer LCM_Denominator=1;
                for (int j=0; j<Matrix_Description.cols(); j++){
                    Integer Denominator=denominator(Matrix_Description[i][j]);
                    LCM_Denominator=lcm(LCM_Denominator,Denominator);
                }
                for (int j=0; j<Matrix_Description.cols(); j++){
                    Matrix_Description_Integer[i][j]=Matrix_Description[i][j]*LCM_Denominator;
                }
            }
            
            //The hermit normal form of the Matrix description is computed
            Matrix<Integer> tfmatrix;
            Matrix<Integer> HNF;
            int kernelDimension;
            HNF=lllHNF(T(Matrix_Description_Integer), tfmatrix, kernelDimension);
            tfmatrix=T(tfmatrix); //contains basis of Lineality_space extended to a bases of hole space and that is a basis of the lattice
            
            //A linear map is computed that maps the lineality space with dimension k to the last k coordinates
            map_lineality_space=tfmatrix;
            tfmatrix=inv(tfmatrix);
            
            //Apply the map to the rays
            Matrix<Rational> Rays=f.give("RAYS");
            Rays=T(Rays);
            Rays=tfmatrix*Rays;
            
            //Project the lineality space down and compute the remaining polyhedral fan
            Matrix <Rational> New_Rays(Rays.rows()-Lineality_Dim, Rays.cols());
            for (int i=0; i<Rays.rows()-Lineality_Dim; i++) {
                for (int j=0; j<Rays.cols(); j++) {
                    New_Rays[i][j]=Rays[i][j];
                }
            }

            New_Rays=T(New_Rays);

            IncidenceMatrix<> Maximal_ConesInc =f.give("MAXIMAL_CONES");
            Vector<Set<int> > Maximal_Cones = incmatrixToVector(Maximal_ConesInc);
	    perl::Object g("fan::PolyhedralFan");
            g.take("RAYS") << New_Rays;
            g.take("MAXIMAL_CONES") << Maximal_Cones;
            g.take("FAN_AMBIENT_DIM") << (Fan_Ambient_Dimension-Lineality_Dim);
            return g;
        }
        
        perl::Object full_dimensional(perl::Object f, Matrix <int> &map_full_dimensional){
            //A basis of the space generated by the rays of f is computed
            Matrix<Rational> Rays=f.give("RAYS");
            Set <int> Basis;
            Basis= basis_rows(Rays);
            
            //The basis is written to the matrix Basis_Subspace
            Matrix<Rational> Basis_Subspace(Basis.size(),Rays.cols());
            int k=0;
            for (Entire< Set<int> >::iterator it=entire(Basis); !it.at_end(); it++) {
                for (int i=0; i<Rays.cols(); i++) {
                    Basis_Subspace[k][i]=Rays[*it][i];
                }
                k++;
            }
           
            //A description of the subspace as the kernel of a map is searched
            Matrix<Rational> Matrix_Description=null_space (Basis_Subspace);
            
            //Make the entries of Matrix_Description integers without changing the kernel
            Matrix<Integer> Matrix_Description_Integer(Matrix_Description.rows(), Matrix_Description.cols());
            for (int i=0; i<Matrix_Description.rows(); i++) {
                Integer LCM_Denominator=1;
                for (int j=0; j<Matrix_Description.cols(); j++){
                    Integer Denominator=denominator(Matrix_Description[i][j]);
                    LCM_Denominator=lcm(LCM_Denominator,Denominator);
                }
                for (int j=0; j<Matrix_Description.cols(); j++){
                    Matrix_Description_Integer[i][j]=Matrix_Description[i][j]*LCM_Denominator;
                }
            }
            
            
            //The hermit normal form of the transposed is computed
            Matrix<Integer> tfmatrix;
            Matrix<Integer> HNF;
            int kernelDimension;
            HNF=lllHNF(T(Matrix_Description_Integer), tfmatrix, kernelDimension);
            
            tfmatrix=T(tfmatrix); //contains basis of the subspace extended to a bases of hole space and that is a basis of the lattice
            
            //NOTE: the above matrix maps the vectors of the standard bases to the basis of the subspace extended to a bases of the hole space, the matrix is a Z-isomorphisem
            
            //The order of the columns is permuted such that the vecotrs e_{1}, ..., e_{k} are mapped to a basis of the subspace by tfmatrix
            Matrix<Integer> tfmatrix_temp(tfmatrix.rows(), tfmatrix.cols());
            for (int i=0; i<tfmatrix_temp.rows(); i++) {
                for (int j=0; j<k; j++) {
                    tfmatrix_temp[i][j]=tfmatrix[i][tfmatrix.cols()-k+j];
                }
            }
            for (int i=0; i<tfmatrix_temp.rows(); i++) {
                for (int j=k; j<tfmatrix.cols(); j++) {
                    tfmatrix_temp[i][j]=tfmatrix[i][j-k];
                }
            }
           
            //The matrix represents a Z-isomorphisem mapping the subspace of dimension k to the vectors e_{1},...,e_{k}
            if (Rays.cols()==0) {
                int Fan_Ambient_Dimension=f.give("FAN_AMBIENT_DIM");
                map_full_dimensional.resize(Fan_Ambient_Dimension,Fan_Ambient_Dimension);
                for (int i=0; i<Fan_Ambient_Dimension; i++) {
                    for (int j=0; j<Fan_Ambient_Dimension; j++) {
                        if (i==j) {
                            map_full_dimensional[i][j]=1;
                        }
                        else{
                            map_full_dimensional[i][j]=0;
                        }
                        
                    }
                }
            }
            else{
                map_full_dimensional=tfmatrix_temp;
            }
            tfmatrix=inv(tfmatrix_temp);
            
            //Computes the immage of the rays under the Z-isomorphism constructed above
            Rays=T(Rays);
            Rays=tfmatrix*Rays;
            
            //Ignore the last coordinates and compute the associated fan
            Matrix <Rational> New_Rays(k, Rays.cols());
            for (int i=0; i<k; i++) {
                for (int j=0; j<Rays.cols(); j++) {
                    New_Rays[i][j]=Rays[i][j];
                }
            }
            New_Rays=T(New_Rays);
            IncidenceMatrix<> Maximal_ConesInc =f.give("MAXIMAL_CONES");
            Vector<Set<int> > Maximal_Cones = incmatrixToVector(Maximal_ConesInc);
            perl::Object g("fan::PolyhedralFan");
            g.take("RAYS") << New_Rays;
            g.take("MAXIMAL_CONES") << Maximal_Cones;
            return g;
        }
        
        result smooth_dim1(perl::Object f){
            result R;
            Matrix<Rational> Rays_Rational=f.give("RAYS");
            //Search the correpsonding integer primitive vectors
            Matrix<Integer> Rays(Rays_Rational.rows(), Rays_Rational.cols());
            for (int i=0; i<Rays.rows(); i++) {
                Integer LCM_Denominator=1;
                for (int j=0; j<Rays.cols(); j++){
                    Integer Denominator=denominator(Rays_Rational[i][j]);
                    LCM_Denominator=lcm(LCM_Denominator,Denominator);
                }
                for (int j=0; j<Rays.cols(); j++){
                    Rays_Rational[i][j]=Rays_Rational[i][j]*LCM_Denominator;
                }
                Integer GCD_Numerator=1;
                for (int j=0; j<Rays.cols(); j++){
                    GCD_Numerator=gcd(numerator(Rays_Rational[i][j]),GCD_Numerator);
                }
                for (int j=0; j<Rays.cols(); j++){
                    Rays[i][j]=Rays_Rational[i][j]/GCD_Numerator;
                }
            }
            Rays=T(Rays);
            
            Rays.resize(Rays.rows(), Rays.cols()-1);
            
            //Compute the Hermit normal form of the Rays
            Matrix<Integer> tfmatrix;
            int kernelDimension;
            Matrix<Integer> HNF=lllHNF(Rays, tfmatrix, kernelDimension);
            
            //Compute the determinant of HNF
            Integer determinant=det(HNF);
            if (determinant==1 or determinant==-1) {
                Array< Set<int> > Bases((Rays.cols()*(Rays.cols()+1))/2);
                int k=0;
                for (int i=0; i<Rays.cols()+1; i++) {
                    for (int j=i+1; j<Rays.cols()+1; j++) {
                        Bases[k]=Bases[k]+i+j;
                        k++;
                    }
                }
                
                
                perl::Object M("matroid::Matroid");
                M.take("BASES") << Bases;
                M.take("N_ELEMENTS") << Rays.cols()+1;
                
                //Prepare the output
                R.is_smooth=1;
                R.matroid=M;
                R.map=inv(tfmatrix);
                return R;
            }
            else{
                R.is_smooth=0;
                return R;
            }
        }
        
        result smooth_dim2(perl::Object f){
            result R;
            R.is_smooth=0;
            
            //Some elementary properties are extracted from f
            IncidenceMatrix<> MaximalCones = f.give("MAXIMAL_CONES");
            int Fan_Ambient_Dimension= f.give("FAN_AMBIENT_DIM");
            
            //The combinatorial graph associated to the fan is computed and some elementary properties are extracted
            perl::Object graph_from_edges = CallPolymakeFunction("graph::graph_from_edges",MaximalCones);
            int N_Nodes=graph_from_edges.give("N_NODES");
            int N_Edges=graph_from_edges.give("N_EDGES");
            
            
            
            //The combinatroics matrix is initialized, it has the following structure:
                //One row for each node of the combinatorial graph and one row for each edge of the combiantorial graph.
                //If the row corresponds to a node, in COLUMN 1 is written which one (by enumerating it from 0 to N_Node-1), and COLUMN 2 is zero. COLUMN 3 contains 1 if the node corresponds to a flat of rank 1 and 2 if it corresponds to a flat of rank 2.
                //If the row corresponds to an edge, COLUMN1 and COLUMN2 indicates which nodes they are connecting. COLUMN 3 is 0 if the edge is not split, 1 if it is split and the splitting point corresponds to a flat of rank 1 and 2 if it is split and the splitting point corrsponds to a flat of rank 2.
            
            Vector <Vector <int> > combinatorics((N_Nodes+N_Edges), Vector<int>(3,0));
            for (int i=0; i<N_Nodes; i++){
                combinatorics[i][0]=i;
                combinatorics[i][1]=0;
                combinatorics[i][2]=0;
            }
            Array< Set<int> > Edges=graph_from_edges.CallPolymakeMethod("EDGES");
            for (int i=0; i<N_Edges; i++){
                int j=0;
                for (Entire< Set<int> >::iterator it=entire(Edges[i]); !it.at_end(); it++){
                    combinatorics[(N_Nodes+i)][j]=*it;
                    j++;
                }
                combinatorics[(N_Nodes+i)][2]=0;
            }
           
            //Tries to distribute the flats in the combinatorial graph such that they fullfill the flat axioms and the corrisponding matroid has a matroid fan Z-isomorphic to the input fan
            bool stop=false; //stop=true when a matroid was found such that the corresponding matroid fan is Z-isomorphic to the input fan
            int k=0;//k indicates the row of the matrix combinatorics that is analised
            while (stop==false) {
                if (k==-1) { //the algorithem has checked all possibilities
                    stop=true;
                }
                else if (k==N_Nodes+N_Edges){ //the algorithem has found a possible distribution of flats
                    //Check if it has the right number of rank 1 flats, since the fan is full dimensional there must be Fan_Ambient_Dimension+1 rank 1 flats
                    int n=0;
                    for (int i=0; i<N_Nodes+N_Edges; i++) {
                        if (combinatorics[i][2]==1){
                            n=n+1;
                        }
                    }
                    if (n==(Fan_Ambient_Dimension+1)){
                        //The functions checks if the combinatorics found corresponds to the flat combinatorics of a matroid. If it corresponds to a matroid and there is a Z-isomorphisem that maps the bergman fan of the matroid to the given one then the functions returns R.is_smooth=1, R.matroid the matroid, R.map the Z-isomorphisem. Otherways R.is_smooth=0 is returned
                        R=is_combiantorial_right(f, combinatorics, N_Nodes, N_Edges, Fan_Ambient_Dimension);
                        if (R.is_smooth==1) {
                            return R;
                        }
                    }
                    k=k-1;
                }
                else{
                    if (k<N_Nodes){ //A node can only be a rank 1 or rank 2 flat
                        if(combinatorics[k][2]<=1){
                            combinatorics[k][2]=combinatorics[k][2]+1;
                            k=k+1;
                        }
                        else{
                            combinatorics[k][2]=0;
                            k=k-1;
                        }
                    }
                    else{ //An edge can be split or not, if it is split the splitting point can be a flat of rank 1 or 2
                        if(combinatorics[k][2]<=2){
                            //Let node 1 and node 2 be the nodes connected by the edge. If the edge is not split one node must be correspond to a rank 1 flat, the other to a rank 2 flat, i.e rank1 --- rank2 or rank 2 ---- rank 1. If the edge is split only the possibilities rank 1 --- rank 2 ----rank1 and rank 2 ---- rank 1 ---- rank 2 are possible
                            if (
                                (combinatorics[k][2]+1==1 and combinatorics[combinatorics[k][0]][2]==2 and combinatorics[combinatorics[k][1]][2]==2)
                                or  (combinatorics[k][2]+1==2 and combinatorics[combinatorics[k][0]][2]==1 and combinatorics[combinatorics[k][1]][2]==1
                                     or  (combinatorics[k][2]+1==3 and ((combinatorics[combinatorics[k][0]][2]==1 and combinatorics[combinatorics[k][1]][2]==2) or (combinatorics[combinatorics[k][0]][2]==2 and  combinatorics[combinatorics[k][1]][2]==1))
                                          ))){
                                         combinatorics[k][2]=combinatorics[k][2]+1;
                                         k=k+1;
                                     }
                            else{
                                combinatorics[k][2]=combinatorics[k][2]+1;
                            }
                            
                        }
                        else{
                            combinatorics[k][2]=0;
                            k=k-1;
                        }
                    }
                }
            }
            return R;
        }
        
        result is_combiantorial_right(perl::Object f, Vector< Vector<int> > combinatorics, int N_Nodes, int N_Edges, int Fan_Ambient_Dimension){
            result R;
            R.is_smooth=0;
            
            //Construct the flats using th distribution of the flats in the combinatorial graph
            Set <int> empty;
            Vector < Set <int> > Flats(N_Nodes+N_Edges, empty);
            
            int n=0;
            for (int i=0; i<N_Nodes+N_Edges; i++){
                if (combinatorics[i][2]==1) {
                    //Insert succsessivly the numbers 1,...,Fan_Ambient_Dimension+1 into the rank 1 flats
                    Flats[i].insert(n);
                    //Reconstruct the rank 2 flats using the rank 1 flats
                    if (i<N_Nodes) {
                        for (int j=N_Nodes; j<N_Nodes+N_Edges; j++){
                            if (combinatorics[j][2]==2 and (combinatorics[j][0]==combinatorics[i][0] or combinatorics[j][1]==combinatorics[i][0])) {
                                Flats[j].insert(n);
                            }
                            if (combinatorics[j][2]==3 and combinatorics[j][0]==combinatorics[i][0]) {
                                Flats[combinatorics[j][1]].insert(n);
                            }
                            if (combinatorics[j][2]==3 and combinatorics[j][1]==combinatorics[i][0]) {
                                Flats[combinatorics[j][0]].insert(n);
                            }
                        }
                    }
                    else{
                        for (int j=0; j<N_Nodes; j++){
                            if (combinatorics[j][2]==2 and (combinatorics[j][0]==combinatorics[i][0] or combinatorics[j][0]==combinatorics[i][1])) {
                                Flats[j].insert(n);
                            }
                        }
                    }
                    n=n+1;
                }
            }
            
            
            //Check that all rank 2 flats are different
            bool rank2_different=true;
            for (int i=0; i<N_Nodes+N_Edges; i++){
                if (combinatorics[i][2]==2){
                    for (int j=0; j<N_Nodes+N_Edges; j++){
                        if (i!=j and combinatorics[j][2]==2 and Flats[i]==Flats[j]) {
                            rank2_different=false;
                        }
                    }
                    
                }
            }
            
            if (rank2_different==true) {
                
                //Check if the flats fullfill the flat axioms
                //The functions returns true if the flats fullfill the flat axioms, false otherwise
                bool test=check_matroid(combinatorics, Flats, Fan_Ambient_Dimension);
                if (test==true){
                    //CHECK IF THE FAN CORRESPONDING TO THE MATROID ABOVE IS Z-ISOMORPHIC TO THE GIVEN ONE
                    //Construct the Rays corresponding to the Flats
                    Matrix<Integer> Rays_Matroid_Fan(Fan_Ambient_Dimension,N_Nodes);
                    for (int j=0; j<N_Nodes; j++){
                        for (int i=0; i<Fan_Ambient_Dimension; i++){
                            if (Flats[j].contains(i)){ 
                                Rays_Matroid_Fan[i][j]=1;
                            }
                        }
                        if (Flats[j].contains(Fan_Ambient_Dimension)) {
                            for (int k=0; k<Fan_Ambient_Dimension; k++){
                                Rays_Matroid_Fan[k][j]=Rays_Matroid_Fan[k][j]-1;
                            }
                        }
                    }
                    
                    //Construct a map that maps the rays of the matroid fan to the rays of the input fan
                    //Extract some of the vectors of the matroid fan such that they are a basis of the hole space and write them into the matrix U
                    Set <int> Basis;
                    Basis= basis_cols(Rays_Matroid_Fan);
                    
                    Entire< Set<int> >::iterator it=entire(Basis);
                    Matrix<Integer> U(Basis.size(),Basis.size());
                    for (int j=0; j<Fan_Ambient_Dimension; j++) {
                        for (int i=0; i<Fan_Ambient_Dimension; i++){
                            U[i][j]=Rays_Matroid_Fan[i][*it];
                        }
                        *it++;
                    }
                    
                    U=inv(U);
                    
                    //Take the rays of the input fan and compute the corresponding primitive vactors
                    Matrix<Rational> Rays_Fan=f.give("RAYS");
                    Rays_Fan=T(Rays_Fan);

                    for (int j=0; j<Rays_Fan.cols(); j++) {
                        Integer LCM_Denominator=1;
                        for (int i=0; i<Rays_Fan.rows(); i++){
                            Integer Denominator=denominator(Rays_Fan[i][j]);
                            LCM_Denominator=lcm(LCM_Denominator,Denominator);
                        }
                        for (int i=0; i<Rays_Fan.rows(); i++){
                            Rays_Fan[i][j]=Rays_Fan[i][j]*LCM_Denominator;
                        }
                        Integer GCD_Numerator=1;
                        for (int i=0; i<Rays_Fan.rows(); i++){
                            GCD_Numerator=gcd(numerator(Rays_Fan[i][j]),GCD_Numerator);
                        }
                        for (int i=0; i<Rays_Fan.rows(); i++){
                            Rays_Fan[i][j]=Rays_Fan[i][j]/GCD_Numerator;
                        }
                    }

                    //Take the rays of the input fan that correspond to the rays of the matroid fan extracted above and write them into the matrix V
                    it=entire(Basis);
                    Matrix<Integer> V(Basis.size(),Basis.size());
                    for (int j=0; j<Fan_Ambient_Dimension; j++) {
                        for (int i=0; i<Fan_Ambient_Dimension; i++){
                            V[i][j]=Rays_Fan[i][*it];
                        }
                        *it++;
                    }
                    
                    Matrix<Integer> A=V*U; //searched map
                    
                    //Check if the map is a Z-isomorphisem and maps all rays of the matroid fan to the input fan
                    if ((det(A)==1 or det(A)==-1) and Rays_Fan==A*Rays_Matroid_Fan){
                        
                        
                        //Compute the Basis of the Matroid with the Flats found above
                        Vector < Set<int> > Bases;
                        int n=0;
                        
                        for (int i=0; i<Fan_Ambient_Dimension+1; i++) {
                            for (int j=i+1; j<Fan_Ambient_Dimension+1; j++) {
                                for (int k=j+1; k<Fan_Ambient_Dimension+1; k++) {
                                    bool possible=true;
                                    Set<int> Possible_Base;
                                    Possible_Base.insert(i);
                                    Possible_Base.insert(j);
                                    Possible_Base.insert(k);
                                    for (int l=0; l<N_Nodes+N_Edges; l++) {
                                        if (combinatorics[l][2]==2 and (Flats[l]*Possible_Base).size()==Possible_Base.size()) {
                                            possible=false;
                                        }
                                    }
                                    if (possible==true){
                                        Bases.resize(n+1);
                                        Bases[n]=Possible_Base;
                                        n=n+1;
                                    }
                                }
                            }
                        }
                        
                        perl::Object M("matroid::Matroid");
                        M.take("BASES") << Bases;
                        M.take("N_ELEMENTS") << Fan_Ambient_Dimension+1;

                        R.is_smooth=1;
                        R.matroid=M;
                        R.map=A;
                    }
                }
            }
            return R;
            
        }
        
        int check_matroid(Vector< Vector<int> > &combinatorics,  Vector < Set <int> > &Flats, int Fan_Ambient_Dimension){
            int s=1; //s=0 this is NOT a matroid, s=1 this is a matroid
            //pm::cout<<"Check second property"<<endl;
            //The second property of flats is checked
            for (int i=0; i<combinatorics.size(); i++){
                if (combinatorics[i][2]==2){
                    for (int j=0; j<combinatorics.size(); j++){
                        if (i!=j and combinatorics[j][2]==2) {
                            if((Flats[i]*Flats[j]).size()>1){
                                s=0;
                                goto label1;
                            }
                        }
                    }
                }
            }
            
            //The third property of flats is checked
            for (int n=0; n<=Fan_Ambient_Dimension; n++){
                //Search flats containing n
                Vector < Set<int> > Flats_containing_n; //and in reality setminus n
                for (int i=0; i<combinatorics.size(); i++){
                    if (combinatorics[i][2]==2 and Flats[i].contains(n)){
                        Flats_containing_n|=(Flats[i]-n);
                    }
                }
                
                //Check that they are disjoint
                for (int i=0; i<Flats_containing_n.size(); i++){
                    for(int j=i+1; j<Flats_containing_n.size(); j++){
                        if ((Flats_containing_n[i]*Flats_containing_n[j]).size()!=0) {
                            s=0;
                            goto label1;
                        }
                    }
                }
                
                //Check that the union is the complement of {n}
                Set<int> Union;
                for (int i=0; i<Flats_containing_n.size(); i++){
                    Union=Union+Flats_containing_n[i];
                }
                if (Union.size()!=Fan_Ambient_Dimension){
                    s=0;
                    goto label1;
                }
            }
            
        label1:
            return s;
        }
        
        result smooth_codim1(perl::Object f){
            result R;
            Matrix<Rational> Rays_Rational=f.give("RAYS");
            //Make the Rays integers
            Matrix<Integer> Rays(Rays_Rational.rows(), Rays_Rational.cols());
            for (int i=0; i<Rays.rows(); i++) {
                Integer LCM_Denominator=1;
                for (int j=0; j<Rays.cols(); j++){
                    Integer Denominator=denominator(Rays_Rational[i][j]);
                    LCM_Denominator=lcm(LCM_Denominator,Denominator);
                }
                for (int j=0; j<Rays.cols(); j++){
                    Rays_Rational[i][j]=Rays_Rational[i][j]*LCM_Denominator;
                }
                Integer GCD_Numerator=1;
                for (int j=0; j<Rays.cols(); j++){
                    GCD_Numerator=gcd(numerator(Rays_Rational[i][j]),GCD_Numerator);
                }
                for (int j=0; j<Rays.cols(); j++){
                    Rays[i][j]=Rays_Rational[i][j]/GCD_Numerator;
                }
            }

            Rays=T(Rays);
            Rays.resize(Rays.rows(), Rays.cols()-1);
            
            //Compute the Hermit normal form
            Matrix<Integer> tfmatrix;
            int kernelDimension;
            Matrix<Integer> HNF=lllHNF(Rays, tfmatrix, kernelDimension);
            
            //Compute the determinant of HNF
            Integer determinant=det(HNF);
            
            if (determinant==1 or determinant==-1) {
                IncidenceMatrix<> MaximalCones = f.give("MAXIMAL_CONES");
                int N_Maximal_Cones_Real=MaximalCones.rows();
                int Fan_Ambient_Dimension=f.give("FAN_AMBIENT_DIM");
                int N_Maximal_Cones_Expected=(Fan_Ambient_Dimension*(Fan_Ambient_Dimension+1))/2;
                if (N_Maximal_Cones_Real==N_Maximal_Cones_Expected) {
                    Array< Set<int> > Bases(Rays.cols()+1);
                    Set<int> E;
                    for (int i=0; i<Rays.cols()+1; i++) {
                        E=E+i;
                    }
                    for (int i=0; i<Rays.cols()+1; i++) {
                        Bases[i]=E-i;
                    }
                    perl::Object M("matroid::Matroid");
                    M.take("BASES") << Bases;
                    M.take("N_ELEMENTS") << Rays.cols()+1;

                    R.is_smooth=1;
                    R.matroid=M;
                    R.map=inv(tfmatrix);
                    return R;
                }
                else{
                    R.is_smooth=0;
                    return R;
                }
            }
            else{
                R.is_smooth=0;
                return R;
            }
        }
        
        UserFunction4perl("# @category Matroid fans"
                          "#Takes a weighted fan and returns if it is smooth "
			  "# (i.e. isomorphic to a Bergman fan B(M)/L for some matroid M) or not. "
                          "# The algorithm works for fans of dimension 0,1,2 and "
                          "# codimension 0,1! For other dimensions the algorithm "
                          "# could give an answer but it is not guaranteed. "
                          "# @param WeightedComplex a tropical fan F"
                          "# @return perl::ListReturn containing an integer s, "
                          "# a matroid M and a Matrix<Int> A. If s=1 then F is smooth, the "
                          "# corresponding matroid fan is Z-isomorphic to the matroid fan "
                          "# associated to M. The Z-isomorphism is given by A, i.e. A*F = B(M)/L."
			   "# If s=0, F is not smooth. If s=2 the algorithm is not able to determine "
                          "# if F is smooth or not. "
                          , is_smooth, "is_smooth(WeightedComplex)");
    }
}