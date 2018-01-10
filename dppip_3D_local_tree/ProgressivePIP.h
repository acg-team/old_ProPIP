#ifndef PROGRESSIVEPIP_H_
#define PROGRESSIVEPIP_H_

#include "Alphabet.h"
#include "PhyTree.h"
#include "ModelFactory.h"
#include "Sequence.h"
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "debug.h"
#include "Sequence.h"
#include <Eigen/Dense>
#include <random>

#include <set>
#include <chrono>

#include <Eigen/Sparse>

#define ERR_STATE 0
#define MATCH_STATE 1
#define GAP_X_STATE 2
#define GAP_Y_STATE 3
#define STOP_STATE 4

const char mytable[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, 2, -1, 1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, 0, 0, -1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, 2, -1, 1,
		-1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0,
		-1, -1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

const char mytableAA[256] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, 0, 20, 1, 2, 3, 4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14,
		15, 16, 20, 17, 18, 20, 19, 20, -1, -1, -1, -1, -1, -1, 0, 20, 1, 2, 3,
		4, 5, 6, 7, 20, 8, 9, 10, 11, 20, 12, 13, 14, 15, 16, 20, 17, 18, 20,
		19, 20, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1 };




template <class ALPHABET>
struct ProgressivePIPResult;

template <class ALPHABET>
struct ProgressivePIPResult{
	std::vector<std::pair<std::string,sequence_t<ALPHABET> > > MSA;
	std::vector< std::pair<std::string,std::string> > MSAs;
	sequence_t<ALPHABET> traceback_path;
	score_t score;
	Eigen::VectorXd Pc;
	std::map<std::string,Eigen::VectorXd> fv_map;
	double lk_gap;
	double pc0;
};

class ProgressivePIPException : public swps3_exception {
public:
	ProgressivePIPException(std::string msg) : swps3_exception(msg) {}
};

//=======================================================================================================
//DP-PIP
void dealloc_3D_matrix(long double ***mat,int depth,int height){

	for(int k=(depth-1); k>=0; k--){
	    for(int j=(height-1); j>=0; j--){
	        delete[] mat[k][j];
	    }
	    delete[] mat[k];
	}
	delete[] mat;

	mat=NULL;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
Eigen::VectorXd go_down(PhyTree &tree,std::string &s,int &idx,int ugglyflag){
	Eigen::VectorXd fv;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;

	if(tree.isLeaf()){

		fv=Eigen::VectorXd::Zero(ALPHABET::DIM+1);
		int ii;//=s[idx].value();
		if(ugglyflag==1){
			ii=mytable[(int)s[idx]];
		}else if(ugglyflag==2){
			ii=mytableAA[(int)s[idx]];
		}else{
			printf("go_down not implemented for codon model yet\n");
			exit(EXIT_FAILURE);
		}
		ii=ii<0?ALPHABET::DIM:ii;
		fv[ii]=1.0;
		idx++;

	}else{

		fvL=go_down<ALPHABET>(tree[0],s,idx,ugglyflag);
		fvR=go_down<ALPHABET>(tree[1],s,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

	}

	return fv;
}
//=======================================================================================================
//DP-PIP
double get_scale_factor(double m,double x,double y){

	double s;

	if(m==0.0){
		if(x==0.0){
			s=y;
		}else if(y==0.0){
			s=x;
		}else{
			s= (x>y) ? x : y;
		}
	}else if(x==0.0){
		if(m==0.0){
			s=y;
		}else if(y==0.0){
			s=m;
		}else{
			s= (m>y) ? m : y;
		}
	}else if(y==0.0){
		if(m==0.0){
			s=x;
		}else if(x==0.0){
			s=m;
		}else{
			s= (m>x) ? m : x;
		}
	}else{
		s= (m>x) ? m : x;
		s= (s>y) ? s : y;
	}

	return (-s);
}
//=======================================================================================================
//DP-PIP
void reset_corner(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int h,int w){
	int delta;

	if(up_corner_j>=w){
		delta=up_corner_j-w+1;
		up_corner_j-=delta;
		up_corner_i+=delta;
	}
	if(bot_corner_i>=h){
		delta=bot_corner_i-h+1;
		bot_corner_i-=delta;
		bot_corner_j+=delta;
	}

}
//=======================================================================================================
//DP-PIP
std::string create_col_MSA_gap(int len){

	std::string colMSA (len,'-');

	return colMSA;
}
//=======================================================================================================
//DP-PIP
int index_of_max(double m, double x, double y,double epsilon,std::default_random_engine &generator,std::uniform_real_distribution<double> &distribution){
	double random_number;

	int ERR=-1;

	if(not(std::isinf(m)) & not(std::isinf(x)) & (fabs(m-x)<epsilon)){
		x=m;
	}

	if(not(std::isinf(m)) & not(std::isinf(y)) & (fabs(m-y)<epsilon)){
		y=m;
	}

	if(not(std::isinf(x)) & not(std::isinf(y)) & (fabs(x-y)<epsilon)){
		y=x;
	}

	if(m>x){
		if(m>y){
			return int(MATCH_STATE);
		}else if (y>m){
			return int(GAP_Y_STATE);
		}else{
			if(abs(m-y)<epsilon){
				//m or y
				random_number  = distribution(generator);
				if(random_number < (1.0/2.0) ){
					return int(MATCH_STATE);
				}else{
					return int(GAP_Y_STATE);
				}
			}else{
				perror("ERROR in index_of_max_3\n");
				exit(EXIT_FAILURE);
			}
		}
	}else if (x>m){
		if(x>y){
			return int(GAP_X_STATE);
		}else if (y>x){
			return int(GAP_Y_STATE);
		}else{
			if(abs(x-y)<epsilon){
				//x or y
				random_number  = distribution(generator);
				if(random_number < (1.0/2.0) ){
					return int(GAP_X_STATE);
				}else{
					return int(GAP_Y_STATE);
				}
			}else{
				perror("ERROR in index_of_max_3\n");
				exit(EXIT_FAILURE);
			}
		}
	}else{

		double mx=x;
		if(mx>y){
			//m or x
			random_number  = distribution(generator);
			if(random_number < (1.0/2.0) ){
				return int(MATCH_STATE);
			}else{
				return int(GAP_X_STATE);
			}
		}else if (y>mx){
			return int(GAP_Y_STATE);
		}else{
			if(abs(mx-y)<epsilon){
				//m or x or y
				random_number  = distribution(generator);
				if(random_number < (1.0/3.0)){
					return int(MATCH_STATE);
				}else if(random_number < (2.0/3.0)){
					return int(GAP_X_STATE);
				}else{
					return int(GAP_Y_STATE);
				}
			}else{
				perror("ERROR in index_of_max_3\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	return ERR;
}
//=======================================================================================================
//DP-PIP
double max_of_three(double a, double b, double c,double epsilon){

	//------------------------------------
	if(fabs(a)<epsilon){
		a=-INFINITY;
	}
	if(fabs(b)<epsilon){
		b=-INFINITY;
	}
	if(fabs(c)<epsilon){
		c=-INFINITY;
	}
	//------------------------------------

	if(std::isinf(a) && std::isinf(b) && std::isinf(c)){
		perror("max_of_three_2: all inf\n");
		exit(EXIT_FAILURE);
	}

	if(a>b){
		if(a>c){
			return a;
		}
		return c;
	}else{
		if(b>c){
			return b;
		}
		return c;
	}

}
//=======================================================================================================
//DP-PIP
double compute_nu(double tau,double lambda,double mu){

	if(fabs(mu)<1e-8){
		error("ERROR in compute_nu: mu too small");
	}

	return lambda*(tau+1/mu);
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void check_sequences_tree(PhyTree &tree,const std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &sequences){

	if(tree.isLeaf()){

		const std::string name = tree.getName();

		typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>> ::const_iterator vect_iter;
		vect_iter it = std::find_if(sequences.begin(),sequences.end(),CompareFirst<ALPHABET>(name));


		if (it == sequences.end()){
			error("ERROR sequence name doesn't match any tree leaf");
			exit(EXIT_FAILURE);
		}

	}else{

		check_sequences_tree<ALPHABET>(tree[0],sequences);
		check_sequences_tree<ALPHABET>(tree[1],sequences);

	}

}
//=======================================================================================================
//DP-PIP
bool is_inside(int x0,int y0,int xf,int yf,int xt,int yt){

	if((xt<x0) || (yt>y0) || (xt>xf) || (yt<yf)){

		return false;
	}

	if( (y0-yt)>(xt-x0) ){
		return false;
	}

	return true;
}
//=======================================================================================================
//DP-PIP
void set_indeces_M(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

	if(level==0){
		up_corner_i=0;
		up_corner_j=0;
		bot_corner_i=0;
		bot_corner_j=0;
	}else{
		up_corner_i=1+level-std::min(w-1,level);
		up_corner_j=std::min(w-1,level);
		bot_corner_i=std::min(h-1,level);
		bot_corner_j=1+level-std::min(h-1,level);
	}

}
//=======================================================================================================
//DP-PIP
void set_indeces_X(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

	if(level==0){
		up_corner_i=0;
		up_corner_j=0;
		bot_corner_i=0;
		bot_corner_j=0;
	}else{
		up_corner_i=1+level-1-std::min(w-1,level-1);
		up_corner_j=std::min(w-1,level-1);
		bot_corner_i=std::min(h-1,level);
		bot_corner_j=level-std::min(h-1,level);
	}

}
//=======================================================================================================
//DP-PIP
void set_indeces_Y(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

	if(level==0){
		up_corner_i=0;
		up_corner_j=0;
		bot_corner_i=0;
		bot_corner_j=0;
	}else{
		up_corner_i=level-std::min(w-1,level);
		up_corner_j=std::min(w-1,level);
		bot_corner_i=std::min(h-1,level-1);
		bot_corner_j=1+level-1-std::min(h-1,level-1);
	}

}
//=======================================================================================================
//DP-PIP
void set_indeces_T(int &up_corner_i,int &up_corner_j,int &bot_corner_i, int &bot_corner_j,int level,int h,int w){

	int up_corner_i_x;
	int up_corner_i_y;

	int up_corner_j_x;
	int up_corner_j_y;

	int bot_corner_i_x;
	int bot_corner_i_y;

	int bot_corner_j_x;
	int bot_corner_j_y;

	set_indeces_X(up_corner_i_x,up_corner_j_x,bot_corner_i_x,bot_corner_j_x,level,h,w);

	set_indeces_Y(up_corner_i_y,up_corner_j_y,bot_corner_i_y,bot_corner_j_y,level,h,w);

	int delta_i,delta_j;

	delta_i=bot_corner_i_x-up_corner_i_y;
	delta_j=up_corner_j_y-bot_corner_j_x;

	if(delta_i>delta_j){
		up_corner_i=up_corner_i_y;
		up_corner_j=up_corner_j_y;
		bot_corner_i=up_corner_i_y+delta_i;
		bot_corner_j=up_corner_j_y-delta_i;
	}else{
		up_corner_i=bot_corner_i_x-delta_j;
		up_corner_j=bot_corner_j_x+delta_j;
		bot_corner_i=bot_corner_i_x;
		bot_corner_j=bot_corner_j_x;
	}

}
//=======================================================================================================
//DP-PIP
int get_indices_M(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

	int idx;

	set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

	if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

		int dx,sx;

		dx=nx-up_corner_i+1;

		sx=((dx+1)*dx/2)-1;

		idx=sx+(ny-up_corner_j);
	}else{
		idx=-999;
	}

	return idx;

}
//=======================================================================================================
//DP-PIP
int get_indices_X(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

	int idx;

	set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

	if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

		int dx,sx;

		dx=nx-up_corner_i+1;

		sx=((dx+1)*dx/2)-1;

		idx=sx+(ny-up_corner_j);
	}else{
		idx=-999;
	}

	return idx;

}
//=======================================================================================================
//DP-PIP
int get_indices_Y(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

	int idx;

	set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

	if(is_inside(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,nx,ny)){

		int dx,sx;

		dx=nx-up_corner_i+1;

		sx=((dx+1)*dx/2)-1;

		idx=sx+(ny-up_corner_j);
	}else{
		idx=-999;
	}

	return idx;

}
//=======================================================================================================
//DP-PIP
int get_indices_T(int nx,int ny,int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int m,int h,int w){

	set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

	reset_corner(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w);

	int idx;
	int dx,sx;

	dx=nx-up_corner_i+1;

	sx=((dx+1)*dx/2)-1;

	idx=sx+(ny-up_corner_j);

	return idx;

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void allgaps(PhyTree &tree,std::string &s,int &idx,bool &flag){

	if(tree.isLeaf()){
		char ch=s[idx];

		idx++;

		if(ch!='-'){
			flag=false;
		}

	}else{
		allgaps<ALPHABET>(tree[0],s,idx,flag);
		allgaps<ALPHABET>(tree[1],s,idx,flag);
	}

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double compute_lk_gap_down(PhyTree &tree,std::string &s,Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,int ugglyflag){

	double pr=0;
	double pL=0;
	double pR=0;
	int idx;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	double fv0;

	if(tree.isLeaf()){
		idx=0;
		fv=go_down<ALPHABET>(tree,s,idx,ugglyflag);
		fv0=fv.dot(pi);
		pr=tree.get_iota()-tree.get_iota()*tree.get_beta()+tree.get_iota()*tree.get_beta()*fv0;

		return pr;
	}else{
		idx=0;
		fv=go_down<ALPHABET>(tree,s,idx,ugglyflag);
		fv0=fv.dot(pi);
		pr=tree.get_iota()-tree.get_iota()*tree.get_beta()+tree.get_iota()*tree.get_beta()*fv0;

		bool flagL=true;
		bool flagR=true;
		idx=0;
		allgaps<ALPHABET>(tree[0],s,idx,flagL);
		int ixx=idx;
		allgaps<ALPHABET>(tree[1],s,idx,flagR);
		int len;

		std::string sL;//=stringFromSequence(s);
		len=ixx;
		sL=s.substr(0,len);
		pL=compute_lk_gap_down<ALPHABET>(tree[0],sL,pi,ugglyflag);

		std::string sR;//=stringFromSequence(s);
		sR=s.substr(ixx);
		pR=compute_lk_gap_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
	}

	return pr+pL+pR;
}
//=======================================================================================================
//DP-PIP
bool checkboundary(int up_corner_i,int up_corner_j,int bot_corner_i,int bot_corner_j,int h,int w){

	if( (up_corner_i  >=0) & (up_corner_i  <h) &\
	   (up_corner_j  >=0) & (up_corner_j  <w) &\
	   (bot_corner_i >=0) & (bot_corner_i <h) &\
	   (bot_corner_j >=0) & (bot_corner_j <w)){
		return true;
	}

	return false;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::string create_col_MSA(std::vector<std::pair<std::string,std::string>> &result,int index){
	std::string colMSA;

	for(unsigned int i=0;i<result.size();i++){
		colMSA.append(result.at(i).second,index,1);
	}

	return colMSA;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::vector<std::pair<std::string,std::string>> align_seq_left(	std::vector<std::pair<std::string,std::string>> &MSA_in,
																	sequence_t<ALPHABET> &traceback_path){

	unsigned int idx;

	std::vector<std::pair<std::string,std::string>> MSA_out;

	typedef typename std::vector<std::pair<std::string,std::string>>::iterator vect_iterator;
	for (vect_iterator iter = MSA_in.begin(); iter != MSA_in.end(); iter++){

		std::pair<std::string,std::string> seq = (*iter);

		std::string seq_name=seq.first;
		std::string seq_not_aligned=seq.second;
		std::string seq_aligned(traceback_path.size(),'-');

		idx=0;
		for(unsigned int j=0;j<traceback_path.size();j++){

			if(traceback_path.at(j)==ALPHABET::match){

				seq_aligned[j]=seq_not_aligned.at(idx);
				idx++;

			}else if(traceback_path.at(j)==ALPHABET::gapX){

					seq_aligned[j]=seq_not_aligned.at(idx);
					idx++;

			}else if(traceback_path.at(j)==ALPHABET::gapY){

					seq_aligned[j]='-';

			}else{
					error("ERROR in align_seq");
			}
		}

		MSA_out.push_back(std::make_pair(seq_name,seq_aligned));
	}

	return MSA_out;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::vector<std::pair<std::string,std::string>> align_seq_right(	std::vector<std::pair<std::string,std::string>> &result,
																																		sequence_t<ALPHABET> &traceback_path){
	unsigned int idx;

	std::vector<std::pair<std::string,std::string>> MSA_out;

	typedef typename std::vector<std::pair<std::string,std::string>>::iterator vect_iterator;
	for (vect_iterator iter = result.begin(); iter != result.end(); iter++){

		std::pair<std::string,std::string> seq = (*iter);

		std::string seq_name=seq.first;
		std::string seq_not_aligned=seq.second;
		std::string seq_aligned(traceback_path.size(),'-');

		idx=0;
		for(unsigned int j=0;j<traceback_path.size();j++){

			if(traceback_path.at(j)==ALPHABET::match){

				seq_aligned[j]=seq_not_aligned.at(idx);
				idx++;

			}else if(traceback_path.at(j)==ALPHABET::gapX){

					seq_aligned[j]='-';

			}else if(traceback_path.at(j)==ALPHABET::gapY){

				seq_aligned[j]=seq_not_aligned.at(idx);
					idx++;

			}else{
					error("ERROR in align_seq");
			}
		}

		MSA_out.push_back(std::make_pair(seq_name,seq_aligned));
	}

	return MSA_out;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
static std::vector<std::pair<std::string,std::string>> build_MSA(sequence_t<ALPHABET> traceback_path,std::vector<std::pair<std::string,std::string>> &MSA_L,std::vector<std::pair<std::string,std::string>> &MSA_R){
	std::vector<std::pair<std::string,std::string>> MSA;
	std::vector<std::pair<std::string,std::string>> MSA_L_out;
	std::vector<std::pair<std::string,std::string>> MSA_R_out;

	MSA_L_out=align_seq_left(MSA_L,traceback_path);
	MSA_R_out=align_seq_right(MSA_R,traceback_path);

	for (unsigned int i=0;i<MSA_L_out.size();i++){
		MSA.push_back(MSA_L_out.at(i));
	}

	for (unsigned int i=0;i<MSA_R_out.size();i++){
		MSA.push_back(MSA_R_out.at(i));
	}

	return MSA;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::vector<ProgressivePIPResult<ALPHABET>> extract_N_best(std::vector<ProgressivePIPResult<ALPHABET>> &vec,unsigned int n){

	std::vector<ProgressivePIPResult<ALPHABET>> best;

	for(unsigned int k=0;k<n;k++){

		double score=vec[0].score;
		unsigned int idx=0;
		for(unsigned int i=1;i<vec.size();i++){
			if(vec[i].score>score){
				idx=i;
				score=vec[i].score;
			}
		}

		best.push_back(vec[idx]);
		vec.erase(vec.begin() + idx);

	}

	return best;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double compute_lk_down(PhyTree &tree,std::string &s,Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,int ugglyflag){

	double pr;
	int idx;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	double fv0;

	if(tree.isLeaf()){

		idx=0;
		fv=go_down<ALPHABET>(tree,s,idx,ugglyflag);
		fv0=fv.dot(pi);
		pr=tree.get_iota()*tree.get_beta()*fv0;

		return pr;

	}else{

		idx=0;
		fv=go_down<ALPHABET>(tree,s,idx,ugglyflag);
		fv0=fv.dot(pi);
		pr=tree.get_iota()*tree.get_beta()*fv0;

		bool flagL=true;
		bool flagR=true;
		idx=0;
		allgaps<ALPHABET>(tree[0],s,idx,flagL);
		int ixx=idx;
		allgaps<ALPHABET>(tree[1],s,idx,flagR);

		int len;
		if(flagR){
			std::string sL;//=stringFromSequence(s);
			len=ixx;
			sL=s.substr(0,len);
			return pr + compute_lk_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
		}

		if(flagL){
			std::string sR;//=stringFromSequence(s);
			sR=s.substr(ixx);
			return pr + compute_lk_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
		}

	}

	return pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double compute_pr_gap_all_edges_s(	PhyTree &tree,
						std::string &sL,
						std::string &sR,
						Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,int ugglyflag){

	double fv0;
	double pr;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	int idx;

	idx=0;
	fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);

	idx=0;
	fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

	fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

	fv0=fv.dot(pi);

//	pr=tree.get_iota()*fv0;
	if(tree.getParent()==NULL){
		pr=(tree.get_iota()*fv0);
	}else{
		pr=(tree.get_iota() - tree.get_iota()*tree.get_beta() + tree.get_iota()*tree.get_beta()*fv0);
	}


	double pL,pR;
	pL=compute_lk_gap_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
	pR=compute_lk_gap_down<ALPHABET>(tree[1],sR,pi,ugglyflag);

	pr=pr+pL+pR;

	//*************************************************************************************************
	PhyTree *p_tree=&tree;

	while(p_tree->getParent()!=NULL){

		fv=(p_tree->get_Pr()*fv);

		fv0=fv.dot(pi);

		if(p_tree->getParent()->getParent()==NULL){
			pr+=(p_tree->getParent()->get_iota()*fv0);
		}else{
			pr+=(p_tree->getParent()->get_iota() - p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta() + p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);
		}

		p_tree=p_tree->getParent();
	}
	//*************************************************************************************************

	return pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double compute_pr_gap_local_tree_s(PhyTree &tree,
						std::string &sL,
						std::string &sR,
						Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,int ugglyflag){

	double fv0;
	double pr;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	int idx;

	idx=0;
	fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);

	idx=0;
	fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

	fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

	fv0=fv.dot(pi);

	pr=tree.get_iota()*fv0;

	double pL,pR;
	pL=compute_lk_gap_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
	pR=compute_lk_gap_down<ALPHABET>(tree[1],sR,pi,ugglyflag);

	pr=pr+pL+pR;

	return pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_M_all_edges_s_opt(double valM,
									double valX,
									double valY,
									double nu,
									PhyTree &tree,
									std::string &sL,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									int m,
									std::map<std::string,double> &lkM,
									int ugglyflag){


	double pr;
	double val;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(sL);
	s.append(sR);

	MapIterator it=lkM.find(s);
	if(it == lkM.end()){

		//------------------------------------------------------------------------------------------
		//DOWN
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		double fv0;
		int idx;

		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);
		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;
		//------------------------------------------------------------------------------------------
		//UP
		PhyTree *p_tree=&tree;

		while(p_tree->getParent()!=NULL){

			fv=(p_tree->get_Pr()*fv);

			fv0=fv.dot(pi);

			pr+=(p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);

			p_tree=p_tree->getParent();
		}
		//------------------------------------------------------------------------------------------

		pr=log(pr);

		lkM[s]=pr;

	}else{
		pr=it->second;
	}

	val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

	return val;
}
//=======================================================================================================
//DP-PIP
//@SB
template <class ALPHABET>
double computeLK_M_all_edges_s_opt_SB(PhyTree &tree,
									std::string &sL,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									std::map<std::string,double> &lkM,int ugglyflag){


	double log_pr;
	double pr;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(sL);
	s.append(sR);

	MapIterator it=lkM.find(s);
	if(it == lkM.end()){

		//------------------------------------------------------------------------------------------
		//DOWN
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		double fv0;
		int idx;

		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);
		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;
		//------------------------------------------------------------------------------------------
		//UP
		PhyTree *p_tree=&tree;

		while(p_tree->getParent()!=NULL){

			fv=(p_tree->get_Pr()*fv);

			fv0=fv.dot(pi);

			pr+=(p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);

			p_tree=p_tree->getParent();
		}
		//------------------------------------------------------------------------------------------

		log_pr=log(pr);

		lkM[s]=log_pr;

	}else{
		log_pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	return log_pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_M_local_tree_s_opt(double valM,
									double valX,
									double valY,
									double nu,
									PhyTree &tree,
									std::string &sL,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									int m,
									std::map<std::string,double> &lkM,
									int ugglyflag){


	double pr;
	double val;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	double fv0;
	int idx;

	s.append(sL);
	s.append(sR);

	MapIterator it=lkM.find(s);
	if(it == lkM.end()){

		//------------------------------------------------------------------------------------------


		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);
		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;

		//------------------------------------------------------------------------------------------

		pr=log(pr);

		lkM[s]=pr;

	}else{
		pr=it->second;
	}
	//------------------------------------------------------------------------------------------

#ifdef VERBOSE
	std::cout<<"prM("<<sL<<":"<<sR<<")"; printf(" %18.16lf\n",pr);
#endif

	val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

	return val;
}
//=======================================================================================================
//DP-PIP
//@SB
template <class ALPHABET>
double computeLK_M_local_tree_s_opt_SB(PhyTree &tree,
									std::string &sL,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									std::map<std::string,double> &lkM,int ugglyflag){


	double log_pr;
	double pr;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(sL);
	s.append(sR);

	MapIterator it=lkM.find(s);
	if(it == lkM.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		double fv0;
		int idx;

		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);
		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;
		//------------------------------------------------------------------------------------------

		log_pr=log(pr);

		lkM[s]= log_pr;

	}else{
		log_pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	return log_pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_X_all_edges_s_opt(double valM,
									double valX,
									double valY,
									double nu,
									PhyTree &tree,
									std::string &sL,
									std::string &col_gap_R,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									int m,
									std::map<std::string,double> &lkX,
									int ugglyflag){


	double pr;
	double val;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(sL);
	s.append(col_gap_R);

	MapIterator it=lkX.find(s);
	if(it == lkX.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		int idx;
		double fv0;

		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],col_gap_R,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;

		double pL;

		//------------------------------------------------------------------------------------------
		pL=compute_lk_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
		//------------------------------------------------------------------------------------------

		pr+=pL;
		//------------------------------------------------------------------------------------------


		//******************************************************************************************
		PhyTree *p_tree=&tree;

		while(p_tree->getParent()!=NULL){

			fv=(p_tree->get_Pr()*fv);

			fv0=fv.dot(pi);

			pr+=(p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);

			p_tree=p_tree->getParent();
		}
		//******************************************************************************************

		pr=log(pr);

		lkX[s]=pr;

	}else{
		pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

	return val;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_X_all_edges_s_opt_SB(PhyTree &tree,
									std::string &sL,
									std::string &col_gap_R,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									std::map<std::string,double> &lkX,int ugglyflag){


	double log_pr;
	double pr;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(sL);
	s.append(col_gap_R);

	MapIterator it=lkX.find(s);
	if(it == lkX.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		int idx;
		double fv0;

		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],col_gap_R,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;

		double pL;

		//------------------------------------------------------------------------------------------
		pL=compute_lk_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
		//------------------------------------------------------------------------------------------

		pr+=pL;
		//------------------------------------------------------------------------------------------


		//******************************************************************************************
		PhyTree *p_tree=&tree;

		while(p_tree->getParent()!=NULL){

			fv=(p_tree->get_Pr()*fv);

			fv0=fv.dot(pi);

			pr+=(p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);

			p_tree=p_tree->getParent();
		}
		//******************************************************************************************


		log_pr=log(pr);

		lkX[s]=log_pr;

	}else{
		log_pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	return log_pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_X_local_tree_s_opt(double valM,
									double valX,
									double valY,
									double nu,
									PhyTree &tree,
									std::string &sL,
									std::string &col_gap_R,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									int m,
									std::map<std::string,double> &lkX,
									int ugglyflag){


	double pr;
	double val;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	int idx;
	double fv0;

	s.append(sL);
	s.append(col_gap_R);

	MapIterator it=lkX.find(s);
	if(it == lkX.end()){

		//------------------------------------------------------------------------------------------


		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],col_gap_R,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;

		double pL;

		//------------------------------------------------------------------------------------------
		pL=compute_lk_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
		//------------------------------------------------------------------------------------------

		pr+=pL;
		//------------------------------------------------------------------------------------------

		pr=log(pr);

		lkX[s]=pr;

	}else{
		pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

	return val;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_X_local_tree_s_opt_SB(PhyTree &tree,
									std::string &sL,
									std::string &col_gap_R,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									std::map<std::string,double> &lkX,int ugglyflag){


	double log_pr;
	double pr;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(sL);
	s.append(col_gap_R);

	MapIterator it=lkX.find(s);
	if(it == lkX.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		int idx;
		double fv0;

		idx=0;
		fvL=go_down<ALPHABET>(tree[0],sL,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],col_gap_R,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;

#ifdef VERBOSE
		std::cout<<"prX("<<sL<<")"; printf(" %18.16lf\n",log_pr);
#endif

		double pL;

		//------------------------------------------------------------------------------------------
		pL=compute_lk_down<ALPHABET>(tree[0],sL,pi,ugglyflag);
		//------------------------------------------------------------------------------------------

		pr+=pL;
		//------------------------------------------------------------------------------------------

		log_pr=log(pr);

		lkX[s]=log_pr;

	}else{
		log_pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	return log_pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_Y_all_edges_s_opt(double valM,
									double valX,
									double valY,
									double nu,
									PhyTree &tree,
									std::string &col_gap_L,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									int m,
									std::map<std::string,double> &lkY,
									int ugglyflag){


	double pr;
	double val;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(col_gap_L);
	s.append(sR);

	MapIterator it=lkY.find(s);
	if(it == lkY.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		double fv0;
		int idx=0;
		fvL=go_down<ALPHABET>(tree[0],col_gap_L,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		//------------------------------------------------------------------------------------------

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;

		double pR;

		//------------------------------------------------------------------------------------------
		pR=compute_lk_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
		//------------------------------------------------------------------------------------------

		pr+=pR;

		//*******************************************************************************************
		PhyTree *p_tree=&tree;

		while(p_tree->getParent()!=NULL){

			fv=(p_tree->get_Pr()*fv);

			fv0=fv.dot(pi);

			pr+=(p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);

			p_tree=p_tree->getParent();
		}
		//*******************************************************************************************

		pr=log(pr);

		lkY[s]=pr;

	}else{
		pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);

	return val;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_Y_all_edges_s_opt_SB(PhyTree &tree,
									std::string &col_gap_L,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									std::map<std::string,double> &lkY,
									int ugglyflag){


	double log_pr;
	double pr;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(col_gap_L);
	s.append(sR);

	MapIterator it=lkY.find(s);
	if(it == lkY.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		double fv0;
		int idx=0;
		fvL=go_down<ALPHABET>(tree[0],col_gap_L,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		//------------------------------------------------------------------------------------------

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;


		double pR;

		//------------------------------------------------------------------------------------------
		pR=compute_lk_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
		//------------------------------------------------------------------------------------------

		pr+=pR;

		//*******************************************************************************************
		PhyTree *p_tree=&tree;

		while(p_tree->getParent()!=NULL){

			fv=(p_tree->get_Pr()*fv);

			fv0=fv.dot(pi);

			pr+=(p_tree->getParent()->get_iota()*p_tree->getParent()->get_beta()*fv0);

			p_tree=p_tree->getParent();
		}
		//*******************************************************************************************


		log_pr=log(pr);

		lkY[s]=log_pr;

	}else{
		log_pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	return log_pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_Y_local_tree_s_opt(double valM,
									double valX,
									double valY,
									double nu,
									PhyTree &tree,
									std::string &col_gap_L,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									int m,
									std::map<std::string,double> &lkY,
									int ugglyflag){


	double pr;
	double val;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;
	Eigen::VectorXd fvL;
	Eigen::VectorXd fvR;
	Eigen::VectorXd fv;
	double fv0;
	int idx=0;

	s.append(col_gap_L);
	s.append(sR);

	MapIterator it=lkY.find(s);
	if(it == lkY.end()){

		//------------------------------------------------------------------------------------------

		fvL=go_down<ALPHABET>(tree[0],col_gap_L,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		//------------------------------------------------------------------------------------------

		fv0=fv.dot(pi);

			pr=tree.get_iota()*tree.get_beta()*fv0;

			double pR;

			//------------------------------------------------------------------------------------------
			pR=compute_lk_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
			//------------------------------------------------------------------------------------------

			pr+=pR;

			pr=log(pr);

		lkY[s]=pr;

	}else{
		pr=it->second;
	}
	//------------------------------------------------------------------------------------------


	val=-log(double(m))+log(nu)+pr+max_of_three(valM,valX,valY,(double)DBL_EPSILON);


	return val;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
double computeLK_Y_local_tree_s_opt_SB(PhyTree &tree,
									std::string &col_gap_L,
									std::string &sR,
									Eigen::Matrix<score_t,ALPHABET::DIM+1,1> &pi,
									std::map<std::string,double> &lkY,
									int ugglyflag){


	double log_pr;
	double pr;
	//------------------------------------------------------------------------------------------
	typedef typename std::map<std::string,double>::iterator MapIterator;
	std::string s;

	s.append(col_gap_L);
	s.append(sR);

	MapIterator it=lkY.find(s);
	if(it == lkY.end()){

		//------------------------------------------------------------------------------------------
		Eigen::VectorXd fvL;
		Eigen::VectorXd fvR;
		Eigen::VectorXd fv;
		double fv0;
		int idx=0;
		fvL=go_down<ALPHABET>(tree[0],col_gap_L,idx,ugglyflag);

		idx=0;
		fvR=go_down<ALPHABET>(tree[1],sR,idx,ugglyflag);

		fv=(tree.get_left_child()->get_Pr()*fvL).cwiseProduct(tree.get_right_child()->get_Pr()*fvR);
		//------------------------------------------------------------------------------------------

		fv0=fv.dot(pi);

		pr=tree.get_iota()*tree.get_beta()*fv0;


			double pR;

			//------------------------------------------------------------------------------------------
			pR=compute_lk_down<ALPHABET>(tree[1],sR,pi,ugglyflag);
			//------------------------------------------------------------------------------------------

			pr+=pR;

		log_pr=log(pr);

		lkY[s]=log_pr;

	}else{
		log_pr=it->second;
	}
	//------------------------------------------------------------------------------------------

	return log_pr;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void forward_SB(	PhyTree &tree,
									long double*** LogM,
									long double*** LogX,
									long double*** LogY,
									int h,int w,int d,
									double nu,
									ProgressivePIPResult<ALPHABET> &result_L,
									ProgressivePIPResult<ALPHABET> &result_R,
									const ModelFactory<ALPHABET> *model_factory,
									std::map<std::string,double> &lkM,
									std::map<std::string,double> &lkX,
									std::map<std::string,double> &lkY,
									bool local,
									int ugglyflag){

	int up_corner_i;
	int up_corner_j;
	int bot_corner_i;
	int bot_corner_j;

	int lw;

	std::string sLs;
	std::string sRs;
	std::string col_gap_Ls;
	std::string col_gap_Rs;

	double log_pr,Log_pr;
	double valM;
	double valX;
	double valY;

	int coordSeq_1;
	int coordSeq_2;
	int coordTriangle_this_i;
	int coordTriangle_this_j;
	int coordTriangle_prev_i;
	int coordTriangle_prev_j;

	double pc0;

	double log1;
	double log2;
	double log3;
	double e1,e2,e3,etot;

	//*******************************************************************************************************//
	col_gap_Ls=create_col_MSA_gap(result_L.MSAs.size());
	col_gap_Rs=create_col_MSA_gap(result_R.MSAs.size());

	Eigen::Matrix<score_t,ALPHABET::DIM+1,1> pi = model_factory->getPiPIP();

	if(local){
		pc0=compute_pr_gap_local_tree_s<ALPHABET>(tree,
										col_gap_Ls,
										col_gap_Rs,
										pi,
										ugglyflag);
	}else{
		pc0=compute_pr_gap_all_edges_s<ALPHABET>(tree,
										col_gap_Ls,
										col_gap_Rs,
										pi,
										ugglyflag);
	}

	LogM[0][0][0]=nu*(pc0-1.0);

	double scale_factor;

	for(int m=1;m<d;m++){

		//***************************************************************************************
		//***************************************************************************************
		set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){

				coordTriangle_this_i=i;
				coordSeq_1=coordTriangle_this_i-1;
				coordTriangle_prev_i=coordTriangle_this_i-1;
				sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1);

				for(int j=0;j<=lw;j++){

					coordTriangle_this_j=up_corner_j-j;
					coordSeq_2=coordTriangle_this_j-1;
					coordTriangle_prev_j=coordTriangle_this_j-1;
					sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2);

					valM=LogM[m-1][coordTriangle_prev_i][coordTriangle_prev_j];
					valX=LogX[m-1][coordTriangle_prev_i][coordTriangle_prev_j];
					valY=LogY[m-1][coordTriangle_prev_i][coordTriangle_prev_j];

					if(local){
						log_pr=computeLK_M_local_tree_s_opt_SB<ALPHABET>(tree,sLs,sRs,pi,lkM,ugglyflag);
					}else{
						log_pr=computeLK_M_all_edges_s_opt_SB<ALPHABET>(tree,sLs,sRs,pi,lkM,ugglyflag);
					}

					if( (valM==0.0) && (valX==0.0) && (valY==0.0) ){
						log1=0.0;
						Log_pr=0.0;
					}else{

						scale_factor=get_scale_factor(valM,valX,valY);


						e1 = (valM==0.0) ? 0.0 : exp(valM+scale_factor);
						e2 = (valX==0.0) ? 0.0 : exp(valX+scale_factor);
						e3 = (valY==0.0) ? 0.0 : exp(valY+scale_factor);

						etot = e1 + e2 + e3;

						log1=log(etot) - scale_factor;

						if(log_pr>=0.0){
							printf("log_pr > 0: %lf\n",log_pr);
							printf("scale factor %lf\n",scale_factor);
							exit(1);
						}
						if(log1>=0.0){
							printf("log1 > 0: %lf\n",log1);
							printf("scale factor %lf\n",scale_factor);
							printf("M %lf X %lf Y %lf\n",valM,valX,valY);
							exit(1);
						}

					}


					Log_pr = log_pr + log1;

					if ( std::isinf(Log_pr) || std::isnan(Log_pr) ){
						printf("ERROR isinf or isnan (M)\n");
						printf("log_pr %lf ;log1 %lf\n",log_pr,log1);
						printf("scale factor %lf\n",scale_factor);
						printf("M %lf X %lf Y %lf\n",valM,valX,valY);
						printf("e1 %lf e2 %lf e3 %lf\n",e1,e2,e3);
						printf("etot %lf\n",etot);
						exit(EXIT_FAILURE);
					}

					LogM[m][coordTriangle_this_i][coordTriangle_this_j]=Log_pr;
				}
				lw++;
			}

		}
		//***************************************************************************************
		//***************************************************************************************
		set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){

				coordTriangle_this_i=i;
				coordTriangle_prev_i=coordTriangle_this_i-1;
				coordSeq_1=coordTriangle_this_i-1;
				sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1);

				for(int j=0;j<=lw;j++){

					coordTriangle_this_j=up_corner_j-j;
					coordTriangle_prev_j=coordTriangle_this_j;

					valM=LogM[m-1][coordTriangle_prev_i][coordTriangle_prev_j];
					valX=LogX[m-1][coordTriangle_prev_i][coordTriangle_prev_j];
					valY=LogY[m-1][coordTriangle_prev_i][coordTriangle_prev_j];


					if(local){
						log_pr=computeLK_X_local_tree_s_opt_SB<ALPHABET>(tree,sLs,col_gap_Rs,pi,lkX,ugglyflag);
					}else{
						log_pr=computeLK_X_all_edges_s_opt_SB<ALPHABET>(tree,sLs,col_gap_Rs,pi,lkX,ugglyflag);
					}

					if( (valM==0.0) && (valX==0.0) && (valY==0.0) ){
						log2=0.0;
						Log_pr=0.0;
					}else{

						scale_factor=get_scale_factor(valM,valX,valY);


						e1 = (valM==0.0) ? 0.0 : exp(valM+scale_factor);
						e2 = (valX==0.0) ? 0.0 : exp(valX+scale_factor);
						e3 = (valY==0.0) ? 0.0 : exp(valY+scale_factor);

						etot = e1 + e2 + e3;

						log2=log(etot) - scale_factor;


						if(log_pr>=0.0){
							printf("log_pr > 0:  %lf\n",log_pr);
							printf("scale factor %lf\n",scale_factor);
							exit(1);
						}
						if(log2>=0.0){
							printf("log2 > 0: %lf\n",log2);
							printf("scale factor %lf\n",scale_factor);
							printf("M %lf X %lf Y %lf\n",valM,valX,valY);
							exit(1);
						}

					}



					Log_pr = log_pr + log2;

					if ( std::isinf(Log_pr) || std::isnan(Log_pr) ){
						exit(EXIT_FAILURE);
					}

					LogX[m][coordTriangle_this_i][coordTriangle_this_j]=Log_pr;
				}
				lw++;
			}

		}
		//***************************************************************************************
		//***************************************************************************************
		set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){
				coordTriangle_this_i=i;
				coordTriangle_prev_i=coordTriangle_this_i;
				for(int j=0;j<=lw;j++){

					coordTriangle_this_j=up_corner_j-j;
					coordTriangle_prev_j=coordTriangle_this_j-1;
					coordSeq_2=coordTriangle_this_j-1;
					sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2);

					valM=LogM[m-1][coordTriangle_prev_i][coordTriangle_prev_j];
					valX=LogX[m-1][coordTriangle_prev_i][coordTriangle_prev_j];
					valY=LogY[m-1][coordTriangle_prev_i][coordTriangle_prev_j];


					if(local){
						log_pr=computeLK_Y_local_tree_s_opt_SB<ALPHABET>(tree,col_gap_Ls,sRs,pi,lkY,ugglyflag);
					}else{
						log_pr=computeLK_Y_all_edges_s_opt_SB<ALPHABET>(tree,col_gap_Ls,sRs,pi,lkY,ugglyflag);
					}

					if( (valM==0.0) && (valX==0.0) && (valY==0.0) ){
						log3=0.0;
						Log_pr=0.0;
					}else{

						scale_factor=get_scale_factor(valM,valX,valY);


						e1 = (valM==0.0) ? 0.0 : exp(valM+scale_factor);
						e2 = (valX==0.0) ? 0.0 : exp(valX+scale_factor);
						e3 = (valY==0.0) ? 0.0 : exp(valY+scale_factor);

						etot = e1 + e2 + e3;

						log3=log(etot) - scale_factor;

						if(log_pr>=0.0){
							printf("log_pr > 0: %lf\n",log_pr);
							printf("scale factor %lf\n",scale_factor);
							printf("M %lf X %lf Y %lf\n",valM,valX,valY);
							exit(1);
						}
						if(log3>=0.0){
							printf("log3 > 0: %lf\n",log3);
							printf("scale factor %lf\n",scale_factor);
							printf("M %lf X %lf Y %lf\n",valM,valX,valY);
							exit(1);
						}

					}


					Log_pr = log_pr + log3;

					if ( std::isinf(Log_pr) || std::isnan(Log_pr) ){
						exit(EXIT_FAILURE);
					}

					LogY[m][coordTriangle_this_i][coordTriangle_this_j]=Log_pr;
				}
				lw++;
			}

		}
		//***************************************************************************************
		//***************************************************************************************
	}

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
ProgressivePIPResult<ALPHABET> backward_SB(PhyTree &tree,
																								long double*** LogM,
																								long double*** LogX,
																								long double*** LogY,
																								int d,int h,int w,
																								ProgressivePIPResult<ALPHABET> &result_L,
																								ProgressivePIPResult<ALPHABET> &result_R,
																								const ModelFactory<ALPHABET> *model_factory,
																								double nu,
																								std::map<std::string,double> &lkM,
																								std::map<std::string,double> &lkX,
																								std::map<std::string,double> &lkY,
																								bool local,
																								std::uniform_real_distribution<double> &distribution,
																								std::default_random_engine &generator,
																								int ugglyflag){

	double m;
	double pm,px,py;
	double Z;
	double lk;
	int STATUS;
	double log_pr;
	std::string sLs;
	std::string sRs;
	std::string col_gap_Ls;
	std::string col_gap_Rs;
	int coordSeq_1,coordSeq_2;
	int prev_i,prev_j;
	double number;
	int depth;

	//----------------------------------------------------------------------------------------------------------------------------//
	std::cout<<"left: "<<std::endl;
	for(unsigned int i=0;i<result_L.MSAs.size();i++){
		std::cout<<"key:"<<result_L.MSAs.at(i).first <<" val: "<<result_L.MSAs.at(i).second;std::cout<<"\n";
	}

	std::cout<<"right: "<<std::endl;
	for(unsigned int i=0;i<result_R.MSAs.size();i++){
		std::cout<<"key:"<<result_R.MSAs.at(i).first <<" val: "<<result_R.MSAs.at(i).second;std::cout<<"\n";
	}
	//----------------------------------------------------------------------------------------------------------------------------//

	col_gap_Ls=create_col_MSA_gap(result_L.MSAs.size());
	col_gap_Rs=create_col_MSA_gap(result_R.MSAs.size());

	Eigen::Matrix<score_t,ALPHABET::DIM+1,1> pi = model_factory->getPiPIP();

	std::string traceback_path;

	//*********************************************************************//
	number = distribution(generator);
	STATUS=index_of_max(	LogM[d-2][h-1][w-1],
														LogX[d-1][h-1][w-1],
														LogY[d-1][h-1][w-1],
														(double)DBL_EPSILON,
														generator,
														distribution);
	//*********************************************************************//
	coordSeq_1 = h;
	coordSeq_2 = w;
	prev_i = h-1;
	prev_j = w-1;
	//*********************************************************************//
	switch (STATUS) {
		case MATCH_STATE:
			coordSeq_1--;
			coordSeq_2--;
			prev_i--;
			prev_j--;
			traceback_path.append("1");

			depth=d-2;

			break;
		case GAP_X_STATE:
			coordSeq_1--;
			prev_i--;
			traceback_path.append("2");

			depth=d-1;

			break;
		case GAP_Y_STATE:
			coordSeq_2--;
			prev_j--;
			traceback_path.append("3");

			depth=d-1;

			break;
		default:
			error("ERROR in backward\n");
			break;
	}
	//*********************************************************************//

	m = 0;

	lk=0.0;

	bool cond=true;

	while(cond){

		m = m + 1;

		//**********************************************
		switch (STATUS) {
			case MATCH_STATE:
				//-----------------------------------------------------------------------------------------
				sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1-1);
				sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2-1);

				if(local){
					log_pr=computeLK_M_local_tree_s_opt_SB<ALPHABET>(tree,sLs,sRs,pi,lkM,ugglyflag);
				}else{
					log_pr=computeLK_M_all_edges_s_opt_SB<ALPHABET>(tree,sLs,sRs,pi,lkM,ugglyflag);
				}
				//-----------------------------------------------------------------------------------------
				break;
			case GAP_X_STATE:
				//-----------------------------------------------------------------------------------------
				sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1-1);

				if(local){
					log_pr=computeLK_X_local_tree_s_opt_SB<ALPHABET>(tree,sLs,col_gap_Rs,pi,lkX,ugglyflag);
				}else{
					log_pr=computeLK_X_all_edges_s_opt_SB<ALPHABET>(tree,sLs,col_gap_Rs,pi,lkX,ugglyflag);
				}
				//-----------------------------------------------------------------------------------------
				break;
			case GAP_Y_STATE:
				//-----------------------------------------------------------------------------------------
				sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2-1);

				if(local){
					log_pr=computeLK_Y_local_tree_s_opt_SB<ALPHABET>(tree,col_gap_Ls,sRs,pi,lkY,ugglyflag);
				}else{
					log_pr=computeLK_Y_all_edges_s_opt_SB<ALPHABET>(tree,col_gap_Ls,sRs,pi,lkY,ugglyflag);
				}
				//-----------------------------------------------------------------------------------------
				break;
			default:
				error("ERROR in backward\n");
				break;
		}
		//**********************************************

		for(int dd=depth-1;dd>=0;dd--){
				pm = LogM[dd][prev_i][prev_j];// + log_pr;
				if(pm<0.0){
					break;
				}
		}
		for(int dd=depth-1;dd>=0;dd--){
				px = LogX[dd][prev_i][prev_j];// + log_pr;
				if(px<0.0){
					break;
				}
		}
		for(int dd=depth-1;dd>=0;dd--){
				py = LogY[dd][prev_i][prev_j];// + log_pr;
				if(py<0.0){
					break;
				}
		}


		//**********************************************
		Z = fabs(pm) + fabs(px) + fabs(py);
		pm = fabs(pm)/Z;
		px = fabs(px)/ Z;
		py= fabs(py)/Z;
		//**********************************************
		//**********************************************
		double Temperature = 0.1;
		pm = (pm==0.0) ? 0.0 : exp(-(1-pm)/Temperature);
		px = (px==0.0) ? 0.0 : exp(-(1-px)/Temperature);
		py = (py==0.0) ? 0.0 : exp(-(1-py)/Temperature);
		Z = pm + px + py;
		pm = pm/Z;
		px = px/ Z;
		py= py/Z;
		//**********************************************

		number = distribution(generator);

		if (number < pm) {
			STATUS = (int)MATCH_STATE;
			traceback_path.append("1");
			prev_i--;
			prev_j--;
			coordSeq_1--;
			coordSeq_2--;
		} else if (number < (pm + px)) {
			STATUS = (int)GAP_X_STATE;
			traceback_path.append("2");
			prev_i--;
			coordSeq_1--;
		} else {
			STATUS = (int)GAP_Y_STATE;
			traceback_path.append("3");
			prev_j--;
			coordSeq_2--;
		}
		//**********************************************

		depth--;

		lk += (log(nu) - log(m) + log_pr);

		if(prev_i==0 && prev_j==0){
			cond=false;
		}

		//**********************************************

	};

	m=m+1;

	//**********************************************
	switch (STATUS) {
		case MATCH_STATE:
			//-----------------------------------------------------------------------------------------
			sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1-1);
			sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2-1);

			if(local){
				log_pr=computeLK_M_local_tree_s_opt_SB<ALPHABET>(tree,sLs,sRs,pi,lkM,ugglyflag);
			}else{
				log_pr=computeLK_M_all_edges_s_opt_SB<ALPHABET>(tree,sLs,sRs,pi,lkM,ugglyflag);
			}
			//-----------------------------------------------------------------------------------------
			break;
		case GAP_X_STATE:
			//-----------------------------------------------------------------------------------------
			sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1-1);

			if(local){
				log_pr=computeLK_X_local_tree_s_opt_SB<ALPHABET>(tree,sLs,col_gap_Rs,pi,lkX,ugglyflag);
			}else{
				log_pr=computeLK_X_all_edges_s_opt_SB<ALPHABET>(tree,sLs,col_gap_Rs,pi,lkX,ugglyflag);
			}
			//-----------------------------------------------------------------------------------------
			break;
		case GAP_Y_STATE:
			//-----------------------------------------------------------------------------------------
			sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2-1);

			if(local){
				log_pr=computeLK_Y_local_tree_s_opt_SB<ALPHABET>(tree,col_gap_Ls,sRs,pi,lkY,ugglyflag);
			}else{
				log_pr=computeLK_Y_all_edges_s_opt_SB<ALPHABET>(tree,col_gap_Ls,sRs,pi,lkY,ugglyflag);
			}
			//-----------------------------------------------------------------------------------------
			break;
		default:
			error("ERROR in backward\n");
			break;
	}
	//**********************************************
	if (STATUS==1) {
		lk += (log(nu) - log(m) + log_pr);
	} else if (STATUS==2) {
		lk += (log(nu) - log(m) + log_pr);
	} else {
		lk += (log(nu) - log(m) + log_pr);
	}

	//flip left-right
	traceback_path=std::string(traceback_path.rbegin(),traceback_path.rend());


	ProgressivePIPResult<ALPHABET> result;

	result.traceback_path=sequenceFromStringPIP<ALPHABET>(traceback_path);
	result.score=lk;
	result.MSAs=build_MSA(result.traceback_path,result_L.MSAs,result_R.MSAs);

	//----------------------------------------------------------------------------------------------------------------------------//
	std::cout<<"trace: "<<traceback_path<<std::endl;
	for(unsigned int i=0;i<result.MSAs.size();i++){
		std::cout<<"key:"<<result.MSAs.at(i).first <<" val: "<<result.MSAs.at(i).second;std::cout<<"\n";
	}
	//----------------------------------------------------------------------------------------------------------------------------//



	return result;
}
//=======================================================================================================
//DP-PIP
void fill_scores(double lk,double &max_lk,double &prev_max_lk,int &level_max_lk,int &last_d,int m,bool &flag_exit,double *scores,int &counter,bool CENTER,int num_subopt,int &index0){

	if (lk>max_lk){
		prev_max_lk=max_lk;
		max_lk=lk;
		level_max_lk=m;

		if (CENTER){
			if(std::isinf(prev_max_lk)){
				scores[0]=max_lk;
				counter=1;
				index0=m;
			}else{
				if(num_subopt>2){
					scores[0]=prev_max_lk;
					scores[1]=max_lk;
					counter=2;
					index0=m-1;
				}else{
					scores[0]=max_lk;
					counter=1;
					index0=m;
				}
			}
		}else{
			scores[0]=max_lk;
			counter=1;
			index0=m;
		}

	}else{

		if (counter>=num_subopt){
			flag_exit=true;
			last_d=m;
			return;
		}

		scores[counter]=lk;
		counter=counter+1;

	}

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
static std::vector<ProgressivePIPResult<ALPHABET>> compute_DP3D_PIP(ProgressivePIPResult<ALPHABET> &result_L,ProgressivePIPResult<ALPHABET> &result_R,PhyTree &tree,const ModelFactory<ALPHABET> *model_factory,double tau,double nu,double gamma_rate,int num_subopt,bool local,int ugglyflag,std::ostream* &out_score){

	//@gamma_distribution
	double lambda_gamma=model_factory->lambda*gamma_rate;
	double mu_gamma=model_factory->mu*gamma_rate;

	if(local){
		tau=tree.computeLength();
		nu=lambda_gamma*(tau+1/mu_gamma);
		tree.set_tau(tau);
		tree.set_nu(nu);
		//@gamma_distribution
		tree.set_iota_local(tau,mu_gamma);
		tree.set_beta_local(tau,mu_gamma);
	}

	int up_corner_i;
	int up_corner_j;
	int bot_corner_i;
	int bot_corner_j;
	int lw;
	int h,w;

	h=get_length_seq_s<ALPHABET>(result_L.MSAs)+1;
	w=get_length_seq_s<ALPHABET>(result_R.MSAs)+1;
	int d=(h-1)+(w-1)+1;

	Eigen::Matrix<score_t,ALPHABET::DIM+1,1> pi = model_factory->getPiPIP();

	double pc0;
	//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
	std::string sLs;
	std::string sRs;
	std::string col_gap_Ls;
	std::string col_gap_Rs;
	//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
	col_gap_Ls=create_col_MSA_gap(result_L.MSAs.size());
	col_gap_Rs=create_col_MSA_gap(result_R.MSAs.size());
	//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
//	std::cout<<"random generator ON\n";
	unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);
	//.-----.------.------.//
//	std::cout<<"random generator OFF\n";
//	unsigned seed = 0;
//	std::default_random_engine generator(seed);
//	std::uniform_real_distribution<double> distribution(0.0,1.0);
	//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

	double epsilon=DBL_EPSILON;

	//***************************************************************************************
	//***************************************************************************************
	if(local){
		pc0=compute_pr_gap_local_tree_s<ALPHABET>(tree,
										col_gap_Ls,
										col_gap_Rs,
										pi,
										ugglyflag);
	}else{
		pc0=compute_pr_gap_all_edges_s<ALPHABET>(tree,
										col_gap_Ls,
										col_gap_Rs,
										pi,
										ugglyflag);
	}

	//***************************************************************************************
	//***************************************************************************************

	double** LogM = new double*[2];
	double** LogX = new double*[2];
	double** LogY = new double*[2];

	int** TR = new int*[d];

	LogM[0] = new double[int((w*(h+1))/2)];
	LogX[0] = new double[int((w*(h+1))/2)];
	LogY[0] = new double[int((w*(h+1))/2)];
	LogM[1] = new double[int((w*(h+1))/2)];
	LogX[1] = new double[int((w*(h+1))/2)];
	LogY[1] = new double[int((w*(h+1))/2)];

	LogM[0][0]=nu*(pc0-1.0);
	LogX[0][0]=nu*(pc0-1.0);
	LogY[0][0]=nu*(pc0-1.0);

	TR[0] = new int[1]();
	TR[0][0]=STOP_STATE;

	double max_of_3;
	double max_lk=-INFINITY;
	double prev_max_lk=-INFINITY;
	int level_max_lk=INT_MIN;
	double val;
	int m_binary_this;
	int m_binary_prev;

	double valM;
	double valX;
	double valY;

	int idx;

	int coordSeq_1;
	int coordSeq_2;
	int coordTriangle_this_i;
	int coordTriangle_this_j;
	int coordTriangle_prev_i;
	int coordTriangle_prev_j;

	int counter;

	double *scores = new double[num_subopt];

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	bool CENTER = true;
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	int start_depth;
	int depth;


	bool flag_exit=false;
	int last_d=d-1;
	int size_tr,tr_up_i,tr_up_j,tr_down_i,tr_down_j;
	std::map<std::string,double> lkM;
	std::map<std::string,double> lkX;
	std::map<std::string,double> lkY;

	counter=0;
	for(int m=1;m<d;m++){

		if(flag_exit){
			break;
		}

		m_binary_this=m%2;
		m_binary_prev=(m+1)%2;
		//***************************************************************************************
		//***************************************************************************************
		set_indeces_M(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){

				coordTriangle_this_i=i;
				coordSeq_1=coordTriangle_this_i-1;
				coordTriangle_prev_i=coordTriangle_this_i-1;
				sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1);

				for(int j=0;j<=lw;j++){

					coordTriangle_this_j=up_corner_j-j;
					coordSeq_2=coordTriangle_this_j-1;
					coordTriangle_prev_j=coordTriangle_this_j-1;
					sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2);

					idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valM=LogM[m_binary_prev][idx];
					}else{
						valM=-INFINITY;
					}

					idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valX=LogX[m_binary_prev][idx];
					}else{
						valX=-INFINITY;
					}

					idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valY=LogY[m_binary_prev][idx];
					}else{
						valY=-INFINITY;
					}

					if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
						exit(EXIT_FAILURE);
					}

					if(local){
						val=computeLK_M_local_tree_s_opt<ALPHABET>(	valM,
																	valX,
																	valY,
																	nu,
																	tree,
																	sLs,sRs,
																	pi,
																	m,
																	lkM,
																	ugglyflag);
					}else{
						val=computeLK_M_all_edges_s_opt<ALPHABET>(	valM,
																valX,
																valY,
																nu,
																tree,
																sLs,sRs,
																pi,
																m,
																lkM,
																ugglyflag);
					}

					if(std::isinf(val)){
						exit(EXIT_FAILURE);
					}

					if(std::isnan(val)){
						exit(EXIT_FAILURE);
					}

					idx=get_indices_M(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

					LogM[m_binary_this][idx]=val;
				}
				lw++;
			}

		}
		//***************************************************************************************
		//***************************************************************************************
		set_indeces_X(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
		tr_down_i=bot_corner_i;
		tr_down_j=bot_corner_j;
		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){

				coordTriangle_this_i=i;
				coordTriangle_prev_i=coordTriangle_this_i-1;
				coordSeq_1=coordTriangle_this_i-1;
				sLs=create_col_MSA<ALPHABET>(result_L.MSAs,coordSeq_1);

				for(int j=0;j<=lw;j++){

					coordTriangle_this_j=up_corner_j-j;
					coordTriangle_prev_j=coordTriangle_this_j;

					idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valM=LogM[m_binary_prev][idx];
					}else{
						valM=-INFINITY;
					}

					idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valX=LogX[m_binary_prev][idx];
					}else{
						valX=-INFINITY;
					}

					idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valY=LogY[m_binary_prev][idx];
					}else{
						valY=-INFINITY;
					}

					if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
						exit(EXIT_FAILURE);
					}

					if(local){
						val=computeLK_X_local_tree_s_opt<ALPHABET>(valM,
																valX,
																valY,
																nu,
																tree,
																sLs,col_gap_Rs,
																pi,
																m,
																lkX,
																ugglyflag);
					}else{
						val=computeLK_X_all_edges_s_opt<ALPHABET>(valM,
															valX,
															valY,
															nu,
															tree,
															sLs,col_gap_Rs,
															pi,
															m,
															lkX,
															ugglyflag);
					}

					if(std::isinf(val)){
						exit(EXIT_FAILURE);
					}

					if(std::isnan(val)){
						exit(EXIT_FAILURE);
					}

					idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

					LogX[m_binary_this][idx]=val;
				}
				lw++;
			}

		}
		//***************************************************************************************
		//***************************************************************************************
		set_indeces_Y(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
		tr_up_i=up_corner_i;
		tr_up_j=up_corner_j;
		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){
				coordTriangle_this_i=i;
				coordTriangle_prev_i=coordTriangle_this_i;
				for(int j=0;j<=lw;j++){

					coordTriangle_this_j=up_corner_j-j;
					coordTriangle_prev_j=coordTriangle_this_j-1;
					coordSeq_2=coordTriangle_this_j-1;
					sRs=create_col_MSA<ALPHABET>(result_R.MSAs,coordSeq_2);

					idx=get_indices_M(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valM=LogM[m_binary_prev][idx];
					}else{
						valM=-INFINITY;
					}

					idx=get_indices_X(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valX=LogX[m_binary_prev][idx];
					}else{
						valX=-INFINITY;
					}

					idx=get_indices_Y(coordTriangle_prev_i,coordTriangle_prev_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m-1,h,w);
					if(idx>=0){
						valY=LogY[m_binary_prev][idx];
					}else{
						valY=-INFINITY;
					}

					if(std::isinf(valM) && std::isinf(valX) && std::isinf(valY)){
						exit(EXIT_FAILURE);
					}

					if(local){
						val=computeLK_Y_local_tree_s_opt<ALPHABET>(valM,
																valX,
																valY,
																nu,
																tree,
																col_gap_Ls,sRs,
																pi,
																m,
																lkY,
																ugglyflag);
					}else{
						val=computeLK_Y_all_edges_s_opt<ALPHABET>(valM,
															valX,
															valY,
															nu,
															tree,
															col_gap_Ls,sRs,
															pi,
															m,
															lkY,
															ugglyflag);
					}

					if(std::isinf(val)){
						exit(EXIT_FAILURE);
					}


					if(std::isnan(val)){
						exit(EXIT_FAILURE);
					}

					idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

					LogY[m_binary_this][idx]=val;
				}
				lw++;
			}

		}
		//***************************************************************************************
		//***************************************************************************************
		size_tr=int(ceil((tr_down_i-tr_up_i+1)*(tr_up_j-tr_down_j+1+1)/2));
		TR[m] = new int[size_tr](); /*TODO: optimize size TR*/
		memset(TR[m],0,size_tr*sizeof(TR[m][0]));
		set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

		if(checkboundary(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,h,w)){

			lw=0;
			for(int i=up_corner_i;i<=bot_corner_i;i++){
				coordTriangle_this_i=i;
				for(int j=0;j<=lw;j++){
					coordTriangle_this_j=up_corner_j-j;

					double mval;
					double xval;
					double yval;

					idx=get_indices_M(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
					if(idx>=0){
						mval=LogM[m_binary_this][idx];
					}else{
						mval=-INFINITY;
					}

					idx=get_indices_X(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
					if(idx>=0){
						xval=LogX[m_binary_this][idx];
					}else{
						xval=-INFINITY;
					}

					idx=get_indices_Y(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);
					if(idx>=0){
						yval=LogY[m_binary_this][idx];
					}else{
						yval=-INFINITY;
					}

					mval=fabs(mval)<epsilon?-INFINITY:mval;
					xval=fabs(xval)<epsilon?-INFINITY:xval;
					yval=fabs(yval)<epsilon?-INFINITY:yval;

					//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
					int ttrr;

					ttrr=index_of_max(mval,xval,yval,epsilon,generator,distribution);

					idx=get_indices_T(coordTriangle_this_i,coordTriangle_this_j,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,m,h,w);

					if(TR[m][idx]!=0){
						exit(EXIT_FAILURE);
					}

					TR[m][idx]=ttrr;

					//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

					//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//
					if( (coordTriangle_this_i==(h-1)) & (coordTriangle_this_j==(w-1)) ){

						max_of_3=max_of_three(mval,xval,yval,epsilon);

						fill_scores(max_of_3,max_lk,prev_max_lk,level_max_lk,last_d,m,flag_exit,scores,counter,CENTER,num_subopt,start_depth);

					}
					//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//



				}
				lw++;
			}
		}
	}

	std::vector<ProgressivePIPResult<ALPHABET>> result_v;
	for(int k=0;k<num_subopt;k++){
		depth=start_depth+k;

		if(depth>=d){
			break;
		}

		ProgressivePIPResult<ALPHABET> result;
		result.score=scores[k];


		//.................................
		out_score->precision(15);
		*out_score<<result.score<<"\n";
		//.................................


		sequence_t<ALPHABET> traceback_path(depth,ALPHABET::unknow);
		int id1=h-1;
		int id2=w-1;
		for(int lev=depth;lev>0;lev--){
			set_indeces_T(up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,lev,h,w);
			idx=get_indices_T(id1,id2,up_corner_i,up_corner_j,bot_corner_i,bot_corner_j,lev,h,w);
			switch(TR[lev][idx]){
				case MATCH_STATE:
					id1=id1-1;
					id2=id2-1;
					traceback_path[lev-1]=ALPHABET::match;
					break;
				case GAP_X_STATE:
					id1=id1-1;
					traceback_path[lev-1]=ALPHABET::gapX;
					break;
				case GAP_Y_STATE:
					id2=id2-1;
					traceback_path[lev-1]=ALPHABET::gapY;
					break;
				default:
					error("ERROR in alignment_reconstruction !!!");
					exit(EXIT_FAILURE);
			}
		}
		result.traceback_path=traceback_path;
		result.MSAs=build_MSA<ALPHABET>(traceback_path,result_L.MSAs,result_R.MSAs);
		result_v.push_back(result);

	}
	delete scores;
	//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

	free(LogM[1]);
	free(LogM[0]);
	free(LogM);

	free(LogX[1]);
	free(LogX[0]);
	free(LogX);

	free(LogY[1]);
	free(LogY[0]);
	free(LogY);

	for(int i=last_d;i>=0;i--){
		free(TR[i]);
	}
	free(TR);
	//.-----.------.------.-------.------.-------.-------.--------.-------.--------.//

	return result_v;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
static std::vector<ProgressivePIPResult<ALPHABET>> compute_DP3D_PIP_cross(	std::vector<ProgressivePIPResult<ALPHABET>> &result_L,
																																							std::vector<ProgressivePIPResult<ALPHABET>> &result_R,
																																							PhyTree &tree,
																																							const ModelFactory<ALPHABET> *model_factory,
																																							double tau,
																																							double nu,
																																							double gamma_rate,
																																							int num_subopt,
																																							bool local_tree,
																																							bool stoch_backtracking_flag,
																																							int ugglyflag,
																																							std::ostream* &out_score){

	std::vector<ProgressivePIPResult<ALPHABET>> result;
	std::vector<ProgressivePIPResult<ALPHABET>> subopt;
	std::vector<ProgressivePIPResult<ALPHABET>> orig_result_L=result_L;
	std::vector<ProgressivePIPResult<ALPHABET>> orig_result_R=result_R;

	for(unsigned int i=0;i<result_L.size();i++){
		for(unsigned int j=0;j<result_R.size();j++){

//			if(stoch_backtracking_flag){
//				//@SB
//				result=compute_DP3D_PIP_SB<ALPHABET>(orig_result_L.at(i),orig_result_R.at(j),tree,model_factory,tau,nu,gamma_rate,num_subopt,local_tree,ugglyflag);
//			}else{
				result=compute_DP3D_PIP<ALPHABET>(orig_result_L.at(i),orig_result_R.at(j),tree,model_factory,tau,nu,gamma_rate,num_subopt,local_tree,ugglyflag,out_score);
//			}

			orig_result_L=result_L;
			orig_result_R=result_R;

			for(unsigned int z=0;z<result.size();z++){
				subopt.push_back(result.at(z));
			}
		}
	}


	unsigned int N;
	N=result.size() < (unsigned int)num_subopt ? result.size() : (unsigned int)num_subopt;
	std::vector<ProgressivePIPResult<ALPHABET>> result_v=extract_N_best(subopt,N);

	return result_v;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::vector<ProgressivePIPResult<ALPHABET>> compute_DP3D_PIP_tree_cross(PhyTree &tree,
																																						const ModelFactory<ALPHABET> *model_factory,
																																						std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &sequences,
																																						double tau,
																																						double nu,
																																						double gamma_rate,
																																						int num_subopt,
																																						bool local_tree,
																																						bool stoch_backtracking_flag,
																																						int ugglyflag,
																																						std::ostream* &out_score) throw (ProgressivePIPException){

	std::vector<ProgressivePIPResult<ALPHABET>> result;

	if(tree.isLeaf()){

		result.resize(1);
		add_sequence_to_vector_string(result[0].MSAs,tree.getName(),sequences);

	}else{

		if(tree.n_children() != 2){
			error("only bifurcating trees allowed");
		}

		std::vector<ProgressivePIPResult<ALPHABET>>  result_L = compute_DP3D_PIP_tree_cross<ALPHABET>(tree[0],model_factory,sequences,tau,nu,gamma_rate,num_subopt,local_tree,stoch_backtracking_flag,ugglyflag,out_score);
		std::vector<ProgressivePIPResult<ALPHABET>>  result_R = compute_DP3D_PIP_tree_cross<ALPHABET>(tree[1],model_factory,sequences,tau,nu,gamma_rate,num_subopt,local_tree,stoch_backtracking_flag,ugglyflag,out_score);

		result=compute_DP3D_PIP_cross<ALPHABET>(result_L,result_R,tree,model_factory,tau,nu,gamma_rate,num_subopt,local_tree,stoch_backtracking_flag,ugglyflag,out_score);

	}

	return result;
}
//=======================================================================================================

#endif /* PROGRESSIVEPIP_H_ */
