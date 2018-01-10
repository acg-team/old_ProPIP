#ifndef PHYTREE_H_
#define PHYTREE_H_

#include <vector>
#include <string>
#include <iostream>
#include <assert.h>
#include <sstream>
#include <cmath>
#include "main.h"

//=====================
//DP-PIP
#include <Eigen/Core>
#include "Alphabet.h"
#include <Eigen/src/Core/IO.h>
#include <Eigen/Dense>
#include "ModelFactory.h"
#include "Model.h"
#include "Sequence.h"
#include <unsupported/Eigen/MatrixFunctions>
//=====================

class PhyTree {

private:
	std::string name;

	double branch_length;
	double branch_support;

	PhyTree *parent;
	std::vector<PhyTree*> children;

	//================================================================
	//DP-PIP
    double iota;
    double beta;
    Eigen::MatrixXd Pr;

    double tau;
    double nu;


    //================================================================

	void print_prefix(std::string prefix) const {

		//====================================================================================
		//DP-PIP
		std::cout << prefix << "(" << this->branch_length << ") " << this->name << " " <<" iota: " << this->iota <<" beta: "<<this->beta<<std::endl;
		//====================================================================================

		for(std::vector<PhyTree*>::const_iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->print_prefix(prefix+"  ");
		}
	}

	void fixDistancesR() {
		if(cmdlineopts.mldist_flag || cmdlineopts.mldist_gap_flag) {
			if(std::isnan(this->branch_length)) this->branch_length = cmdlineopts.max_dist;
			this->branch_length = std::min(std::max(cmdlineopts.min_dist,this->branch_length),cmdlineopts.max_dist);
		} else {
			if(std::isnan(this->branch_length)) this->branch_length = cmdlineopts.max_pdist;
			this->branch_length = std::min(std::max(cmdlineopts.min_pdist,this->branch_length),cmdlineopts.max_pdist);
		}
		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->fixDistancesR();
		}
	}

	std::string formatNewickR() const {
		if(this->isLeaf()) {
			return this->getName();
		} else {
			std::stringstream newick;
			newick << "(";
			std::vector<PhyTree*>::const_iterator i=this->children.begin();
			newick << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
			for(++i; i < this->children.end(); ++i) {
				newick << "," << (*i)->formatNewickR() << ":" << (*i)->getBranchLength();
			}
			newick << ")";
			return newick.str();
		}
	}

public:

	//================================================================

	PhyTree(std::string name="") {
		this->parent = NULL;
		this->branch_length = 0;
		this->branch_support = 1;
		this->name = name;

		//==============================
		//DP-PIP
		this->iota = 0;
		this->beta = 0;
		this->Pr.resize(0,0);
		this->tau=0;
		this->nu=0;
		//==============================

	}

	~PhyTree() {
		assert(this->parent == NULL);
		for(std::vector<PhyTree*>::reverse_iterator i=this->children.rbegin(); i < this->children.rend(); ++i) {
			PhyTree *child = *i;
			child->parent = NULL;
			child->branch_length = 0;

			delete child;
		}
	}

	PhyTree* copy() {
		PhyTree* out = new PhyTree();
		out->branch_length = this->branch_length;
		out->branch_support = this->branch_support;
		out->name = this->name;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			out->addChild((*i)->copy(),(*i)->branch_length,(*i)->branch_support);
		}
		return out;
	}

	void addChild(PhyTree *child, double branch_length = 0, double branch_support = 1) {
		assert(child != this);
		assert(child->parent == NULL);

		this->children.push_back(child);
		child->parent = this;
		child->branch_length = branch_length;
		child->branch_support = branch_support;
	}

	index_t indexOf() {
		PhyTree *parent = this->parent;
		assert(parent != NULL);

		for(index_t i=0; i < parent->children.size(); ++i) {
			if(parent->children[i] == this) {
				return i;
			}
		}

		assert(false);
		return -1;
	}

	void pluck() {
		assert(this->parent != NULL);

		index_t index = this->indexOf();
		std::vector<PhyTree*>::iterator iter = this->parent->children.begin()+index;

		this->parent->children.erase(iter);
		this->parent = NULL;

		this->branch_length = 0;
		this->branch_support = 1;
	}

	PhyTree* pluckChild(index_t index) {
		std::vector<PhyTree*>::iterator iter = this->children.begin()+index;
		PhyTree *child = *iter;

		this->children.erase(iter);
		child->parent = NULL;
		child->branch_length = 0;
		this->branch_support = 1;

		return child;
	}

	void deleteChild(index_t index) {
		PhyTree *child = this->pluckChild(index);
		delete child;
	}

	void fixDistances() {
		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->fixDistancesR();
		}
	}

	index_t countLeaves() {
		if(this->isLeaf()) {
			return 1;
		} else {
			index_t leaves = 0;
			for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
				leaves += (*i)->countLeaves();
			}
			return leaves;
		}
	}

	double computeLength() {
		if(this->isLeaf()) {
			return 0;
		} else {
			double length = 0;
			for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
				length += (*i)->branch_length + (*i)->computeLength();
			}
			return length;
		}
	}

	std::string getName() const {
		return this->name;
	}

	PhyTree *getParent() {
		return this->parent;
	}

	const PhyTree *getParent() const {
		return this->parent;
	}

	double getBranchLength() const {
		return this->branch_length;
	}

	double getBranchSupport() const {
		return this->branch_support;
	}

	PhyTree& operator[](int i) {
		return *this->children[i];
	}

	const PhyTree& operator[](int i) const {
		return *this->children[i];
	}

	std::vector<PhyTree*>::iterator firstChild() {
		return this->children.begin();
	}

	std::vector<PhyTree*>::iterator lastChild() {
		return this->children.end();
	}

	std::vector<PhyTree*>::const_iterator firstChild() const {
		return this->children.begin();
	}

	std::vector<PhyTree*>::const_iterator lastChild() const {
		return this->children.end();
	}

	index_t n_children() const {
		return this->children.size();
	}

	bool isLeaf() const {
		return this->children.empty();
	}

	void print() const {
		this->print_prefix("");
	}

	std::string formatNewick() const {
		return this->formatNewickR() + ";";
	}

	//=======================================================================================================
	//DP-PIP
	std::vector<PhyTree*> get_children(){

		return this->children;
	}
	//=======================================================================================================
	//DP-PIP
	PhyTree* get_left_child(){

		if(this->n_children()>0){
			return this->children[0];
		}else{
			return NULL;
		}

	}
	//=======================================================================================================
	//DP-PIP
	PhyTree* get_right_child(){

		if(this->n_children()>0){
			return this->children[1];
		}else{
			return NULL;
		}

	}
	//=======================================================================================================
	//DP-PIP
	double get_iota() const {
		return this->iota;
	}
	//=======================================================================================================
	//DP-PIP
	double get_beta() const {
		return this->beta;
	}
	//=======================================================================================================
	//DP-PIP
	static bool empty(Eigen::VectorXd &v){

		if(v.rows() * v.cols() == 0){
			return true;
		}else{
			return false;
		}
	}
	//=======================================================================================================
	//DP-PIP
	void setName(std::string name){
		this->name=name;
	}
	//=======================================================================================================
	//DP-PIP
	void set_missing_node_name(std::string s){

		static int num=0;

		if(this->isLeaf()){

		}else{

			if(this->name.empty()){
				std::string name(s);
				name.append(std::to_string(num++));
				this->setName(name);
			}

		}

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_missing_node_name(s);
		}

	}
	//=======================================================================================================
	//DP-PIP
	const Eigen::MatrixXd& get_Pr(){

		return this->Pr;

	}
	//=======================================================================================================
	//DP-PIP
	void set_iota(double tau,double mu){

		if(fabs(mu)<1e-8){
			error("ERROR in set_iota: mu too small");
		}

		double T=tau+1/mu;

		if(fabs(T)<1e-8){
			error("ERROR in set_iota: T too small");
		}

		if(this->parent==NULL){
			this->iota=(1/mu)/T;
		}else{
			this->iota=this->branch_length/T;
		}

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_iota(tau,mu);
		}

	}
	//=======================================================================================================
	//DP-PIP
	double get_nu_local() const {
		return this->nu;
	}
	//=======================================================================================================
	//DP-PIP
	double get_tau_local() const {
		return this->tau;
	}
	//=======================================================================================================
	//DP-PIP
	void compute_nu_local(double tau,double lambda,double mu){
		this->nu=lambda*(tau+1/mu);
	}
	//=======================================================================================================
	//DP-PIP
	void compute_tau_local(){

		if(this->isLeaf()) {
			this->tau=this->branch_length;
		} else {
			this->get_left_child()->compute_tau_local();
			this->get_right_child()->compute_tau_local();

			this->tau=this->get_left_child()->tau+this->get_right_child()->tau;
		}

	}
	//=======================================================================================================
	//DP-PIP
	void set_iota_local_down(double tau,double mu){

		this->iota=this->branch_length/(tau+1/mu);

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_iota_local_down(tau,mu);
		}

	}
	//=======================================================================================================
	//DP-PIP
	void set_iota_local(double tau,double mu){

		this->iota=(1/mu)/(tau+1/mu);

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_iota_local_down(tau,mu);
		}


	}
	//=======================================================================================================
	//DP-PIP
	void set_beta(double tau,double mu){

		if(fabs(mu)<1e-8){
			error("ERROR in set_beta: mu too small");
		}

		if(this->parent==NULL){
			this->beta=1;
		}else{

			if(fabs(this->branch_length)<1e-8){
				error("ERROR in set_beta: branch_length too small");
			}

			this->beta=(1-exp(-mu*this->branch_length))/(mu*this->branch_length);
		}

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_beta(tau,mu);
		}

	}
	//=======================================================================================================
	//DP-PIP
	void set_beta_down(double tau,double mu){

		if(fabs(mu)<1e-8){
			error("ERROR in set_beta: mu too small");
		}

		if(fabs(this->branch_length)<1e-8){
			error("ERROR in set_beta: branch_length too small");
		}

		this->beta=(1-exp(-mu*this->branch_length))/(mu*this->branch_length);

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_beta_down(tau,mu);
		}

	}
	//=======================================================================================================
	//DP-PIP
	void set_beta_local(double tau,double mu){

		this->beta=1.0;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_beta_down(tau,mu);
		}

	}
	//=======================================================================================================
	//DP-PIP
	void set_tau(double tau){

		this->tau=tau;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_tau(tau);
		}

	}
	//=======================================================================================================
	//DP-PIP
	void set_nu(double nu){

		this->nu=nu;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->set_nu(nu);
		}

	}
	//=======================================================================================================
	//DP-PIP
	void print_local_var(){

		std::cout<<"----------------------\n";
		std::cout<<"Name: "<<this->name<<"\n";
		std::cout<<"tau: "<<this->tau<<"\n";
		std::cout<<"nu: "<<this->nu<<"\n";
		std::cout<<"iota: "<<this->iota<<"\n";
		std::cout<<"beta: "<<this->beta<<"\n";

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->print_local_var();
		}
	}
	//=======================================================================================================
	//DP-PIP
	void print_br(){

		std::cout<<"Name: "<<this->name<<" bl: "<<this->branch_length<<std::endl;

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->print_br();
		}
	}
	//=======================================================================================================
	//DP-PIP
	void printOnlyName(){

		std::cout<<"node name: "<<this->name<<"\n";
		if(this->n_children()==2){
			std::cout<<"      L child: "<<this->get_left_child()->getName()<<"\n";
			std::cout<<"      R child: "<<this->get_right_child()->getName()<<"\n";
		}


		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->printOnlyName();
		}


	}
	//=======================================================================================================
	//DP-PIP
	template <class ALPHABET>
	void initPrPIP(const ModelFactory<ALPHABET> *model_factory,double gamma_rate){

		//@gamma_distribution
		Model<ALPHABET> model = model_factory->getPIPModel(this->branch_length*gamma_rate);

		if(this->Pr.rows()*this->Pr.cols()!=0){
			this->Pr.resize(0,0);
		}
		if(this->parent!=NULL){
			this->Pr=model.P_PIP;
		}

		for(std::vector<PhyTree*>::iterator i=this->children.begin(); i < this->children.end(); ++i) {
			(*i)->initPrPIP<ALPHABET>(model_factory,gamma_rate);
		}

	}

};


PhyTree* midpointRoot(PhyTree *root);

std::vector<std::string> get_tree_order_ancestral(const PhyTree *tree);
std::vector<std::string> get_tree_order(const PhyTree *tree);

#endif /* PHYTREE_H_ */
