/*
 * ModelFactoryJCPIP.cpp
 *
 *  Created on: May 24, 2016
 *      Author: max
 */


#include <iostream>
#include <sstream>
#include <cmath>
#include "debug.h"
#include <Eigen/Dense>
//==============================
//DP-PIP
#include "ModelFactoryJCPIP.h"
//==============================

ModelFactoryJCPIP::ModelFactoryJCPIP(double mu,double lambda,double k){

	double alpha = k/4.0;

	static double data[] = {
		-3*alpha,
		alpha,
		alpha,
		alpha,
		alpha,
		-3*alpha,
		alpha,
		alpha,
		alpha,
		alpha,
		-3*alpha,
		alpha,
		alpha,
		alpha,
		alpha,
		-3*alpha
      };


   this->Q = Eigen::Map<Model<DNA>::Subst>(data);
   Eigen::EigenSolver<Model<DNA>::Subst> solver2(this->Q.transpose());
   Model<DNA>::Freqs sigma2 = solver2.eigenvalues().real();
   Model<DNA>::Subst V2 = solver2.eigenvectors().real();

   Model<DNA>::Freqs::Index izero;
   sigma2.maxCoeff(&izero);

   assert(std::abs(sigma2(izero)) < 1e-8 && "Invalid Q-Matrix");

   this->freqs = V2.col(izero)/V2.col(izero).sum();

   // normalize rate
   this->Q.diagonal().setZero();
   this->Q.diagonal() = -this->Q.rowwise().sum().eval();
   this->Q /= -(this->freqs.transpose() * this->Q.diagonal())(0,0);

   //==================================================================

   Eigen::EigenSolver<Model<DNA>::Subst> solver_PIP(this->Q.transpose());
   Model<DNA>::Freqs sigma_PIP = solver_PIP.eigenvalues().real();
   Model<DNA>::Subst V_PIP = solver_PIP.eigenvectors().real();

   Model<DNA>::Freqs::Index izero_PIP;

   sigma_PIP.maxCoeff(&izero_PIP);

   this->mu=mu;
   this->lambda=lambda;

   for(int i=0;i<this->Q.rows();i++){
	   for(int j=0;j<this->Q.cols();j++){
		   this->Q_PIP(i,j)=this->Q(i,j);
		   this->Q_PIP(DNA::DIM,j)=0.0;
	   }
	   this->Q_PIP(i,DNA::DIM)=mu;
	   this->freqs_PIP(i)=this->freqs(i);
   }
   this->freqs_PIP(DNA::DIM)=0.0;
   this->Q_PIP(DNA::DIM,DNA::DIM)=0.0;

   //normalize rate
   this->Q_PIP.diagonal().setZero();
   this->Q_PIP.diagonal() = -this->Q_PIP.rowwise().sum().eval();

   Eigen::EigenSolver<Model<DNA>::Subst_PIP> solver(this->Q_PIP);
   this->sigma_PIP = solver.eigenvalues().real();
   this->V_PIP = solver.eigenvectors().real();
   this->Vi_PIP = this->V_PIP.inverse();

}







