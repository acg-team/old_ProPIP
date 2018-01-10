/*
 * ModelFactoryHKY85PIP.cpp
 *
 *  Created on: May 13, 2016
 *      Author: max
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include "debug.h"
#include <Eigen/Dense>
//==============================
//DP-PIP
#include "ModelFactoryHKY85PIP.h"
//==============================

ModelFactoryHKY85PIP::ModelFactoryHKY85PIP(double mu,double lambda){
   static double data[] = {
		-1.0000,
		0.1667,
		0.6667,
		0.1667,
		0.1667,
		-1.0000,
		0.1667,
		0.6667,
		0.6667,
		0.1667,
		-1.0000,
		0.1667,
		0.1667,
		0.6667,
		0.1667,
		-1.0000
      };

   //==================================================================
   //DP-PIP


   error("NOT implemented... data matrix is random\n");
   exit(EXIT_FAILURE);


   //Compute QLD
   this->Q = Eigen::Map<Model<DNA>::Subst>(data);

   Eigen::EigenSolver<Model<DNA>::Subst> solver_PIP(this->Q.transpose());
   Model<DNA>::Freqs sigma_PIP = solver_PIP.eigenvalues().real();
   Model<DNA>::Subst V_PIP = solver_PIP.eigenvectors().real();

   Model<DNA>::Freqs::Index izero_PIP;

   sigma_PIP.maxCoeff(&izero_PIP);
   assert(std::abs(sigma_PIP(izero_PIP)) < 1e-8 && "Invalid Q-Matrix");

   score_t col_sum=V_PIP.col(izero_PIP).sum();
   assert(std::abs(col_sum) > 1e-8 && "Division by 0!");
   this->freqs = V_PIP.col(izero_PIP)/col_sum;

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

   score_t norm_fact=-(this->freqs_PIP.transpose() * this->Q_PIP.diagonal())(0,0);
   assert(std::abs(norm_fact) > 1e-8 && "Division by 0!");
   this->Q_PIP /= norm_fact;

   Eigen::EigenSolver<Model<DNA>::Subst_PIP> solver(this->Q_PIP);
   this->sigma_PIP = solver.eigenvalues().real();
   this->V_PIP = solver.eigenvectors().real();
   this->Vi_PIP = this->V_PIP.inverse();

}




