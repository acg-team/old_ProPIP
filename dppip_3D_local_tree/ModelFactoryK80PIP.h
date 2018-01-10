/*
 * ModelFactoryK80PIP.h
 *
 *  Created on: May 13, 2016
 *      Author: max
 */

#ifndef MODELFACTORYK80PIP_H_
#define MODELFACTORYK80PIP_H_

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryK80PIP : public ModelFactory<DNA> {
public:
	ModelFactoryK80PIP(double mu,double lambda,double k);
	virtual ~ModelFactoryK80PIP() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MODELFACTORYK80PIP_H_ */
