/*
 * ModelFactoryHKY85PIP.h
 *
 *  Created on: May 13, 2016
 *      Author: max
 */

#ifndef MODELFACTORYHKY85PIP_H_
#define MODELFACTORYHKY85PIP_H_

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryHKY85PIP : public ModelFactory<DNA> {
public:
	ModelFactoryHKY85PIP(double mu,double lambda);
	virtual ~ModelFactoryHKY85PIP() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MODELFACTORYHKY85PIP_H_ */
