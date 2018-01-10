/*
 * ModelFactoryJCPIP.h
 *
 *  Created on: May 24, 2016
 *      Author: max
 */

#ifndef MODELFACTORYJCPIP_H_
#define MODELFACTORYJCPIP_H_

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryJCPIP : public ModelFactory<DNA> {
public:
	ModelFactoryJCPIP(double mu,double lambda,double k);
	virtual ~ModelFactoryJCPIP() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MODELFACTORYJCPIP_H_ */
