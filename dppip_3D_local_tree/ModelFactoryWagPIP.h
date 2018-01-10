#ifndef MODELFACTORYWAGPIP_H_
#define MODELFACTORYWAGPIP_H_

#include "ModelFactory.h"
#include "main.h"

class ModelFactoryWagPIP : public ModelFactory<AA> {
public:
	ModelFactoryWagPIP(double mu,double lambda);
	virtual ~ModelFactoryWagPIP() {};
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif /* MODELFACTORYWAGPIP_H_ */
