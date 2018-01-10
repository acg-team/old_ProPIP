#include <fstream>

#include "main.h"
#include "ModelFactory.h"
#include "ModelFactoryDarwin.h"
#include "ModelFactoryWag.h"
#include "ModelFactoryEcm.h"
#include "ModelFactoryCustom.h"
#include "ModelFactoryPlusF.h"
//==============================
//DP-PIP
#include "ModelFactoryWagPIP.h"
#include "ModelFactoryHKY85PIP.h"
#include "ModelFactoryK80PIP.h"
#include "ModelFactoryJCPIP.h"
//==============================

//============================================================================================================================
//#ifdef WITH_AA
template<>
ModelFactory<AA> *ModelFactory<AA>::getDefault(const std::map<std::string,sequence_t<AA> > &seqs)
{
	ModelFactory<AA> *model_factory = NULL;

	if(cmdlineopts.cmodel_file != "") {
		std::ifstream qmat_file(cmdlineopts.cmodel_file.c_str());
		model_factory = new ModelFactoryCustom<AA>(qmat_file);
	} else if(cmdlineopts.darwin_flag) {
		model_factory = new ModelFactoryDarwin();
	} else {
		model_factory = new ModelFactoryWag();
	}

	if(cmdlineopts.aafreqs_flag) {
		std::vector<sequence_t<AA> > seq_vector;
		for (std::map<std::string, sequence_t<AA> >::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
			seq_vector.push_back((*it).second);
		}
		ModelFactoryPlusF<AA> *plusF = new ModelFactoryPlusF<AA>(model_factory);
		plusF->estimateFreqs(seq_vector);
		model_factory = plusF;
	}

	return model_factory;
}
//#endif

//#ifdef WITH_AA
//============================================================================================================================
//DP-PIP
template <>
ModelFactory<AA> *ModelFactory<AA>::getPIP(double mu,double lambda)
{
	ModelFactory<AA> *model_factory = NULL;

	if(cmdlineopts.cmodel_file != "") {
		std::cerr<<"ERROR in getPIP: Model not implemented\n";
		exit(EXIT_FAILURE);
	} else if(cmdlineopts.darwin_flag) {
		std::cerr<<"ERROR in getPIP: Model not implemented\n";
		exit(EXIT_FAILURE);
	} else {
		model_factory = new ModelFactoryWagPIP(mu,lambda);
	}

	if(cmdlineopts.aafreqs_flag) {
		std::cerr<<"ERROR in getPIP: Model not implemented\n";
		exit(EXIT_FAILURE);
	}

	return model_factory;
}
//#endif
//============================================================================================================================
//#ifdef WITH_DNA
template<>
ModelFactory<DNA> *ModelFactory<DNA>::getDefault(const std::map<std::string,sequence_t<DNA> > &seqs)
{
	ModelFactory<DNA> *model_factory = NULL;

	if(cmdlineopts.cmodel_file != "") {
		std::ifstream qmat_file(cmdlineopts.cmodel_file.c_str());
		model_factory = new ModelFactoryCustom<DNA>(qmat_file);
	} else {
		error("custom model file necessary for DNA alignments");
	}

	if(cmdlineopts.aafreqs_flag) {
		std::vector<sequence_t<DNA> > seq_vector;
		for (std::map<std::string,sequence_t<DNA> >::const_iterator it2 = seqs.begin(); it2 != seqs.end(); ++it2) {
			seq_vector.push_back((*it2).second);
		}
		ModelFactoryPlusF<DNA> *plusF = new ModelFactoryPlusF<DNA>(model_factory);
		plusF->estimateFreqs(seq_vector);
		model_factory = plusF;
	}

	return model_factory;
}
//#endif
//============================================================================================================================
//DP-PIP
//#ifdef WITH_DNA
template <>
ModelFactory<DNA> *ModelFactory<DNA>::getPIP(double mu,double lambda)
{
	ModelFactory<DNA> *model_factory = NULL;

	if(cmdlineopts.cmodel_file != "") {
		std::cerr<<"ERROR in getPIP: Model not implemented\n";
		exit(EXIT_FAILURE);
	}else if(cmdlineopts.hky85_flag) {
		error("ERROR, model not implemented");
		exit(EXIT_FAILURE);

		model_factory = new ModelFactoryHKY85PIP(mu,lambda);
	}else if(cmdlineopts.k80_flag){
		model_factory = new ModelFactoryK80PIP(mu,lambda,cmdlineopts.k80);
	}else if(cmdlineopts.jc_flag){
			model_factory = new ModelFactoryJCPIP(mu,lambda,cmdlineopts.jc);
	}else {
		std::cerr<<"ERROR in getPIP: Model not implemented\n";
		exit(EXIT_FAILURE);
	}

	if(cmdlineopts.aafreqs_flag) {
		std::cerr<<"ERROR in getPIP: Model not implemented\n";
		exit(EXIT_FAILURE);
	}

	return model_factory;
}
//#endif
//============================================================================================================================
//#ifdef WITH_CODON
template<>
ModelFactory<Codon> *ModelFactory<Codon>::getDefault(const std::map<std::string,sequence_t<Codon> > &seqs)
{
	ModelFactory<Codon> *model_factory = NULL;

	if(cmdlineopts.cmodel_file != "") {
		std::ifstream qmat_file(cmdlineopts.cmodel_file.c_str());
		model_factory = new ModelFactoryCustom<Codon>(qmat_file);
	} else {
		model_factory = new ModelFactoryEcm();
	}

	if(cmdlineopts.aafreqs_flag) {
		std::vector<sequence_t<Codon> > seq_vector;
		for (std::map<std::string,sequence_t<Codon> >::const_iterator it2 = seqs.begin(); it2 != seqs.end(); ++it2) {
			seq_vector.push_back((*it2).second);
		}
		ModelFactoryPlusF<Codon> *plusF = new ModelFactoryPlusF<Codon>(model_factory);
		plusF->estimateFreqs(seq_vector);
		model_factory = plusF;
	}

	return model_factory;
}
//#endif
//============================================================================================================================
//DP-PIP
//#ifdef WITH_CODON
template <>
ModelFactory<Codon> *ModelFactory<Codon>::getPIP(double mu,double lambda)
{
	ModelFactory<Codon> *model_factory = NULL;

	return model_factory;
}
//#endif
//============================================================================================================================
