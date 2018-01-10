#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include "Alphabet.h"
#include "PhyTree.h"

//=======================================================================================================
//DP-PIP
template <class ALPHABET>
struct CompareFirst
{
  CompareFirst(std::string val) : val_(val) {}
  bool operator()(const std::pair<std::string,sequence_t<ALPHABET>>& elem) const {
    return val_ == elem.first;
  }
  private:
    std::string val_;
};
//=======================================================================================================
//DP-PIP

template <class ALPHABET>
void PrintMap(const std::map<sequence_t<ALPHABET>,Eigen::VectorXd>& m){

	typedef typename std::map<sequence_t<ALPHABET>,Eigen::VectorXd>::const_iterator MapIterator;
    for (MapIterator iter = m.begin(); iter != m.end(); iter++){
        std::cout << "Key: " << stringFromSequence(iter->first) << " " << "Values: " << (iter->second).transpose() << std::endl;
    }

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void PrintVector(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &result){

	for(unsigned int i=0;i<result.size();i++){
		std::cout<<"key:"<<result.at(i).first <<" val: "<<stringFromSequence(result.at(i).second);std::cout<<"\n";
	}

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::map<std::string,std::string> convert_to_string(std::map<std::string,sequence_t<ALPHABET>> &seqs){
	std::map<std::string,std::string> seqs_str;

	typedef typename std::map<std::string,sequence_t<ALPHABET>>::iterator vect_iter;

	for (vect_iter it = seqs.begin(); it != seqs.end(); ++it) {
		seqs_str[it->first]=stringFromSequence(it->second);
	}

	return seqs_str;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::vector<std::pair<std::string,sequence_t<ALPHABET>>> convert_to_sequence(const std::map<std::string,std::string> &seqs_str){
	std::vector<std::pair<std::string,sequence_t<ALPHABET>>> seqs;

	for (std::map<std::string, std::string>::const_iterator it = seqs_str.begin(); it != seqs_str.end(); ++it) {
		seqs.push_back(std::make_pair(it->first,sequenceFromStringPIP<ALPHABET>(it->second)));
	}

	return seqs;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
std::vector<std::pair<std::string,sequence_t<ALPHABET>>> convert_to_sequence_s(const std::vector<std::pair<std::string,std::string>> &seqs_str){
	std::vector<std::pair<std::string,sequence_t<ALPHABET>>> seqs;

	for (std::vector<std::pair<std::string,std::string>>::const_iterator it = seqs_str.begin(); it != seqs_str.end(); ++it) {
		seqs.push_back(std::make_pair(it->first,sequenceFromStringPIP<ALPHABET>(it->second)));
	}

	return seqs;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void temporary_copy_vect_to_map(std::map<std::string,sequence_t<ALPHABET>> &seqs,const std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &result){

	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>>::const_iterator vector_iterator;

	for(vector_iterator iter = result.begin(); iter != result.end(); iter++){
		seqs[iter->first] = iter->second;
	}

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void remove_start_codon(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &seqs,bool *any_startStripped,std::map<std::string,bool> &startStripped){



	if(!cmdlineopts.noforcealign_flag) {
		typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>>::iterator vect_iter;
		for (vect_iter iter = seqs.begin(); iter != seqs.end(); iter++){

			std::pair<std::string,sequence_t<ALPHABET>> tmp = (*iter);
			sequence_t<ALPHABET> &seq = tmp.second;

			if(ALPHABET::stripStart != ALPHABET::GAP && seq[0] == ALPHABET::stripStart) {
				seq = seq.substr(1);

				*(iter)=std::make_pair(tmp.first,seq);

				*any_startStripped = true;
				startStripped[iter->first] = true;
			} else {
				startStripped[iter->first] = false;
			}

		}


    }
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void remove_stop_codon(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &seqs,bool *any_endStripped,std::map<std::string,bool> &endStripped){

	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>>::iterator mapIterator;
    for (mapIterator iter = seqs.begin(); iter != seqs.end(); iter++){
		if(!cmdlineopts.noforcealign_flag) {

			std::pair<std::string,sequence_t<ALPHABET>> tmp = (*iter);
			sequence_t<ALPHABET> &seq = tmp.second;

			if(ALPHABET::stripEnd != ALPHABET::GAP && seq[seq.length()-1] == ALPHABET::stripEnd) {
				seq = seq.substr(0,seq.length()-1);

				*(iter)=std::make_pair(tmp.first,seq);

				*any_endStripped = true;
				endStripped[iter->first] = true;
			} else {
				endStripped[iter->first] = false;
			}

		}

	}

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void insert_start_codon(bool any_startStripped,std::map<std::string,bool> &startStripped,std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &MSA){

	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>>::iterator vector_iterator;
    for (vector_iterator iter = MSA.begin(); iter != MSA.end(); iter++){

		std::pair<std::string,sequence_t<ALPHABET>> tmp = (*iter);
		sequence_t<ALPHABET> &aseq = tmp.second;

    	if(any_startStripped) {
			if(startStripped.find(iter->first) != startStripped.end() && startStripped[iter->first]) {
				aseq.insert(aseq.begin(),ALPHABET::X);

				*(iter)=std::make_pair(tmp.first,aseq);

			} else {
				aseq.insert(aseq.begin(),ALPHABET::GAP);

				*(iter)=std::make_pair(tmp.first,aseq);
			}

		}

	}


}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void insert_stop_codon(bool any_endStripped,std::map<std::string,bool> &endStripped,std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &MSA){

	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>>::iterator vector_iterator;
    for (vector_iterator iter = MSA.begin(); iter != MSA.end(); iter++){

		std::pair<std::string,sequence_t<ALPHABET>> tmp = (*iter);
		sequence_t<ALPHABET> &aseq = tmp.second;

		if(any_endStripped) {
			if(endStripped.find(iter->first) != endStripped.end() && endStripped[iter->first]) {
				aseq.insert(aseq.end(),ALPHABET::X);

				*(iter)=std::make_pair(tmp.first,aseq);

			} else {
				aseq.insert(aseq.end(),ALPHABET::GAP);

				*(iter)=std::make_pair(tmp.first,aseq);
			}
		}

    }

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
bool check_uniform_len(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &result){
	unsigned int len;

	if(result.size()==0){
		return false;
	}

	len=result.at(0).second.length();

	for(unsigned int i=1;i<result.size();i++){
		if(len!=result.at(i).second.length()){
			return false;
		}
	}

	return true;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
bool check_uniform_len_s(std::vector<std::pair<std::string,std::string>> &result){
	unsigned int len;

	if(result.size()==0){
		return false;
	}

	len=result.at(0).second.length();

	for(unsigned int i=1;i<result.size();i++){
		if(len!=result.at(i).second.length()){
			return false;
		}
	}

	return true;
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void add_sequence_to_vector(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &result,std::string name,std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &sequences){
	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>> ::iterator vect_iter;

	vect_iter iter = std::find_if(sequences.begin(),sequences.end(),CompareFirst<ALPHABET>(name));
	result.push_back(std::make_pair(iter->first,iter->second));

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
void add_sequence_to_vector_string(std::vector< std::pair<std::string,std::string> > &result,std::string name,std::vector< std::pair<std::string,sequence_t<ALPHABET>> > &sequences){
	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>> ::iterator vect_iter;

	vect_iter iter = std::find_if(sequences.begin(),sequences.end(),CompareFirst<ALPHABET>(name));
	std::string s=stringFromSequence<ALPHABET>(iter->second);
	result.push_back(std::make_pair(iter->first,s));

}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
int get_length_seq(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &result){

	if(result.size()==0){
		return 0;
	}

	if(!check_uniform_len(result)){
		error("ERROR in get_length_seq: non aligned");
	}

	return result.at(0).second.length();
}
//=======================================================================================================
//DP-PIP
template <class ALPHABET>
int get_length_seq_s(std::vector<std::pair<std::string,std::string>> &result){

	if(result.size()==0){
		return 0;
	}

	if(!check_uniform_len_s<ALPHABET>(result)){
		error("ERROR in get_length_seq: non aligned");
	}

	return result.at(0).second.length();
}
//=======================================================================================================

#endif /* SEQUENCE_H_ */
