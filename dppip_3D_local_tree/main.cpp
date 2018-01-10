//===================================================================================================================
//===================================================================================================================
#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <limits>
#include "Fasta.h"
#include "Stockholm.h"
#include "ProgressiveAlignment.h"
#include "ModelFactoryWag.h"
#include "ModelFactoryDarwin.h"
#include "ModelFactoryEcm.h"
#include "ModelFactoryCustom.h"
#include "ModelFactoryPlusF.h"
#include "newick.h"
#include "profile.h"
#include "TreeNJ.h"
#include "FindRoot.h"
#include "CSProfile.h"
#include "Alphabet.h"
#include "main.h"
//===================================================================================================================
//DP-PIP
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <vector>
#include <utility>
#include "ProgressivePIP.h"
#include "gamma.h"
#include "Sequence.h"
//===================================================================================================================
//===================================================================================================================
// ugly global variable
cmdlineopts_t cmdlineopts;
//===================================================================================================================
template <class ALPHABET>
void doAlign(const std::map<std::string,std::string> &seqs_str, std::map<std::string,std::string> &out_aligned_seqs, std::vector<const PhyTree*> &out_all_trees,double &score,int ugglyflag,std::ostream* &out_score);
void parse_inut(int argc, char** argv);
double get_wall_time();
double get_cpu_time();
void write_score(std::ostream* &out_score,double score);
void write_tree(std::ostream* &out_tree,std::vector<const PhyTree *> all_trees);
void write_MSA(std::ostream *out_MSA,std::map<std::string,std::string> &aligned_seqs,std::vector<std::string> &order,std::vector<const PhyTree *> &all_trees);
bool fileExists(const std::string& filename);
void print_usage(struct rusage &usage);
std::string genUniqueFilename(std::string s);
//===================================================================================================================
//===================================================================================================================
int main(int argc, char** argv)
{
    //start timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    //get inut data
    parse_inut(argc,argv);

    std::vector<std::string> input_order;
	std::map<std::string,std::string> seqs = FastaLib(cmdlineopts.sequence_file.c_str()).readAll(input_order);

	std::map<std::string,std::string> aligned_seqs;
	std::vector<const PhyTree *> all_trees;

	//=========================================================
	//open MSA output
	std::ofstream custom_out_MSA;
	std::ostream *out_MSA;
	if(strcmp(cmdlineopts.output_file_MSA.c_str(),"/dev/null")) {
		custom_out_MSA.open(genUniqueFilename(cmdlineopts.output_file_MSA.c_str()));
		out_MSA = &custom_out_MSA;

	} else {
		out_MSA = &std::cout;
	}
	if(!*out_MSA) {
		error("error opening output file");
	}
	//=========================================================
	//open tree output
	std::ostream *out_tree;
	std::ofstream custom_out_tree;

	if(strcmp(cmdlineopts.output_file_tree.c_str(),"/dev/null")) {
		custom_out_tree.open(genUniqueFilename(cmdlineopts.output_file_tree.c_str()));
		out_tree = &custom_out_tree;

	} else {
		out_tree = &std::cout;
	}
	if(!*out_tree) {
		error("error opening output file");
	}
	//=========================================================
	//open lk output
	std::ofstream custom_out_lk;
	std::ostream *out_lk;
	if(strcmp(cmdlineopts.output_file_MSA.c_str(),"/dev/null")){
		std::string filename(cmdlineopts.output_file_MSA);
		filename.append(".lk");

		custom_out_lk.open(genUniqueFilename(filename.c_str()));
		out_lk = &custom_out_lk;

	}
	if(!*out_lk){
		error("error opening output file");
	}
	//=========================================================
	double score;

	if(cmdlineopts.dna_flag){
		doAlign<DNA>(seqs,aligned_seqs,all_trees,score,1,out_lk);
	}
	if(cmdlineopts.aa_flag){
		doAlign<AA>(seqs,aligned_seqs,all_trees,score,2,out_lk);
	}
	if(cmdlineopts.codon_flag){
		doAlign<Codon>(seqs,aligned_seqs,all_trees,score,3,out_lk);
	}

	std::vector<std::string> order = input_order;
	if(!cmdlineopts.inputorder_flag) {
		order = get_tree_order(all_trees.back());
	}

//	write_score(out_lk,score);
	write_tree(out_tree,all_trees);
	write_MSA(out_MSA,aligned_seqs,order,all_trees);

	for(std::vector<const PhyTree*>::const_iterator it=all_trees.begin(); it < all_trees.end(); ++it) {
		delete *it;
	}

	return 0;
}
//===================================================================================================================
//===================================================================================================================
template <class ALPHABET>
void doAlign(const std::map<std::string,std::string> &seqs_str,std::map<std::string,std::string> &out_aligned_seqs, std::vector<const PhyTree*> &out_all_trees,double &score,int ugglyflag,std::ostream* &out_score){

	//TODO:remove seqs
	std::vector< std::pair<std::string,sequence_t<ALPHABET> > > seqs;
	std::map<std::string,sequence_t<ALPHABET>> seq_map;

	/* used to infer the tree with NJ */
	bool pre_aligned;

	//=======================================================
	//DP-PIP
	/* remove start and stop codons */
	bool any_startStripped = false;
	bool any_endStripped = false;
	std::map<std::string,bool> startStripped;
	std::map<std::string,bool> endStripped;
	seqs=convert_to_sequence<ALPHABET>(seqs_str);
	remove_start_codon<ALPHABET>(seqs,&any_startStripped,startStripped);
	remove_stop_codon(seqs,&any_endStripped,endStripped);
	//=======================================================

	//=======================================================
	//DP-PIP
	cmdlineopts.max_dist=100.0;
	cmdlineopts.min_dist=0.0;

	/* load evolutionary model */
	ModelFactory<ALPHABET> *model_factory = ModelFactory<ALPHABET>::getPIP(cmdlineopts.muPIP,cmdlineopts.lambdaPIP);
	//=======================================================

	/* create/read initial guide tree */
	PhyTree* tree = NULL;
	PhyTree* topo = NULL;

	temporary_copy_vect_to_map(seq_map,seqs);

	if(cmdlineopts.topo_file != "") {
		std::ifstream tree_str(cmdlineopts.topo_file.c_str());
		topo = newick_parser::parse_newick(&tree_str);
	}

	if(cmdlineopts.tree_file != ""){
		std::ifstream tree_str(cmdlineopts.tree_file.c_str());
		tree = newick_parser::parse_newick(&tree_str);
	}else{

		//=================================================================
		//DP-PIP
		pre_aligned = false;
		tree = TreeNJ<ALPHABET>(seq_map,pre_aligned,model_factory,topo);
	}

	out_all_trees.push_back(tree->copy());

	//=====================================================
	//DP_PIP
	check_sequences_tree<ALPHABET>(*tree,seqs);
	//=====================================================

	//=====================================================
	tree->set_missing_node_name("V");
	//=====================================================

	//=====================================================
	//DP_PIP
	ProgressivePIPResult<ALPHABET> resultPIP;
	std::vector< ProgressivePIPResult<ALPHABET> > resultPIP_cross;
	//=====================================================

	//=====================================================
	//=====================================================
	//=====================================================
	//DP_PIP
	//@gamma distribution
	int num_cat;
	bool gamma_median=false;
	double* gamma_rr; /* substitution rates defined by the discrete gamma distribution. */
	double* gamma_r_proba; /* probabilities of the substitution rates defined by the discrete gamma distribution. */
	if(cmdlineopts.gamma_flag){
		//@gamma distribution
		num_cat=4;

		double alpha=cmdlineopts.alpha;
		double beta=cmdlineopts.beta;

		gamma_rr = new double[num_cat];
		gamma_r_proba = new double[num_cat];

        DiscreteGamma(gamma_r_proba,gamma_rr,alpha,beta,num_cat,gamma_median);

	}else{
		num_cat=1;

		gamma_rr = new double[num_cat];
		gamma_r_proba = new double[num_cat];

		gamma_rr[0]=1.0;
		gamma_r_proba[0]=1.0;
	}
	//=====================================================
	//=====================================================
	//=====================================================



	//=====================================================

	double gamma_rate;

	//@gamma distribution
	score=0.0;
	/* further rounds of alignment followed by estimation of
	 * improved tree from induced pairwise distances */
	for(int i = 0; i < cmdlineopts.iters; ++i) {

		//@gamma distribution
		// gamma distribution: num_cat > 1 otherwise 1
		for(int j=0;j<num_cat;j++){

			//@gamma distribution
			gamma_rate=gamma_rr[j];

		//=================================================
		//DP_PIP
		//@gamma distribution
		tree->initPrPIP<ALPHABET>(model_factory,gamma_rate);

		//=======================================================
		//DP-PIP
		/* normalized Poisson intensity */
		double nu;

		/* total branch lengths */
		double tau;
		if(cmdlineopts.local_tree_flag){
			tau=0.0;
			nu=0.0;
		}else{
			tau=tree->computeLength();
			nu=compute_nu(tau,model_factory->lambda,model_factory->mu);
			tree->set_tau(tau);
			tree->set_nu(nu);
			tree->set_iota(tau,model_factory->mu);
			tree->set_beta(tau,model_factory->mu);
		}
		//=======================================================

		//@gamma distribution+cross
		resultPIP_cross=compute_DP3D_PIP_tree_cross<ALPHABET>(	*tree,
																															model_factory,
																															seqs,tau,nu,
																															gamma_rate,
																															cmdlineopts.n_sub_sol,
																															cmdlineopts.local_tree_flag,
																															cmdlineopts.stoch_backtracking_flag,
																															ugglyflag,
																															out_score);

		resultPIP=resultPIP_cross.front();
		//------------------------------------------------------------------------------------------

		//@gamma distribution
		score += resultPIP.score * gamma_r_proba[j];

		//=================================================================
		}
	}

	delete tree;

	//==================================================================
	//DP_PIP
	// gamma distribution
	delete gamma_rr;
	delete gamma_r_proba;
	//==================================================================

	if(topo){
		delete topo;
	}

	delete model_factory;

	//=================================================================
	//DP-PIP
	resultPIP.MSA=convert_to_sequence_s<ALPHABET>(resultPIP.MSAs);

	/* re-insert start/stop codons */
	insert_start_codon(any_startStripped,startStripped,resultPIP.MSA);
	insert_stop_codon(any_endStripped,endStripped,resultPIP.MSA);
	temporary_copy_vect_to_map(seq_map,resultPIP.MSA);

	out_aligned_seqs=convert_to_string<ALPHABET>(seq_map);
	//=================================================================

}
//===================================================================================================================
//DP-PIP
std::string genUniqueFilename(std::string s){
	int version=1;
	std::string filename=s;
	while(fileExists(filename)){
		filename = s + "#" + std::to_string(version);
		version++;
	}

	return filename;
}
//==============================================================================
template <class ALPHABET>
double get_expected_seq_length(std::vector<std::pair<std::string,sequence_t<ALPHABET>>> &seqs){

	typedef typename std::vector<std::pair<std::string,sequence_t<ALPHABET>>>::const_iterator vector_iterator;

	double len=0.0;
	for(vector_iterator iter = seqs.begin(); iter != seqs.end(); iter++){
		len+= (iter->second).length();
	}

	return len/seqs.size();
}
//==============================================================================
//DP-PIP
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
//==============================================================================
//DP-PIP
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
//==============================================================================
//DP-PIP
void parse_inut(int argc, char** argv){

	//command line parsing
	try {

		TCLAP::CmdLine cmd(
				"ProGraphMSA, fast multiple sequence alignment", ' ', "1.0");

		TCLAP::ValueArg<std::string> outMSAputArg("o", "outputMSAfile", "Output MSA-file name",
				false, "/dev/null", "filename");
		cmd.add(outMSAputArg);

		TCLAP::ValueArg<std::string> outTreeputArg("O", "outputTreefile", "Output tree-file name",
				false, "/dev/null", "filename");
		cmd.add(outTreeputArg);

		TCLAP::ValueArg<std::string> treeArg("t", "tree", "initial guide tree", false, "", "newick file");
		cmd.add(treeArg);

		TCLAP::ValueArg<std::string> topoArg("", "topology", "topology of initial guide tree (branch lengths will be estimated)", false, "", "newick file");
		cmd.add(topoArg);

		TCLAP::SwitchArg codonArg("","codon","align DNA sequence based on a codon model");
		cmd.add(codonArg);

		TCLAP::SwitchArg dnaArg("","dna","align DNA sequence");
		cmd.add(dnaArg);

		TCLAP::SwitchArg aaArg("","aa","align AA sequence");
		cmd.add(aaArg);

		TCLAP::SwitchArg fastaArg("f","fasta","output fasta format (instead of stockholm)");
		cmd.add(fastaArg);

		TCLAP::ValueArg<double> indelArg("g", "indel_rate", "insertion/deletion rate", false, 0.0093359375, "rate");
		cmd.add(indelArg);

		TCLAP::ValueArg<double> gapextArg("e", "gap_ext", "gap extension probability", false, 0.6119140625, "probability");
		cmd.add(gapextArg);

		TCLAP::ValueArg<double> endindelArg("E", "end_indel_prob", "probability of mismatching sequence ends (set to -1 to disable this feature)", false, 0.12, "probability");
		cmd.add(endindelArg);

		TCLAP::ValueArg<double> edgehlArg("l", "edge_halflife", "edge half-life", false, 0.3, "distance");
		cmd.add(edgehlArg);

		TCLAP::ValueArg<double> altspliceArg("s", "altsplice_prob", "alternative splicing probability", false, 0.328125, "probability");
		cmd.add(altspliceArg);

		TCLAP::ValueArg<double> cutdistArg("x", "cutoff_dist", "cutoff value for pairwise distance estimation", false, 2.2, "distance");
		cmd.add(cutdistArg);

		TCLAP::ValueArg<double> mindistArg("d", "min_dist", "minimum distance for alignment", false, 0.05, "distance");
		cmd.add(mindistArg);

		TCLAP::ValueArg<double> maxdistArg("D", "max_dist", "maximum distance for alignment", false, 2.2, "distance");
		cmd.add(maxdistArg);

		TCLAP::ValueArg<double> minpdistArg("p", "min_pdist", "minimum p-distance (divergence) for alignment", false, 0.05, "distance");
		cmd.add(minpdistArg);

		TCLAP::ValueArg<double> maxpdistArg("P", "max_pdist", "maximum p-distance (divergence) for alignment", false, 0.8, "distance");
		cmd.add(maxpdistArg);

		TCLAP::SwitchArg noforcealignArg("A","no_force_align","do not force alignment of start/stop codons");
		cmd.add(noforcealignArg);

		TCLAP::ValueArg<double> reprateArg("", "repeat_indel_rate", "insertion/deletion rate for repeat units (per site)", false, 0.1, "rate");
		cmd.add(reprateArg);

		TCLAP::ValueArg<double> repextArg("", "repeat_indel_ext", "repeat indel extension probability", false, 0.3, "probability");
		cmd.add(repextArg);

		TCLAP::SwitchArg repalignArg("", "repalign", "re-align detected tandem repeat units");
		cmd.add(repalignArg);

		TCLAP::MultiSwitchArg repeatsArg("R", "repeats", "use T-Reks to identify tandem repeats");
		cmd.add(repeatsArg);

		TCLAP::ValueArg<std::string> readrepsArg("", "read_repeats", "read repeats from file", false, "", "T-Reks format output");
		cmd.add(readrepsArg);

		TCLAP::ValueArg<std::string> trdoutputArg("", "trd_output", "write TR detector output to file", false, "", "filename");
		cmd.add(trdoutputArg);

		TCLAP::ValueArg<std::string> customtrcmdArg("", "custom_tr_cmd", "custom command for detecting tandem-repeats", false, "", "command");
		cmd.add(customtrcmdArg);

		TCLAP::MultiSwitchArg rerootArg("r", "reroot", "reroot tree on all branches and minimize gap parsimony (specify twice for heuristic root search)");
		cmd.add(rerootArg);

		TCLAP::MultiSwitchArg wlsrefineArg("W","wls_refine","refine guide tree with weighted least-squares");
		cmd.add(wlsrefineArg);

		TCLAP::SwitchArg earlyrefArg("","early_refinement","perform early refinement");
		cmd.add(earlyrefArg);

		TCLAP::ValueArg<std::string> csArg("c", "cs_profile", "library of context-sensitive profiles", false, "", "file");
		cmd.add(csArg);

		TCLAP::SwitchArg darwinArg("w","darwin","use darwin's model of evolution (instead of WAG)");
		cmd.add(darwinArg);

		TCLAP::ValueArg<std::string> profileArg("", "profile_out", "output ancestral sequence profiles", false, "", "file");
		cmd.add(profileArg);

		TCLAP::ValueArg<std::string> cmodelArg("", "custom_model", "custom substitution model in qmat format", false, "", "file");
		cmd.add(cmodelArg);

		TCLAP::SwitchArg aafreqsArg("F", "estimate_aafreqs", "estimate equilibrium amino acid frequencies from input data");
		cmd.add(aafreqsArg);

		TCLAP::ValueArg<double> pcountArg("C", "aafreqs_pseudocount", "pseudo-count for estimating equilibrium amino acid frequencies", false, 1000.0, "count");
		cmd.add(pcountArg);

		TCLAP::SwitchArg nwdistArg("a","nwdist","estimate initial distance tree from NW alignments");
		cmd.add(nwdistArg);

		TCLAP::SwitchArg mldistArg("m","mldist","use ML distances");
		cmd.add(mldistArg);

		TCLAP::SwitchArg mldistgapArg("M","mldist_gap","use ML distances with gaps");
		cmd.add(mldistgapArg);

		TCLAP::SwitchArg inputorderArg("I","input_order","output sequences in input order (default: tree order)");
		cmd.add(inputorderArg);

		//TCLAP::SwitchArg onlytreeArg("T","only_tree","only output final tree instead of MSA");
		//cmd.add(onlytreeArg);

		TCLAP::ValueArg<int> iterArg("i","iterations","number of iterations re-estimating guide tree", false, 2, "iterations");
		cmd.add(iterArg);

		TCLAP::SwitchArg alltreesArg("","all_trees","output all intermediate guide trees");
		cmd.add(alltreesArg);

		TCLAP::SwitchArg ancestralArg("","ancestral_seqs","output all ancestral sequences");
		cmd.add(ancestralArg);

		TCLAP::UnlabeledValueArg<std::string> seqArg("sequences", "input sequences", true, "", "fasta file");
		cmd.add(seqArg);

		//==========================================================================================================
		//DP-PIP
		TCLAP::ValueArg<double> lambdaPIPArg("", "lambda", "insertion rate", false,1.0, "rate");
		cmd.add(lambdaPIPArg);

		TCLAP::ValueArg<double> muPIPArg("", "mu", "deletion rate", false,1.0, "rate");
		cmd.add(muPIPArg);

		//TODO: without aram
		TCLAP::ValueArg<double> hky85Arg("","hky85","TsToTvRatio",false,1.0, "ratio");
		cmd.add(hky85Arg);

		TCLAP::ValueArg<double> k80Arg("", "k80", "TsToTvRatio",false,1.0, "ratio");
		cmd.add(k80Arg);

		TCLAP::ValueArg<double> JCArg("", "jc", "rate", false,1.0, "rate");
		cmd.add(JCArg);

		TCLAP::SwitchArg PIPArg("","pip","use PIP model of evolution");
		cmd.add(PIPArg);

		TCLAP::SwitchArg GammaArg("","gamma","use gamma distribution (4 classes)");
		cmd.add(GammaArg);

		TCLAP::ValueArg<double> alphaArg("", "alpha", "alpha", false,1.0, "value");
		cmd.add(alphaArg);

		TCLAP::ValueArg<double> betaArg("", "beta", "beta", false,1.0, "value");
		cmd.add(betaArg);

		TCLAP::ValueArg<int> CrossArg("", "cross", "cross sub-optimal solutions",false,1, "num sub-optimal");
		cmd.add(CrossArg);

		TCLAP::SwitchArg StochBackArg("","SB","stochastic backtracking");
		cmd.add(StochBackArg);

		TCLAP::SwitchArg LocalTreeArg("","local","use local/global tree");
		cmd.add(LocalTreeArg);
		//==========================================================================================================

		cmd.parse(argc, argv);

		cmdlineopts.output_file_MSA = outMSAputArg.getValue();
		cmdlineopts.output_file_tree = outTreeputArg.getValue();
		cmdlineopts.sequence_file = seqArg.getValue();
		cmdlineopts.tree_file = treeArg.getValue();
		cmdlineopts.topo_file = topoArg.getValue();
		cmdlineopts.cs_file = csArg.getValue();
		cmdlineopts.cmodel_file = cmodelArg.getValue();
		cmdlineopts.profile_file = profileArg.getValue();
		cmdlineopts.fasta_flag = fastaArg.getValue();
		cmdlineopts.indel_rate = indelArg.getValue();
		cmdlineopts.end_indel_prob = endindelArg.getValue();
		cmdlineopts.gapext_prob = gapextArg.getValue();
		cmdlineopts.edge_halflife = edgehlArg.getValue();
		cmdlineopts.altsplice_prob = altspliceArg.getValue();
		cmdlineopts.pseudo_count = pcountArg.getValue();
		cmdlineopts.darwin_flag = darwinArg.getValue();
		cmdlineopts.repeat_rate = reprateArg.getValue();
		cmdlineopts.repeatext_prob = repextArg.getValue();
		cmdlineopts.repalign_flag = repalignArg.getValue();
		cmdlineopts.repeats_flag = repeatsArg.getValue();
		cmdlineopts.readreps_file = readrepsArg.getValue();
		cmdlineopts.trdout_file = trdoutputArg.getValue();
		cmdlineopts.customtr_cmd = customtrcmdArg.getValue();
		cmdlineopts.reroot_flag = rerootArg.getValue();
		cmdlineopts.earlyref_flag = earlyrefArg.getValue();
		cmdlineopts.nwdist_flag = nwdistArg.getValue();
		cmdlineopts.mldist_flag = mldistArg.getValue();
		cmdlineopts.mldist_gap_flag = mldistgapArg.getValue();
		//cmdlineopts.onlytree_flag = onlytreeArg.getValue();
		cmdlineopts.alltrees_flag = alltreesArg.getValue();
		cmdlineopts.ancestral_flag = ancestralArg.getValue();
		cmdlineopts.wlsrefine_flag = wlsrefineArg.getValue();
		cmdlineopts.aafreqs_flag = aafreqsArg.getValue();
		cmdlineopts.noforcealign_flag = noforcealignArg.getValue();
		cmdlineopts.inputorder_flag = inputorderArg.getValue();
		cmdlineopts.min_dist = mindistArg.getValue();
		cmdlineopts.max_dist = maxdistArg.getValue();
		cmdlineopts.min_pdist = minpdistArg.getValue();
		cmdlineopts.max_pdist = maxpdistArg.getValue();
		cmdlineopts.cutoff_dist = cutdistArg.getValue();
		cmdlineopts.iters = iterArg.getValue();

		//==========================================================================================================
		//DP-PIP

		cmdlineopts.hky85_flag=hky85Arg.isSet();
		cmdlineopts.k80_flag=k80Arg.isSet();
		cmdlineopts.jc_flag=JCArg.isSet();

		cmdlineopts.hky85 = hky85Arg.getValue();
		//cmdlineopts.k80_flag = k80Arg.getValue();
		//cmdlineopts.jc_flag = JCArg.getValue();

		cmdlineopts.lambdaPIP = lambdaPIPArg.getValue();
		assert(cmdlineopts.lambdaPIP > 1e-8 && "lambda > 0");

		cmdlineopts.muPIP = muPIPArg.getValue();

		assert(cmdlineopts.muPIP > 1e-8 && "mu > 0");

		//assert(cmdlineopts.lambdaPIP > cmdlineopts.muPIP && "lambda > mu");

		cmdlineopts.k80 = k80Arg.getValue();
		assert(cmdlineopts.k80 > 1e-8 && "tsTotvRatio > 0");

		cmdlineopts.jc = JCArg.getValue();
		assert(cmdlineopts.jc > 1e-8 && "jc_rate > 0");

		cmdlineopts.pip_flag = PIPArg.getValue();

		cmdlineopts.gamma_flag = GammaArg.getValue();
		cmdlineopts.alpha = alphaArg.getValue();
		cmdlineopts.beta = betaArg.getValue();

		if(CrossArg.isSet()){
			cmdlineopts.n_sub_sol  = CrossArg.getValue();

			cmdlineopts.n_sub_sol = (cmdlineopts.n_sub_sol  <=0) ? 1 : cmdlineopts.n_sub_sol ;

		}else{
			cmdlineopts.n_sub_sol  = 1;
		}

		cmdlineopts.local_tree_flag = false;
		cmdlineopts.local_tree_flag = LocalTreeArg.getValue();

		cmdlineopts.stoch_backtracking_flag = false;
		cmdlineopts.stoch_backtracking_flag = StochBackArg.getValue();

		//==========================================================================================================


		cmdlineopts.codon_flag = false;
		cmdlineopts.dna_flag = false;
		cmdlineopts.aa_flag = false;

		if(codonArg.isSet()){
			cmdlineopts.codon_flag=true;
		}else if(dnaArg.isSet()){
			cmdlineopts.dna_flag=true;
		}else if(aaArg.isSet()){
			cmdlineopts.aa_flag=true;
		}else{
			error("ERROR, no model selected ('codon', 'dna', 'aa')");
		}

//#ifdef WITH_CODON
		/* scale default parameters for codon distances */
		if(codonArg.isSet()) {
			if(!indelArg.isSet()) {
				cmdlineopts.indel_rate /= 2.6;
			}
			if(!edgehlArg.isSet()) {
				cmdlineopts.edge_halflife *= 2.6;
			}
			if(!maxdistArg.isSet()) {
				cmdlineopts.max_dist = 5.0;
			}
			if(!cutdistArg.isSet()) {
				cmdlineopts.cutoff_dist = 5.0;
			}
		}
//#endif

		/* do not iterate when guide tree provided */
		if(cmdlineopts.tree_file != "" && !iterArg.isSet()) {
			cmdlineopts.iters = 0;
		}

	}catch (TCLAP::ArgException &e){
		std::cerr << "Command line error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

}
//==============================================================================
//DP-PIP
void write_score(std::ostream* &out_score,double score){

	out_score->precision(15);
	*out_score<<score;

}
//==============================================================================
//DP-PIP
void write_tree(std::ostream* &out_tree,std::vector<const PhyTree *> all_trees){

	if(cmdlineopts.alltrees_flag) {
		for(std::vector<const PhyTree*>::const_iterator it=all_trees.begin(); it < all_trees.end(); ++it) {
			*out_tree << (*it)->formatNewick() << std::endl;
		}
	}else{
		*out_tree << all_trees.back()->formatNewick() << std::endl;
	}

}
//==============================================================================
void write_MSA(std::ostream *out_MSA,std::map<std::string,std::string> &aligned_seqs,std::vector<std::string> &order,std::vector<const PhyTree *> &all_trees){

	if(cmdlineopts.fasta_flag) {
		write_fasta(aligned_seqs,order,*out_MSA);
	} else {
		if(!cmdlineopts.alltrees_flag) {
			write_stockholm(aligned_seqs,order,*all_trees.back(),*out_MSA);
		} else {
			write_stockholm(aligned_seqs,order,*all_trees.back(),all_trees,*out_MSA);
		}
	}

}
//==============================================================================
// Function: fileExists
/**
    Check if a file exists
@param[in] filename - the name of the file to check

@return    true if the file exists, else false

*/
bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}




