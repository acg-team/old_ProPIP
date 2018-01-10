# ProPIP

## Usage

USAGE:

Brief USAGE: 
   ./progressivePIP  [--local] [--SB] [--cross <num sub-optimal>] [--beta
                     <value>] [--alpha <value>] [--gamma] [--pip] [--jc
                     <rate>] [--k80 <ratio>] [--hky85 <ratio>] [--mu
                     <rate>] [--lambda <rate>] [--ancestral_seqs]
                     [--all_trees] [-i <iterations>] [-I] [-M] [-m] [-a]
                     [-C <count>] [-F] [--custom_model <file>]
                     [--profile_out <file>] [-w] [-c <file>]
                     [--early_refinement] [-W] ...  [-r] ... 
                     [--custom_tr_cmd <command>] [--trd_output <filename>]
                     [--read_repeats <T-Reks format output>] [-R] ... 
                     [--repalign] [--repeat_indel_ext <probability>]
                     [--repeat_indel_rate <rate>] [-A] [-P <distance>] [-p
                     <distance>] [-D <distance>] [-d <distance>] [-x
                     <distance>] [-s <probability>] [-l <distance>] [-E
                     <probability>] [-e <probability>] [-g <rate>] [-f]
                     [--aa] [--dna] [--codon] [--topology <newick file>]
                     [-t <newick file>] [-O <filename>] [-o <filename>]
                     [--] [--version] [-h] <fasta file>


Where: 

   --local
     use local/global tree

   --SB
     stochastic backtracking

   --cross <num sub-optimal>
     cross sub-optimal solutions

   --beta <value>
     beta

   --alpha <value>
     alpha

   --gamma
     use gamma distribution (4 classes)

   --pip
     use PIP model of evolution

   --jc <rate>
     rate

   --k80 <ratio>
     TsToTvRatio

   --hky85 <ratio>
     TsToTvRatio

   --mu <rate>
     deletion rate

   --lambda <rate>
     insertion rate

   --ancestral_seqs
     output all ancestral sequences

   --all_trees
     output all intermediate guide trees

   -i <iterations>,  --iterations <iterations>
     number of iterations re-estimating guide tree

   -I,  --input_order
     output sequences in input order (default: tree order)

   -M,  --mldist_gap
     use ML distances with gaps

   -m,  --mldist
     use ML distances

   -a,  --nwdist
     estimate initial distance tree from NW alignments

   -C <count>,  --aafreqs_pseudocount <count>
     pseudo-count for estimating equilibrium amino acid frequencies

   -F,  --estimate_aafreqs
     estimate equilibrium amino acid frequencies from input data

   --custom_model <file>
     custom substitution model in qmat format

   --profile_out <file>
     output ancestral sequence profiles

   -w,  --darwin
     use darwin's model of evolution (instead of WAG)

   -c <file>,  --cs_profile <file>
     library of context-sensitive profiles

   --early_refinement
     perform early refinement

   -W,  --wls_refine  (accepted multiple times)
     refine guide tree with weighted least-squares

   -r,  --reroot  (accepted multiple times)
     reroot tree on all branches and minimize gap parsimony (specify twice
     for heuristic root search)

   --custom_tr_cmd <command>
     custom command for detecting tandem-repeats

   --trd_output <filename>
     write TR detector output to file

   --read_repeats <T-Reks format output>
     read repeats from file

   -R,  --repeats  (accepted multiple times)
     use T-Reks to identify tandem repeats

   --repalign
     re-align detected tandem repeat units

   --repeat_indel_ext <probability>
     repeat indel extension probability

   --repeat_indel_rate <rate>
     insertion/deletion rate for repeat units (per site)

   -A,  --no_force_align
     do not force alignment of start/stop codons

   -P <distance>,  --max_pdist <distance>
     maximum p-distance (divergence) for alignment

   -p <distance>,  --min_pdist <distance>
     minimum p-distance (divergence) for alignment

   -D <distance>,  --max_dist <distance>
     maximum distance for alignment

   -d <distance>,  --min_dist <distance>
     minimum distance for alignment

   -x <distance>,  --cutoff_dist <distance>
     cutoff value for pairwise distance estimation

   -s <probability>,  --altsplice_prob <probability>
     alternative splicing probability

   -l <distance>,  --edge_halflife <distance>
     edge half-life

   -E <probability>,  --end_indel_prob <probability>
     probability of mismatching sequence ends (set to -1 to disable this
     feature)

   -e <probability>,  --gap_ext <probability>
     gap extension probability

   -g <rate>,  --indel_rate <rate>
     insertion/deletion rate

   -f,  --fasta
     output fasta format (instead of stockholm)

   --aa
     align AA sequence

   --dna
     align DNA sequence

   --codon
     align DNA sequence based on a codon model

   --topology <newick file>
     topology of initial guide tree (branch lengths will be estimated)

   -t <newick file>,  --tree <newick file>
     initial guide tree

   -O <filename>,  --outputTreefile <filename>
     Output tree-file name

   -o <filename>,  --outputMSAfile <filename>
     Output MSA-file name

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <fasta file>
     (required)  input sequences



## Building from source

For building progressivePIP from source you need:
CMake   >=2.8           (http://www.cmake.org)
tclap   >=1.1.0         (http://tclap.sourceforge.net)
Eigen   2.0.x or 3.0.x  (http://eigen.tuxfamily.org)

on Debian/Ubuntu you can install these programs/libraries with:
sudo apt-get install cmake libtclap-dev libeigen2-dev


Then perform the following command to configure/build/install progressivePIP:
cd BUILD
ccmake .. (press "c" to configure and "g" to generate the Makefile, see below
           for additional configuration options)
make progressivePIP
make install


Additional CMake configuration options (in "ccmake .."):

EIGEN2_INCLUDE_DIR: set this to the path, where Eigen is installed, if you use
                    Eigen 3.0.x or if Eigen has been installed at a non-default
                    location (default location: /usr/include/eigen2)

WITH_EIGEN3:        set this to ON, if you want to compile progressivePIP with
                    Eigen 3.0.x

CMAKE_CXX_FLAGS:    add options for the C++ compiler, like optimization flags,
                    or additional search paths for include files (-I <path>)

WITH_SSE:           disable this option, if you build progressivePIP for a machine
                    that does not support SSE2



 
