This repository contains the scripts used compare predictability of intonation
in infant-directed speech (IDS) and adult-directed speech (ADS), using the
data from ManyBabies project (The ManyBabies Consortium, 2017; https://osf.io/re95x/), 
and as reported in :

   Rasanen, O. Kakouros, S. & Soderstrom, M.: 
   "Is infant-directed speech interesting because it is surprising? – 
   Linking properties of IDS to statistical learning and attention at the 
   prosodic level". Submitted for publication.

Any use or derivation of this code should cite the above manuscript. 

Please send any questions, comments or bug reports to Okko Rasanen
(firstname.surname@aalto.fi).

Date of this readme: March 22, 2018

------------------------------------------------------------------------

THIRD-PARTY COMPONENTS THAT ARE NEEDED TO RUN THE EXPERIMENTS FROM SCRACTH:

1) requires YAAPT (Zahorian & Hu, 2008) for F0 estimation. 
Version 4.0 was used for the results reported in the manuscript.
Download from: http://ws2.binghamton.edu/zahorian/yaapt.htm

2) Requires ManyBabies IDS/ADS stimuli
Download the IDS and ADS samples from https://osf.io/xbv95/files/
--> OSF Storage --> "IDS Clips Soderstrom" and "ADS Clips Soderstrom",
and add the corersponding paths to the filedir1 and filedir2 in the 
main scripts (see below). 

(optional)
3) To run the LSTM-based statistical model in addition to the mixed-order
Markov chains (MOCM), you also need Python (tested with v. 2.7) with
Keras (tested with v. 2.0.8) and Theano or TensorFlow backend (tested only 
with Theano v. 0.9.0.dev) installed. The script idsads_LSTM.py in this main 
folder is called from the MATLAB through system() functionality to carry 
out model training and likelihood estimation. 
Whether LSTM is used is determined by a flag in the main functions.

------------------------------------------------------------------------

INCLUDED ALGORITHMS THAT HAVE ALREADY BEEN PUBLISHED ELSEWHERE:

The package includes oscillator-based syllabification algorithm
published in Rasanen, Doyle & Frank (2018): Pre-linguistic
segmentation of speech into syllable-like units. Cognition, 171, 130–150. 
This algorithm is in /thetaOscillator folder
Direct download is also available at https://github.com/orasanen/thetaOscillator


The package also includes MATLAB implementation of a mixed-order 
Markov chain model (MOCM) (Saul & Pereira, 1997) written by O. Rasanen. 
This is the same implementation as used in Rasanen & Saarinen (2016):
"Sequence prediction with sparse distributed hyperdimensional coding applied 
to the analysis of mobile phone use patterns", 
IEEE Transactions on Neural Networks and Learning Systems, 27, 1878–1889.
The algorithm is in /MixedOrderMarkov folder

------------------------------------------------------------------------

HOW TO USE THESE SCRIPTS

The functions to run the main experiments reported in the paper are in this
root directory, including 

- predictability_IDS_ADS_main.m	(the main experiment)

- predictability_IDS_ADS_main_revision.m 
(the main experiment with varying proportion of IDS in the training)

Please change the path variables in those files to point to the location of
ManyBabies data in your computing environment 

(and to Python binary, for LSTM computing in case MATLAB system() does not
load proper environmental variables on your platform)

-----

The scripts used for the main result analysis and plotting can be found in folder 
/result_analysis_and_plotting/

-----

Results from the experiments are stored in /results

The result files corresponding to the results reported in the paper can be
found in named sub-folders under /results, separately for the main experiment
and separately for varying proportion of IDS data in the training which was
carried out for a revised version of the article (March, 2018).

Running the main analyses in idsads_analyze_results.m will also create
.csv files that can be analyzed outside MATLAB (in the paper, 
logistic regression experiments were conducted in SPSS).

-----

Temporary files (features) for the experiments are stored in /saves

Auxiliary functions for various stuff related to data reading, plotting,
run-time progres bars etc. are in /aux folder. 










