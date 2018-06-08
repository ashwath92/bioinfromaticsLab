# Algorithms_in_Bioinformatics

Interfaces and materials for the programming course "Algorithms in Bioinformatics"

INSTRUCTIONS FOR THE LAB COURSE AT THE CHAIR OF BIOINFORMATICS, ALBERT LUDWIGS UNIVERSITY OF FREIBURG:
Abstract base classes for the algorithms are located in directory prakt.
Minimal examples to get you started with your implementations are located in the base directory. Examples for pytest unittests are located in directory test. Calling ```pytest``` from the base directory will find and execute these tests. Please keep the script and class names, they are used by the tests to import your derived classes.


Description of additional files present:

The base classes provide a common interface for unittests. The test case descriptions are present in the file 'test_case_description.txt'

Note that the packages needed to run the programs are given in 'requiredPackages.txt'.
Also, that the friend.py contains a class which has a few functions which are common to 2 or more programs. There are functions to validate sequences, to parse Fasta files, and to get substitution matrix scores. 

Execution instructions:

All the main programs are located directly in the base directory, and should be called from there.

NEEDLEMAN-WUNSCH:

The Needleman Wunsch program takes 5 parameters, 4 of which are mandatory.
1. Fasta file name 1: Relative path of the first fasta file
2. Fasta file name 2: Relative path of the second fasta file
3. Substitution matrix: 2 possible values can be given: 'blosum62' or 'pam250'
4. Gap cost: cost of opening a gap
5. Complete traceback or 1 Random traceback: Optional argument which indicates whether the algorithm should return all optimal tracebacks or only a single one. Possible values:
   -c: complete traceback, i.e. all possible optimal tracebacks
   -r: 1 random optimal traceback
   Default: 1 random optimal traceback

Execution:
Given below are 3 execution strings from the command line: one with the -c flag, one with the -r flag, and the last one without any flag (default behaviour: same as -r).
The first 2 execution strings use the blosum 62 matrix, and the last one uses the pam250 matrix. (Note: this should be executed from the folder which contains needleman_wunsch.py)
Examples in Windows:
python needleman_wunsch.py \data\sequences\seq1.fasta \data\sequences\seq2.fasta blosum62 -6 -c
python needleman_wunsch.py \data\sequences\seq11.fasta \data\sequences\seq22.fasta blosum62 -6 -r
python needleman_wunsch.py \data\sequences\seq1.fasta \data\sequences\seq2.fasta pam250 -8

Examples in Unix:
python needleman_wunsch.py /data/sequences/seq1.fasta /data/sequences/seq2.fasta blosum62 -6 -c
python needleman_wunsch.py /data/sequences/seq11.fasta /data/sequences/seq22.fasta blosum62 -6 -r
python needleman_wunsch.py /data/sequences/seq1.fasta /data/sequences/seq2.fasta pam250 -6

GOTOH:

The Gotoh program takes 6 parameters, 5 of which are mandatory.

1. Fasta file name 1: Relative path of the first fasta file
2. Fasta file name 2: Relative path of the second fasta file
3. Substitution matrix: 2 possible values can be given: 'blosum62' or 'pam250'
4. Gap opening cost
5. Gap extension cost
6. Complete traceback or 1 Random traceback: Optional argument which indicates whether the algorithm should return all optimal tracebacks or only a single one. Possible values:
   -c: complete traceback, i.e. all possible optimal tracebacks
   -r: 1 random optimal traceback
   Default: 1 random optimal traceback

Execution:
Given below are 3 execution strings from the command line: one with the -c flag, one with the -r flag, and the last one without any flag (default behaviour: same as -r).
The first 2 execution strings use the blosum 62 matrix, and the last one uses the pam250 matrix. (Note: this should be executed from the folder which contains gotoh.py)
Examples in Windows:
python gotoh.py \data\sequences\seq1.fasta \data\sequences\seq2.fasta blosum62 -8 -1 -c
python gotoh.py \data\sequences\seq11.fasta \data\sequences\seq22.fasta blosum62 -8 -1 -r
python gotoh.py \data\sequences\seq1.fasta \data\sequences\seq2.fasta pam250 -8 -1

Examples in Unix:
python gotoh.py /data/sequences/seq1.fasta /data/sequences/seq2.fasta blosum62 -8 -1 -c
python gotoh.py /data/sequences/seq11.fasta /data/sequences/seq22.fasta blosum62 -8 -1 -r
python gotoh.py /data/sequences/seq1.fasta /data/sequences/seq2.fasta pam250 -8 -1


UWPGMA (UPGMA and WPGMA):
This program takes the following 4 parameters:
1. The path to the Fasta file.
2. Substitution matrix: 2 possible values can be given: 'blosum62' or 'pam250'
3. Gap cost (Needleman Wunsch linear gap)
4. Clustering: Type of phylogenetic tree clustering to be used: UPGMA or WPGMA

Execution:
2 sample execution strings each are given below for Unix and Windows. The first one in each case uses UPGMA and the second one uses WPGMA.

PLEASE NOTE that the output may change each time the programme is run because:
1. Needleman-Wunsch, which is used by the program, gives any 1 random traceback which may be different each time the program is run. 
2. The conversion of similarity to distance uses randomizatioin to calculate S_rand, which might result in slightly different distance scores each time.

Examples in Windows:
python uwpgma.py \data\sequences\multiple.fasta pam250 -8 UPGMA
python uwpgma.py \data\sequences\multiple.fasta blosum62 -6 WPGMA

Examples in Unix:
python uwpgma.py /data/sequences/multiple.fasta pam250 -8 UPGMA
python uwpgma.py /data/sequences/multiple.fasta blosum62 -6 WPGMA


Feng-Doolittle:
The Feng-Doolittle program takes the following 4 parameters (everything is identical to UWPGMA execution-wise):
1. The path to the Fasta file.
2. Substitution matrix: 2 possible values can be given: 'blosum62' or 'pam250'
3. Gap cost (Needleman Wunsch linear gap)
4. Clustering: Type of phylogenetic tree clustering to be used: UPGMA or WPGMA

Execution:
2 sample execution strings each are given below for Unix and Windows. The first one in each case uses UPGMA and the second one uses WPGMA.

PLEASE NOTE that the output may change each time the programme is run because:
1. Needleman-Wunsch, which is used by the program, gives any 1 random traceback which may be different each time the program is run. 
2. The conversion of similarity to distance uses randomizatioin to calculate S_rand, which might result in slightly different distance scores each time.

Examples in Windows:
python feng_doolittle.py \data\sequences\multiple.fasta pam250 -8 UPGMA
python feng_doolittle.py \data\sequences\multiple.fasta blosum62 -6 WPGMA

Examples in Unix:
python feng_doolittle.py /data/sequences/multiple.fasta pam250 -8 UPGMA
python feng_doolittle.py /data/sequences/multiple.fasta blosum62 -6 WPGMA


NUSSINOV:
The Nussinov program takes only 1 parameter: 
1. The path to the Fasta file.

Execution in Windows: 
python nussinov_alg.py \data\sequences\nuss_1.fasta

Execution in Unix:
python nussino_alg.py /data/sequences/nuss_1.fasta


SUM OF PAIRS:
This program was created for testing if the correct sum of pairs score is obtained for the example given in the guideline.
It takes the alignments given in the guideline by default (for PAM or BLOSUM based on the command line parameter). This
is hardcoded
The 2 command-line parameters it takes are:
1. Substitution matrix: 2 possible values can be given: 'blosum62' or 'pam250'
2. Gap cost (Needleman Wunsch linear gap)

Examples Execution:
python sumOfPairs.py pam250 -2
python sumOfPairs.py blosum62 -2
