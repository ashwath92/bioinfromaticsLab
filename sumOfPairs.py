#-------------------------------------------------------------------------------
# Name:        sumOfPairs
# Purpose:     Used to test if the sumOfPairs function used in Feng-Doolittle
#              works as expected for inputs in testing guideline
# Author:      Ashwath Sampath
#
# Created:     09-02-2018
# Copyright:   (c) Ashwath Sampath 2018
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from friend import FriendClass
import argparse

class sumOfPairs():
    """This is a standalone class for calculating the sum of pairs score"""
    def sumOfPairsScore(self,MSA,gapOpenCost,subsMat):
        """ This function calculates the sum of pairs score of the multiple sequence
        alignment given in the first argument as a list."""
        sumOfPairs = 0
        # Calls the getSubsMatScore function of the FriendClass to calculate
        # substitution matrix scores.
        for i,alignment1 in enumerate(MSA):
            for j,alignment2 in enumerate(MSA):
                if i == j or i > j:
                    continue
                else:
                    # All the alignments are of the same length.
                    for index in range(len(alignment2)):
                        if alignment1[index] == 'X' and alignment2[index] == 'X':
                            continue
                        elif ((alignment1[index] == 'X' and alignment2[index] != 'X')
                        or (alignment1[index] != 'X' and alignment2[index] == 'X')):
                            sumOfPairs += gapOpenCost
                        # if both of the alignments have amino-acid characters
                        else:
                            sumOfPairs += fr.getSubsMatScore(alignment1[index],alignment2[index],subsMat,gapOpenCost)
        return sumOfPairs

if __name__ == '__main__':

    sop = sumOfPairs()
    fr = FriendClass()
    parser = argparse.ArgumentParser()
    parser.add_argument("subsMatrixType", choices=["pam250", "blosum62"],
        help="Choose if you want to use a PAM250 or BLOSUM62 substitution"
        " matrix for calculating match/mismatch score")
    parser.add_argument("gapOpenCost", type=int,
        help="Specify the cost of opening a gap")
    args=parser.parse_args()

    # Hard code the alignments for testing, as Needleman Wunsch otherwise
    # gives randomized alignments.
    MSA_PAM = ["---MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA",
    "MTAMEESQSDISLELPLSQETFSGLWKLLPPEDIL-PSP-HCMDDLLL-PQDVEEFF-E--G----P--SE-A",
    "-----EP--------PLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEEFF-E--G----P--SE-A"]

    MSA_BLOSUM = ["---MEEPQSDPSVEPPLSQETFSDLWKLL-PENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA",
    "MTAMEESQSDISLELPLSQETFSGLWKLLPPE-DIL-PSP-HCMDDLLL-PQDVEEFF-E--GP--S----E-A",
    "-----EP--------PLSQETFSDLWKLL-PENNVLSPLPSQAMDDLMLSPDDIEEFF-E--GP--S----E-A"]

    MSA = MSA_PAM if args.subsMatrixType == "pam250" else MSA_BLOSUM
    MSA =[ele.replace('-','X') for ele in MSA]
    print("Substitution matrix:",args.subsMatrixType)
    print ("Gap cost:", args.gapOpenCost)
    print("Sum of pairs score:", sop.sumOfPairsScore(MSA,args.gapOpenCost, args.subsMatrixType))

