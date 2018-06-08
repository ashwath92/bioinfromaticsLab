#-------------------------------------------------------------------------------
# Name:        Feng Doolittle Concrete class
# Purpose:     This class has functions to get a multiple sequence alignment
#              from a FASTA file containing multiple sequences. It implements
#              the popular progressive alignment algorithm devised by Feng and
#              Doolittle.
#              It imports the 'Xpgma' class to construct a guide tree, and the
#              'NeedlemanWunsch' class to get pairwise alignments. The
#              'FriendClass' class has some common functions to validate sequences,
#               get substitution matrix scores and parse FASTA files.
#
# Author:      Ashwath Sampath
#
# Created:     28-01-2018
# Copyright:   (c) Ashwath Sampath 2018
#-------------------------------------------------------------------------------
from prakt.fd import FengDoolittleBase
from needleman_wunsch import NeedlemanWunsch
from friend import FriendClass
from uwpgma import Xpgma
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from Bio import Phylo
import numpy as np
import argparse
import random
import math
import re
import os
import sys
from io import StringIO
from decimal import Decimal

@FengDoolittleBase.register
class FengDoolittle(FengDoolittleBase):
    """This concrete class implements the abstract class FengDoolittleBase
    and contains functions to execute the Feng-Doolittle algorithm."""

    def __init__(self,subsMat,gapOpenCost):
        self.subsMat = subsMat
        self.gapOpenCost = gapOpenCost
        self.MSA = []
        self.group1 = []
        self.aligned = {}

    def pairwiseAlignment(self,s1,s2,subsMat,gapOpenCost):
        # Get the similarity using the Needleman Wunsch algorithm
        nw = NeedlemanWunsch()
        gma = Xpgma()
        traceback, optimalScore = nw.buildMatrices(s1,s2,subsMat,gapOpenCost)
        alignment_strings = nw.getAlignmentsFromTracebacks(s1,s2,traceback)
        num_alignments = len(alignment_strings)
        # Get only one random optimal alignment
        randomNum = random.randint(0,num_alignments-1)
        alignment = alignment_strings[randomNum] # alignment is a list: ['','','']
        distance = gma.similarityToDistance(optimalScore,s1,s2,nw,alignment,subsMat,gapOpenCost)
        return distance,alignment

    def addEquivalentGaps(self,seq1,seq2):
        """ This function searches for gaps in seq1 and adds gaps at the
        same locations in seq2, while moving the existing characters at
        those positions forward. """
        for i,c in enumerate(seq1):
            if c == '-':
                seq2 = seq2[:i] + c + seq2[i:]
        return seq2

    def alignAndCombineGroups(self,group1,group2):
        """This function aligns 2 groups/a sequence and a group/2 sequences"""
        nw=NeedlemanWunsch()
        gma=Xpgma()
        minDistAlignment = ["","",""]
        # distance is always between 0 and 1. So, 1.5 will always be larger than
        # any distance (it is just a random number larger than the accepted range
        # of values which has no special meaning.
        minDist = 1.5
        # Get the pairwise Needleman-Wunsch alignment between every pair of sequences,
        # and find the one with the minimum distance. Note: we get 1 random alignment
        # from Needleman Wunsch here. Note: similarity score is converted to distance
        # in the pairwiseAlignment function.
        for id1,seq1 in enumerate(group1):
            for id2,seq2 in enumerate(group2):
                dist,alignment = self.pairwiseAlignment(seq1,seq2,self.subsMat,self.gapOpenCost)
                # We don't need the 3rd part of the alignment array, which indicates
                # whether it is a match or a mismatch
                del alignment[2]
                if dist < minDist:
                    minDist = dist
                    # minDistAlignment is the alignment with minimum distance.
                    # It is a list of 2 sequences.
                    minDistAlignment = alignment
                    minIndex1 = id1
                    minIndex2 = id2
        # Assign the alignments of the 2 sequences with min. distance back to their
        # original position in their respective groups.
        group1[minIndex1] = minDistAlignment[0]
        group2[minIndex2] = minDistAlignment[1]
        # After aligning one sequence in group 1 with one sequence in group 2,
        # we need to add gaps in the other sequences in the same group where
        # there is a gap in the minimum distance sequences.
        # Do this for group1
        for id1,seq1 in enumerate(group1):
            if id1 != minIndex1:
                changed1 = self.addEquivalentGaps(group1[minIndex1],seq1)
                group1[id1] = changed1
        # and for group2 as well.
        for id2,seq2 in enumerate(group2):
            if id2 != minIndex2:
                changed2  = self.addEquivalentGaps(group2[minIndex2],seq2)
                group2[id2] = changed2

        # Combine all the sequences into a single group (extend group 1).
        group1.extend(group2)
        return group1

    def createGroups(self,groups,seqClusterMap):
        """ This function forms groups out of the individual sequences
        and calls a function to align these groups together"""

        # if the groups varible contains 2 leaf nodes from the Newick tree, (nodes are
        # C1,C2,...), do the following. Eg. groups = (C1,C2)
        #
        # We can count the number of sequences by simply counting the number of Cs.
        if groups.count('C') == 2:
            # gp1 and gp2 are sequence lists.
            gp1 = []
            gp2 = []
            # Split the sequences after stripping off (, ) and ,.
            # groupList[0] and [1] are 1st and 2nd Newick tree (leaf) node names resp.
            # g1 and g2 are the corresponding sequences, obtained from the
            # dictionary seqClusterMap.
            groupList = groups[1:len(groups)-1].split(',')
            g1 = seqClusterMap[groupList[0]]
            g2 = seqClusterMap[groupList[1]]
            # Put the 1st sequence g1 in list gp1 and the 2nd sequence g2 in list gp2
            gp1.append(g1)
            gp2.append(g2)

            # Align sequence in gp1 with sequence in gp2, we get an alignment of
            # 2 sequences.
            combinedAlignment = self.alignAndCombineGroups(gp1,gp2)
            # Replace all gap symbols in the combinedAlignment by the symbol X.
            combinedAlignment = [ele.replace('-','X') for ele in combinedAlignment]

            # Append the alignment to the final Multiple sequence Alignment, self.MSA
            # if this is the first set of alignments to be seen (if MSA is empty).
            if self.MSA == []:
                for ele in combinedAlignment:
                    self.MSA.append(ele)
                # Add already seen Newick tree nodes as a string key to the 'self.aligned'
                # dict, which contains nodes already aligned. The values of the dict
                # is just the current values in the self.MSA variable.
                # Delete all existing key-value pairs, start afresh each time.
                self.aligned.clear()
                self.aligned['({},{})'.format(groupList[0],groupList[1])] = self.MSA
                return
            # If there are already alignments present in self.MSA
            else:
                # The 2nd group is now a concat of the 1st and 2nd sequence nodes
                groupList[1] = '({},{})'.format(groupList[0],groupList[1])
                # The first group is made up of the nodes in the MSA string, and gp1 has
                # the corresponding sequences in a list.
                gp1 = self.MSA
                groupList[0] = ''.join([i for i in self.aligned.keys()])
                # gp2 is the list of alignments obtained before from
                # combinedAlignment before the if condition
                gp2 = combinedAlignment
                # Combine the new gp1 (existing value in MSA and gp2 (the new
                # alignment of 2 groups)
                combinedAlignment = self.alignAndCombineGroups(gp1,gp2)
                # and replace gaps by X'es
                combinedAlignment = [ele.replace('-','X') for ele in combinedAlignment]
                # Overwrite the existing MSA list with the new alignment
                self.MSA.clear()
                for ele in combinedAlignment:
                    self.MSA.append(ele)
                # Similarly, erase the aligned dictionary, and put in the new key.
                self.aligned.clear()
                self.aligned['({},{})'.format(groupList[0],groupList[1])] = self.MSA
                return

        # There are 3 or more sequences
        if groups.count('C') > 2:
            # If the group in the groups string already exists in aligned, don't
            # do anything -- these groups are already aligned.
            if self.aligned.get(groups) != None:
                return
            # If the groups are not aligned
            # gp1 and gp2 are sequence lists.
            gp1 = []
            gp2 = []
            groupList = []

            # Case where the right group has already been aligned. Get the leftmost
            # comma and form 2 groups: before and after the comma. They are inserted
            # into groupList. Here, the left group will be a single sequence.
            leftmostComma = groups.find(',')
            groupList.append(groups[1:leftmostComma])
            groupList.append(groups[leftmostComma + 1:groups.rfind(')')])
            # The following if condchecks if the right-hand side group already
            # exists in self.aligned.
            if self.aligned.get(groupList[1]) != None:
                # Get the sequence corresponding to the left hand side node in g1,
                # and append it to the list gp1. Get gp2 directly from aligned.
                g1 = seqClusterMap[groupList[0]]
                gp1.append(g1)
                gp2 = self.aligned.get(groupList[1])
                # Get the combined alignment, replace gaps with X'es
                combinedAlignment = self.alignAndCombineGroups(gp1,gp2)
                combinedAlignment = [ele.replace('-','X') for ele in combinedAlignment]
                # Erase the MSA list and put in the new alignment
                self.MSA.clear()
                for ele in combinedAlignment:
                    self.MSA.append(ele)
                # Similarly, erase the aligned dictionary, and put in the new key.
                self.aligned.clear()
                self.aligned['({},{})'.format(groupList[0],groupList[1])] = self.MSA
                return

            # Case where the left group has already been aligned. Get the rightmost
            # comma and form 2 groups: before and after the comma. They are inserted
            # into groupList. Here, the right group will be a single sequence.
            groupList = []
            rightMostComma = groups.rfind(',')
            # The group before the comoma is in groupList[0], the 1st group.
            groupList.append(groups[1:rightMostComma])
            # The node after the comma is in groupList[1]
            groupList.append(groups[rightMostComma + 1:groups.rfind(')')])
            # The following if condition checks if the left-hand side group already
            # exists in self.aligned.
            if self.aligned.get(groupList[0]) != None:
                # Get the sequence corresponding to the right hand side node in g2,
                # (from the seqClusterMap dictionary) and append it to the list gp2.
                # Get gp1 directly from aligned.
                g2 = seqClusterMap[groupList[1]]
                gp2.append(g2)
                gp1 = self.aligned.get(groupList[0])
                # Get the combined alignment, replace gaps with X'es
                combinedAlignment = self.alignAndCombineGroups(gp1,gp2)
                combinedAlignment = [ele.replace('-','X') for ele in combinedAlignment]
                # Overwrite the MSA list with the new alignments (all alignments so far
                # including the latest one)
                self.MSA.clear()
                for ele in combinedAlignment:
                    self.MSA.append(ele)
                # Similarly, erase the aligned dictionary, and put in the new key.
                self.aligned.clear()
                self.aligned['({},{})'.format(groupList[0],groupList[1])] = self.MSA
                return

    def processNewickString(self,newick,seqClusterMap):
        """ This function reads the Newick string without distances, and
        calls the function to create groups"""
        openBrackets = closeBrackets = 0
        openInd = []
        closeInd = []
        # Iterate through the Newick string
        for i,item in enumerate(newick):
            # ( indicates item a pairing begins
            if item == '(':
                openBrackets += 1
                # index of ( is appended
                openInd.append(i)
            # ) indicates the end of a pairing
            if item == ')':
                closeBrackets += 1
                openBrackets -= 1
                # Get the index of the corresponding (
                start = openInd[openBrackets]
                end = i
                # 'groups' contains the string from ( upto )
                groups = newick[start:end+1]

                del(openInd[openBrackets])
                # Use the first group to create further groups and combine them
                # until all sequences are aligned together
                self.createGroups(groups,seqClusterMap)

    def sumOfPairs(self):
        """ This function gets the sum of pairs score of the multiple sequence
        alignment."""
        sumOfPairs = 0
        # Calls the FriendClass's getSubsMatScore function.
        fr = FriendClass()
        for i,alignment1 in enumerate(self.MSA):
            for j,alignment2 in enumerate(self.MSA):
                if i == j or i > j:
                    continue
                else:
                    # All the alignments are of the same length.
                    for index in range(len(alignment2)):
                        if alignment1[index] == 'X' and alignment2[index] == 'X':
                            continue
                        elif ((alignment1[index] == 'X' and alignment2[index] != 'X')
                        or (alignment1[index] != 'X' and alignment2[index] == 'X')):
                            sumOfPairs += self.gapOpenCost
                        # if both of the alignments have amino-acid characters
                        else:
                            sumOfPairs += fr.getSubsMatScore(alignment1[index],alignment2[index],self.subsMat,self.gapOpenCost)
        return sumOfPairs

    def printer(self,newick,SOP):
        """ This function prints the multiple sequence alignment and its sum of
        pairs score."""
        print("The guide tree in Newick format is: {}".format(newick))

        self.MSA = [ele.replace('X','-') for ele in self.MSA]
        for i,ele in enumerate(self.MSA):
            print ('{} <- Sequence {}'.format(ele,i+1))

        print("\nThe sum of pairs score of the multiple sequence alignment is {}".format(SOP))
    def run(self,
            seq_fasta_file,
            subst_matrix_fn,
            cost_gap_open,
            clustering):
            """This is the main function which parses the fasta file,
            calls functions to create the UPGMA and WPGMA trees, and
            calls another function to print the final result"""

            self.subsMat = subst_matrix_fn
            self.gapOpenCost = cost_gap_open
            fr = FriendClass()
            # Parse the fasta file. Get 2 sequences out of them
            # record is a list containing the sequences and ids in
            # the fasta file.
            record = fr.parseMultSequenceFastaFile(seq_fasta_file)
            # if there is a problem with the fasta file, list(SeqIO.parse) returns an empty list
            if len(record) == 0:
                print("You have a problem with your FASTA file.  Hint: check if the first character is '>'")
                sys.exit(1) # error code 1

            # If gap cost is positive, take the additive inverse, return the
            # negative version of the same value.
            if cost_gap_open > 0:
                print ("Your gap cost is positive. I assume you want it to be negative, I have added a minus")
                cost_gap_open = - cost_gap_open

            # The number of sequences is obtained from the length of the list 'record'
            num_sequences = len(record)
            # Get the pairwise similarities using Needleman-Wunsch.

            ids = []
            s = []
            for i in range(0,num_sequences):
                ids.append(record[i].id)
                s.append(str(record[i].seq)) # convert from Bio.Seq.Seq to str
                # Make sure s doesn't contain non-amino acid characters
                if fr.validateAminoSequence(s[i]) == 0:
                    print("You have invalid character(s) in your FASTA file")
                    sys.exit(11)         # error code 11

            # Call the function in Xpgma which will call other functions to create
            # a UPGMA/WPGMA tree based on the 'clustering' argument that is sent.
            # 2 Newick format outputs string are returned: with and without distances.
            # Only the one without distances will be used. The 3rd Newick format,
            # contains the original IDs in the fasta file, and is only needed for display
            gma = Xpgma()
            newick, newickNoDistance,distanceMatrix,newickIds = gma.UandWpgma(ids,s,num_sequences,subst_matrix_fn,cost_gap_open,clustering)
            # seqClusterMap is a dict with cluster names as keys and the
            # corresponding sequences as values.
            seqClusterMap = {}
            cl = 0
            for seq in s:
                seqClusterMap['C'+str(cl)] = seq
                cl += 1
            # Call a function which will read and parse the Newick string, and
            # will internally call other functions to create groups and get
            # the final multiple sequence alignment.
            self.processNewickString(newickNoDistance,seqClusterMap)
            SOP = self.sumOfPairs()
            self.printer(newickIds,SOP)
            return SOP, newick,newickNoDistance

if __name__ == '__main__':
    # Run Needleman-Wunsch with some parameters
    # Assign a dummy value to fd.subsMat and fd.gapOpenCost, these will be
    # overwritten in the run function.
    fd = FengDoolittle('dummy',-1)
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help='Specify the path to the'
    ' fasta file + the file name')
    parser.add_argument("subsMatrixType", choices=["pam250", "blosum62"],
        help="Choose if you want to use a PAM250 or BLOSUM62 substitution"
        " matrix for calculating match/mismatch score")
    parser.add_argument("gapOpenCost", type=int,
        help="Specify the cost of opening a gap")
    parser.add_argument("clustering", choices=["UPGMA","WPGMA"],
        help="Choose if you want to use a PAM250 or BLOSUM62 substitution"
        " matrix for calculating match/mismatch score")
    args=parser.parse_args()
    SOP, newick,newickNoDistance = fd.run(args.fasta, args.subsMatrixType, args.gapOpenCost, args.clustering.upper())
