#-------------------------------------------------------------------------------
# Name:         Needleman Wunsch Concrete class
# Purpose:      This class has functions to implement the Needleman-Wunsch
#               algorithm, which uses linear gap costs.
#               It returns the pairwise alignment of 2 sequences from 2 different
#               FASTA files. It expects 1 sequence per file.
#               It imports the 'FriendClass' class, has some common functions to
#               validate sequences, get substitution matrix scores and parse FASTA files.
#
# Author:      Ashwath Sampath
#
# Created:     08-12-2017
# Copyright:   (c) Ashwath Sampath 2018
#-------------------------------------------------------------------------------


from prakt.nw import NeedlemanWunschBase
from friend import FriendClass
from Bio import SeqIO
import numpy as np
import argparse
import random
import copy
import os
import sys

@NeedlemanWunschBase.register
class NeedlemanWunsch(NeedlemanWunschBase):
    """This concrete class implements the abstract class NeedlemanWunschBase
    and contains functions to execute the whole Needleman-Wunsch algorithm."""

    def buildMatrices(self, s1,s2,subst_matrix_fn,gap_cost):
        """ This function creates the Needleman-Wunsch matrix, taking the sequences,
        the type of substitution matrix and the gap opening cost as arguments. It also
        cretes a traceback matrix which can be used in a later function to compute the
        optimal alignments"""
        s1_length = len(s1)
        s2_length = len(s2)
        nw_matrix = np.zeros((s1_length + 1, s2_length + 1), dtype=int)
        traceback = np.zeros((s1_length, s2_length),dtype=int)
        fr = FriendClass()

        #Initialize
        for i in range(1, s1_length+1):
            nw_matrix[i,0] = nw_matrix[i-1,0] + gap_cost
        for j in range(1, s2_length+1):
            nw_matrix[0,j] = nw_matrix[0,j-1] + gap_cost
        # sequence 1 is on the left and sequence 2 is on top
        for i in range(1, s1_length+1):
            for j in range(1, s2_length+1):
                # Cost of inserting a gap into seq 1
                seq1_gap = nw_matrix[i, j-1] + gap_cost
                # Cost of inserting a gap into seq 2
                seq2_gap = nw_matrix[i-1, j] + gap_cost
                # Cost of a match/mismatch
                # i-1, j-1 to index the strings as the i and j loops start with value 1
                substitution = nw_matrix[i-1, j-1] + fr.getSubsMatScore(s1[i-1], s2[j-1], subst_matrix_fn,gap_cost)

                nw_matrix[i][j] = max(seq1_gap, seq2_gap, substitution)
                # Store which direction we came from, we need this for traceback
                # traceback is a s1.length x s2.length matrix, so we need to index
                # from [0][0], so we use [i-1][j-1]
                """ We add 1 whenever the value was caluclated from seq1_gap, we add 2
                when the value was calculated from seq2_gap, and add 4 when the value
                was calculated from a substitution. We get the values 5,6,7 when the
                value came from 2 or 3 directions (i.e. combinations of seq1_gap, seq2_gap
                and substitutions (1,2 and 4). Note that there are three ifs and not elifs,
                so all 3 have conditions are checked, and the values are added."""

                if seq1_gap == max(seq1_gap, seq2_gap, substitution):
                    traceback[i-1][j-1] += 1
                if seq2_gap == max(seq1_gap, seq2_gap, substitution):
                    traceback[i-1][j-1] += 2
                if substitution == max(seq1_gap, seq2_gap, substitution):
                    traceback[i-1][j-1] += 4
        optimalScore = nw_matrix[s1_length][s2_length]
        return traceback, optimalScore


    def getAlignmentsFromTracebacks(self,s1,s2,traceback):
        """This function takes as input the traceback matrix with values from 1 to 7 indicating the direction from which
        a cell is calculated. (1=from left, 2=from top, 3=from left and top, 4=from diagonal, 5=from left and diagonal,
        6=from top and diagonal, 7=from diagonal, top and left. It returns a list of lists containing the alignment."""

        indices_list=[[]]
        trace_list =[[]]
        i = len(traceback) - 1
        j = len(traceback[0]) - 1

        indices_list[0] = [i,j]
        trace_list[0] = ["","",""]
        indices_duplicate = copy.deepcopy(indices_list) # A copy of indices list is needed for going through the for loop below
        #completed_counter = 0
        while True:
            completed_counter = 0 #This counter will be set to the number of tracebacks found.
            for index, [i,j] in enumerate(indices_duplicate):

                if i == -1 and j == -1:
                    # We reach here only when we have got the complete sequence
                    completed_counter +=1 #increment indicates that we have got 1 more complete traceback
                    continue

                if i ==-1 and j >= 0:
                    # We reach here only when s1 has reached the beginning of the sequence
                    trace_list[index][0] += '-'
                    trace_list[index][1] += s2[j]
                    trace_list[index][2] += ' '
                    indices_list[index][1] -= 1
                    continue

                if i >= 0 and j == -1:
                    # We reach here only when s2 has reached the beginning of the sequence
                    trace_list[index][0] += s1[i]
                    trace_list[index][1] += '-'
                    trace_list[index][2] += ' '
                    indices_list[index][0] -= 1
                    continue

                if traceback[i][j] == 1:                #case 1: we get the traceback from the left
                    trace_list[index][0] += '-'
                    trace_list[index][1] += s2[j]
                    trace_list[index][2] += ' '
                    # index of i stays as it is. Check if index of j is already less than 0. If yes, don't do anything, otherwise decrement it.

                    if indices_list[index][1] >= 0:
                        indices_list[index][1] -= 1

                elif traceback[i][j] == 2:                #case 2: we get the traceback from top
                    trace_list[index][0] += s1[i]
                    trace_list[index][1] += '-'
                    trace_list[index][2] += ' '
                    # index of j stays as it is. Check if index of i is already less than 0. If yes, don't do anything, otherwise decrement it.

                    if indices_list[index][0] >= 0:
                        indices_list[index][0] -= 1


                elif traceback[i][j] == 3:                #case 3: we get 2 tracebacks: from top and left
                    # we need to split the traceback sublist (and indices sublist) into 2 equal lists. deepcopy is used, because the
                    # normal shallow copy will result in both copies being updated whenever one is updated.
                    trace_list.append(copy.deepcopy(trace_list[index]))
                    indices_list.append(copy.deepcopy(indices_list[index]))
                    # treat traceback[index] as the list where the traceback has come from the left (decrease j)
                    trace_list[index][0] += '-'
                    trace_list[index][1] += s2[j]
                    trace_list[index][2] += ' '

                    # treat traceback[second] as the list where the traceback has come from the top (decrease i)
                    # second will store the index of the newly duplicated list (it will always be at the end because that's how append works)
                    second = len(trace_list) - 1
                    trace_list[second][0] += s1[i]
                    trace_list[second][1] += '-'
                    trace_list[second][2] += ' '

                    if indices_list[index][1] >= 0:
                        indices_list[index][1] -= 1
                    if indices_list[second][0] >= 0:
                        indices_list[second][0] -= 1

                elif traceback[i][j] == 4:                #case 4: we get the traceback from the diagonal
                    trace_list[index][0] += s1[i]
                    trace_list[index][1] += s2[j]
                    if s1[i] == s2[j]:
                        trace_list[index][2] += '*'
                    else:
                        trace_list[index][2] += ':'

                    # index of j and i need to be decremented if they are not already less than 0

                    if indices_list[index][0] >= 0 and indices_list[index][1] >= 0:
                        indices_list[index][0] -= 1
                        indices_list[index][1] -= 1


                elif traceback[i][j] == 5:                #case 3: we get 2 tracebacks: from diagonal and left
                    # split trace_list and indices_list into 2 equal lists
                    trace_list.append(copy.deepcopy(trace_list[index]))
                    indices_list.append(copy.deepcopy(indices_list[index]))
                    # treat traceback[index] as the list where the traceback has come from the left (decrease j)
                    trace_list[index][0] += '-'
                    trace_list[index][1] += s2[j]
                    trace_list[index][2] += ' '

                    # treat traceback[second] as the list where the traceback has come from the diagonal (decrease i and j)
                    #second will store the index of the newly duplicated list (it will always be at the end because that's how append works)
                    second = len(trace_list) - 1
                    trace_list[second][0] += s1[i]
                    trace_list[second][1] += s2[j]
                    if s1[i] == s2[j]:
                        trace_list[second][2] += '*'
                    else:
                        trace_list[second][2] += ':'

                    if indices_list[index][1] >= 0:
                        indices_list[index][1] -= 1
                    if indices_list[second][0] >= 0 and indices_list[second][1] >= 0:
                        indices_list[second][0] -= 1
                        indices_list[second][1] -= 1

                elif traceback[i][j] == 6:                #case 6: we get 2 tracebacks: from top and diagonal
                    # split trace_list and indices_list into 2 equal lists
                    trace_list.append(copy.deepcopy(trace_list[index])) # we need to split the traceback sublist into 2 lists
                    indices_list.append(copy.deepcopy(indices_list[index]))
                    # treat traceback[index] as the list where the traceback has come from the top (decrease i)
                    trace_list[index][0] += s1[i]
                    trace_list[index][1] += '-'
                    trace_list[index][2] += ' '

                    # treat traceback[second] as the list where the traceback has come from the left (decrease j)
                    #second will store the index of the newly duplicated list (it will always be at the end because that's how append works)
                    second = len(trace_list) - 1
                    trace_list[second][0] += s1[i]
                    trace_list[second][1] += s2[j]
                    if s1[i] == s2[j]:
                        trace_list[second][2] += '*'
                    else:
                        trace_list[second][2] += ':'

                    if indices_list[index][1] >= 0:
                        indices_list[index][0] -= 1
                    if indices_list[second][0] >= 0 and indices_list[second][0] >= 0:
                        indices_list[second][0] -= 1
                        indices_list[second][1] -= 1

                elif traceback[i][j] == 7:                #case 6: we get 3 tracebacks: from top, left and diagonal
                    # split trace_list and indices_list into 3 equal lists
                    trace_list.append(copy.deepcopy(trace_list[index])) # we need to split the trace_list sublist into 3 lists
                    trace_list.append(copy.deepcopy(trace_list[index]))
                    indices_list.append(copy.deepcopy(indices_list[index])) #first copy
                    indices_list.append(copy.deepcopy(indices_list[index])) #second copy
                    # treat traceback[index] as the list where the traceback has come from the top (decrease i)
                    trace_list[index][0] += s1[i]
                    trace_list[index][1] += '-'
                    trace_list[index][2] += ' '

                    # treat traceback[second] as the list where the traceback has come from the left (decrease j)
                    #second will store the index of the newly duplicated list (it will always be at the end because that's how append works)
                    second = len(trace_list) - 1
                    trace_list[second][0] += '-'
                    trace_list[second][1] += s2[j]
                    trace_list[second][2] += ' '

                    # treat traceback[third] as the list where the traceback has come from the diagonal (decrease i and j)
                    third = len(trace_list) - 2
                    trace_list[third][0] += s1[i]
                    trace_list[third][1] += s2[j]
                    if s1[i] == s2[j]:
                        trace_list[third][2] += '*'
                    else:
                        trace_list[third][2] += ':'

                    if indices_list[index][0] >= 0:
                        indices_list[index][0] -= 1
                    if indices_list[second][1] >= 0:
                        indices_list[second][1] -= 1
                    if indices_list[third][0] >= 0 and indices_list[third][1] >= 0:
                        indices_list[third][0] -= 1
                        indices_list[third][1] -= 1

            # indices_duplicate, the for loop variable, needs to store the updated value of indices_list before the next loop starts
            indices_duplicate = copy.deepcopy(indices_list)
            # when the number of indices (same as no. of tracebacks) is equal to the 'done counter', which is incremented once for
            # each traceback, we can break out of the while(True) infinite loop
            if completed_counter == len(indices_duplicate):
                break
        # As trace_list contains all the strings (S1, S2 and connect) in the opposite order, they need to be reversed.
        alignment_strings = [[string[::-1] for string in trace] for trace in trace_list]
        return alignment_strings

    def printer(self,alignments,num_alignments, optimalScore,complete_traceback,id1,id2,s1,s2,subsMat):
        """This function prints the alignment of string 1, a string showing * for matches and : for mismatches, and the
        alignment of string along with the optimal score. It takes a list of 3 strings (alignment of string1, alignment
        of string2 and the string for indicating match/mismatches),  the number of alignments, the optimal score,
        the complete_traceback flag vairable and the ids and sequence of the input sequences."""


        # Print only 80 characters per line from the input file (just to follow the Fasta specification).
        # Description of method used is given below (before printing the alignments).
        lineSize = 80
        len_s1 = len(s1)
        len_s2 = len(s2)

        print("Input sequences:")
        numLinesS1 = ((len_s1 // lineSize) + 1) if len_s1 % lineSize != 0 else (len_s1 // lineSize)
        for line in range(0,numLinesS1):
            print("{} \t<- Sequence 1 ({})".format(s1[lineSize * line : lineSize * (line + 1)],id1))
            if line != numLinesS1 - 1:
                print("Continued below \n")
        numLinesS2 = ((len_s2 // lineSize) + 1) if len_s2 % lineSize != 0 else (len_s2 // lineSize)
        for line in range(0,numLinesS2):
            print("{} \t<- Sequence 2 ({})".format(s2[lineSize * line : lineSize * (line + 1)],id2))
            if line != numLinesS2 - 1:
                print("Continued below \n")
        print("Substitution Matrix used: {}".format(subsMat))
        print("Number of optimal alignments between sequence 1 (with id {}) and sequence 2 (with id {}) is:{}".format(id1,id2,num_alignments))
        for index, alignment in enumerate(alignments):
            if complete_traceback == False:
                print ("One Random alignment (with alignment score {}):".format(optimalScore))
            else:
                print("Alignment {} (with alignment score {}):\n".format(index+1,optimalScore))

            # We need to print only 80 characters per line.
            # If number of characters in the ouput is exactly 80 (or a multiple of 80), the number of lines is the integer division
            # of length//linesize, for eg. for 80=> 80//80 = 1. If the number of characters in the ouptut is not exactly 80 or a
            # multiple of 80, number of lines is the integer division of length//linesize + 1. For eg, 159 => 159//80 + 1 = 2
            # This is implemented here using the ternary operator.

            alignmentLength = len(alignment[0])
            numLines = ((len(alignment[0]) // lineSize) + 1) if alignmentLength % lineSize != 0 else (len(alignment[0]) // lineSize)
            for line in range(0,numLines):
                print(alignment[0][lineSize * line : lineSize * (line + 1)] + "\t<- Sequence 1 ("+ id1+")")
                print(alignment[2][lineSize * line : lineSize * (line + 1)])
                print(alignment[1][lineSize * line : lineSize * (line + 1)] + "\t<- Sequence 2 ("+ id2+")")
                if line != numLines - 1:
                    print("Continued below \n")

    def run(self,
            seq1_fasta_file,
            seq2_fasta_file,
            subst_matrix_fn,
            cost_gap_open,
            complete_traceback):
            """This is the main function which parses fasta files,
            calls functions to create Needleman Wunsch and traceback
            matrices, and calls another function to print the final result"""

            fr = FriendClass()
            # Parse the fasta files. Get sequences out of them,
            # record1 and record2 are lists containing the sequences and ids in
            # fasta file 1 and 2 respectively.
            record1,record2 = fr.parseFastaFiles(seq1_fasta_file,seq2_fasta_file)
            # if there is a problem with the fasta files, list(SeqIO.parse) returns an empty list
            if len(record1) == 0 or len(record2) == 0:
                print("You have a problem with one of your FASTA files.  Hint: check if the first character is '>'")
                sys.exit(1) # error code 1
            id1 = record1[0].id
            s1 = str(record1[0].seq) # convert from Bio.Seq.Seq to str
            # Make sure s1 doesn't contain non-amino acid characters
            fr = FriendClass()
            if fr.validateAminoSequence(s1) == 0:
                print("You have invalid character(s) in your 1st file")
                sys.exit(11)         # error code 11
            id2 = record2[0].id
            s2 = str(record2[0].seq) # convert from Bio.Seq.Seq to str
            #Make sure s2 doesn't contain any non-amino acid characters
            if fr.validateAminoSequence(s2) == 0:
                print("You have invalid character(s) in your 2nd file")
                sys.exit(12)        # error code 12

            # If gap cost is positive, take the additive inverse, return the
            # negative version of the same value.
            if cost_gap_open > 0:
                print ("Your gap cost is positive. I assume you want it to be negative, I have added a minus")
                cost_gap_open = - cost_gap_open

            (traceback, optimalScore) = self.buildMatrices(s1,s2,subst_matrix_fn,cost_gap_open)
            alignment_strings = self.getAlignmentsFromTracebacks(s1,s2,traceback)
            num_alignments = len(alignment_strings)

            if complete_traceback == False:
                randomNum = random.randint(0,num_alignments-1)
                alignment_strings = [alignment_strings[randomNum]]

            # Call a function which prints the 3 strings on the console
            self.printer(alignment_strings,num_alignments,optimalScore,complete_traceback,id1,id2,s1,s2,subst_matrix_fn)
            return (id1,
                    s1,
                    id2,
                    s2,
                    optimalScore,
                    alignment_strings,
                    num_alignments)



if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta1", help='Specify the path to the first'
    ' fasta file + the file name')
    parser.add_argument("fasta2", help='Specify the path to the second'
    ' fasta file + the file name')
    parser.add_argument("subsMatrixType", choices=["pam250", "blosum62"],
        help="Choose if you want to use a PAM250 or BLOSUM62 substitution"
        " matrix for calculating match/mismatch score")
    parser.add_argument("gapOpenCost", type=int,
        help="Specify the cost of opening a gap")
    parser.add_argument("-c", action="store_true", default=False)
    parser.add_argument("-r", action="store_true", default=False)
    args=parser.parse_args()
    (id1, s1, id2, s2, optimalScore, alignment_strings, num_alignments) = nw.run(args.fasta1, args.fasta2,
        args.subsMatrixType, args.gapOpenCost, args.c)
