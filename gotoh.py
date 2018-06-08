#-------------------------------------------------------------------------------
# Name:         Gotoh Concrete class
# Purpose:      This class has functions to implement the Gotoh algorithm, which
#               uses affine gap costs. It returns the pairwise alignment of 2
#               sequences from 2 different FASTA files. It expects 1 sequence per file.
#               It imports the 'FriendClass' class, has some common functions to
#               validate sequences, get substitution matrix scores and parse FASTA files.
#
# Author:      Ashwath Sampath
#
# Created:     13-12-2017
# Copyright:   (c) Ashwath Sampath 2018
#-------------------------------------------------------------------------------


from prakt.gt import GotohBase
from friend import FriendClass
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import random
import copy
import argparse
import numpy as np
import sys

@GotohBase.register
class Gotoh(GotohBase):
    """This concrete class implements the abstract class GotohBase
    and contains functions to execute the whole Gotoh algorithm."""

    def buildMatrices(self, s1,s2,subst_matrix_fn,gap_open_cost, gap_extend_cost, g):
        """ This function creates the Needleman-Wunsch matrix, taking the sequences,
        the type of substitution matrix and the gap opening cost as arguments. It also
        cretes a traceback matrix which can be used in a later function to compute the
        optimal alignments"""
        s1_length = len(s1)
        s2_length = len(s2)
        # assign a high negative number to infinity, which will be used in initialization of P and Q matrices
        inf = -60000
        D = np.zeros((s1_length + 1, s2_length + 1), dtype=int)
        # P matrix is used to extend gaps in Sequence 2
        P = np.zeros((s1_length + 1, s2_length + 1), dtype=int)
        # Q matrix is used to extend gaps in Sequence 1
        Q = np.zeros((s1_length + 1, s2_length + 1), dtype=int)
        traceback = np.zeros((s1_length, s2_length),dtype=int)
        fr=FriendClass()

        #Initialize
        D[0,1] = D[1,0] = g
        P[0,1] = inf
        Q[1,0] = inf
        for i in range(2, s1_length+1):
            D[i,0] = D[i-1,0] + gap_extend_cost
            # P does not need to be initialized in the 1st column. These values are not used in the algorithm
            Q[i,0] = inf
        for j in range(2, s2_length+1):
            D[0,j] = D[0,j-1] + gap_extend_cost
            P[0,j] = inf
            # Q does not need to be initialized in the 1st column. These values are not used in the algorithm
        # sequence 1 is on the left and sequence 2 is on top
        for i in range(1, s1_length+1):
            for j in range(1, s2_length+1):
                #D_i-1,j + g
                # Update P[i,j] -> we can either extend the gap from the previous row in P or create a new gap in Seq 1, which
                # means that we need to take the previous row's value in D (We don't take into account different values of j)
                P[i,j] = max(D[i-1,j] + g, P[i-1,j] + gap_extend_cost)
                #Next, update Q[i,j] -> we can either extend the gap from the previous col in Q or create a new gap in Seq 1, which
                # means that we need to take the previous col's value in D (We don't take into account different values of i)
                Q[i,j] = max(D[i,j-1] + g, Q[i,j-1] + gap_extend_cost)
                # Finally, update D[i,j]: it is the max of the substitution score (match/mismatch), and the resp. P[i,j] and Q[i,j],
                # which correspond to gap extension in seq 2 and seq 1 respectively
                substitution = D[i-1, j-1] + fr.getSubsMatScore(s1[i-1], s2[j-1], subst_matrix_fn,gap_extend_cost)
                D[i,j] = max(substitution,P[i,j],Q[i,j])

        optimalScore = D[s1_length][s2_length]
        return D,P,Q,optimalScore

    def tracebackInP(self,P,D,i2,j2,beta,g):
        """This recursive function checks if the current P[i,j] value came from P[i-1,j] or D[i-1,j]. If it came from P[i-1,j], it
        keeps recursing till at some point, it has to have come from D[i,j]. This is the base case for the recursion.
        It returns the row index (i) in D[i,j] in this case"""
        if P[i2,j2] == P[i2-1,j2] + beta:
            i2 -= 1
            return self.tracebackInP(P,D,i2,j2,beta,g)
        if P[i2,j2] == D[i2-1,j2] + g:
            i2 -= 1
            return i2

    def tracebackInQ(self,Q,D,i2,j2,beta,g):
        """This recursive function checks if the current Q[i,j] value came from Q[i,j-1] or D[i,j-1]. If it came from Q[i,j-1], it
        keeps recursing till at some point, it is guaranteed to have come from D[i,j]. This is the base case for the recursion.
        It returns the column index (j) in D[i,j] in this case" """
        if Q[i2,j2] == Q[i2,j2-1] + beta:
            j2 -= 1
            return self.tracebackInQ(Q,D,i2,j2,beta,g)
        if Q[i2,j2] == D[i2,j2-1] + g:
            j2 -= 1
            return j2

    def getAlignmentsFromTracebacks(self,s1,s2,subst_matrix_fn,D,P,Q,alpha,beta,g):
        """This function takes as input the matrices D, P and Q created in an earlier function. It
        computes the traceback by essentially reversing the process of building the matrices.
         It returns a list of lists containing the alignment."""

        indices_list=[[]]
        trace_list =[[]]
        fr=FriendClass()
        # Set i and j to the index of the last row and column of the ndarray D respectively
        i = D.shape[0] - 1
        j = D.shape[1] - 1
        indices_list[0] = [i,j]
        trace_list[0] = ["","",""]
        indices_duplicate = copy.deepcopy(indices_list) # A copy of indices list is needed for going through the for loop below
        while True:
            completed_counter = 0 #This counter will be set to the number of tracebacks found.
            for index, [i,j] in enumerate(indices_duplicate):

                if i == 0 and j == 0:
                    # We reach here only when we have got the complete sequence
                    completed_counter +=1 #increment indicates that we have got 1 more complete traceback
                    continue

                if i == 0 and j >= 0:
                    # We reach here only when s1 has reached the beginning of the sequence
                    trace_list[index][0] += '-'
                    trace_list[index][1] += s2[j]
                    trace_list[index][2] += ' '
                    indices_list[index][1] -= 1
                    continue

                if i >= 0 and j == 0:
                    # We reach here only when s2 has reached the beginning of the sequence
                    trace_list[index][0] += s1[i]
                    trace_list[index][1] += '-'
                    trace_list[index][2] += ' '
                    indices_list[index][0] -= 1
                    continue
                # indicates that the value in D[i,j] came from P, Q and D
                if D[i,j] == P[i,j] and D[i,j] == Q[i,j] and D[i,j] == D[i-1,j-1] + fr.getSubsMatScore(s1[i-1], s2[j-1], subst_matrix_fn,beta):
                    trace_list.append(copy.deepcopy(trace_list[index])) # we need to split the trace_list sublist into 3 lists
                    trace_list.append(copy.deepcopy(trace_list[index]))
                    indices_list.append(copy.deepcopy(indices_list[index])) #first copy
                    indices_list.append(copy.deepcopy(indices_list[index])) #second copy
                    # treat traceback[index] as the list where the traceback has come from P[i,j]
                    # i2 is the index after calling tracebackInP: after traversing through the
                    # P matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length, and closing the gap.
                    i2 = self.tracebackInP(P,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_i gives the number of the gap we have to insert in sequence 2, and in
                    # trace_list[index][1] . It is also used to decrement the corresponding indices_list row index.
                    diff_i = i - i2
                    indices_list[index][0] -= diff_i
                    trace_list[index][0] += s1[i2:i]
                    trace_list[index][1] += '-' * diff_i
                    trace_list[index][2] += ' ' * diff_i

                    # treat traceback[second] as the list where the traceback has come from Q[i,j]
                    # second will store the index of the newly duplicated list
                    # (it will always be at the end because that's how append works)
                    second = len(trace_list) - 1
                    # j2 is the index after calling tracebackInQ: after traversing through the
                    # Q matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length and closing the gap.
                    j2 = self.tracebackInQ(Q,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_j gives the number of gaps we have to insert in sequence 1, and in
                    # trace_list[second][1] . It is also used to decrement the corresponding indices_list col index.
                    diff_j = j - j2
                    indices_list[second][1] -= diff_j
                    trace_list[second][0] += '-' * diff_j
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between j2 and j as we always add to
                    # trace_list in the reverse order
                    trace_list[second][1] += s2[j2:j][::-1]
                    trace_list[second][2] += ' ' * diff_j

                    # treat traceback[third] as the list where the traceback has come from D[i-1,j-1] (decrease i and j)
                    third = len(trace_list) - 2
                    trace_list[third][0] += s1[i-1]
                    trace_list[third][1] += s2[j-1]
                    if s1[i-1] == s2[j-1]:
                        trace_list[third][2] += '*'
                    else:
                        trace_list[third][2] += ':'
                    indices_list[third][0] -= 1
                    indices_list[third][1] -= 1

                # The value in D[i,j] came from P and Q, not from D.
                elif D[i,j] == P[i,j] and D[i,j] == Q[i,j]:
                    trace_list.append(copy.deepcopy(trace_list[index])) # we need to split the trace_list sublist into 3 lists
                    indices_list.append(copy.deepcopy(indices_list[index])) #copy
                    # treat traceback[index] as the list where the traceback has come from P[i,j]
                    # i2 is the index after calling tracebackInP: after traversing through the
                    # P matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length, and closing the gap.
                    i2 = self.tracebackInP(P,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_i gives the number of the gap we have to insert in sequence 2, and in
                    # trace_list[index][1] . It is also used to decrement the corresponding indices_list row index.
                    diff_i = i - i2
                    indices_list[index][0] -= diff_i
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between i2 and i as we always add to
                    # trace_list in the reverse order
                    trace_list[index][0] += s1[i2:i][::-1]
                    trace_list[index][1] += '-' * diff_i
                    trace_list[index][2] += ' ' * diff_i

                    # treat traceback[second] as the list where the traceback has come from Q[i,j]
                    # second will store the index of the newly duplicated list
                    # (it will always be at the end because that's how append works)
                    second = len(trace_list) - 1
                    # j2 is the index after calling tracebackInQ: after traversing through the
                    # Q matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length and closing the gap.
                    j2 = self.tracebackInQ(Q,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_j gives the number of gaps we have to insert in sequence 1, and in
                    # trace_list[second][1] . It is also used to decrement the corresponding indices_list col index.
                    diff_j = j - j2
                    indices_list[second][1] -= diff_j
                    trace_list[second][0] += '-' * diff_j
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between j2 and j as we always add to
                    # trace_list in the reverse order
                    trace_list[second][1] += s2[j2:j][::-1]
                    trace_list[second][2] += ' ' * diff_j

                # D[i,j] came from P and D
                elif D[i,j] == P[i,j] and D[i,j] == D[i-1,j-1] + fr.getSubsMatScore(s1[i-1], s2[j-1], subst_matrix_fn,beta):
                    trace_list.append(copy.deepcopy(trace_list[index])) # we need to split the trace_list sublist into 3 lists
                    indices_list.append(copy.deepcopy(indices_list[index])) #copy
                    # treat traceback[index] as the list where the traceback has come from P[i,j]
                    # i2 is the index after calling tracebackInP: after traversing through the
                    # P matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length, and closing the gap.
                    i2 = self.tracebackInP(P,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_i gives the number of the gap we have to insert in sequence 2, and in
                    # trace_list[index][1] . It is also used to decrement the corresponding indices_list row index.
                    diff_i = i - i2
                    indices_list[index][0] -= diff_i
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between i2 and i as we always add to
                    # trace_list in the reverse order
                    trace_list[index][0] += s1[i2:i][::-1]
                    trace_list[index][1] += '-' * diff_i
                    trace_list[index][2] += ' ' * diff_i

                    # treat traceback[second] as the list where the traceback has come from D[i-1,j-1] (decrease i and j)
                    second = len(trace_list) - 1
                    trace_list[second][0] += s1[i-1]
                    trace_list[second][1] += s2[j-1]
                    if s1[i-1] == s2[j-1]:
                        trace_list[second][2] += '*'
                    else:
                        trace_list[second][2] += ':'
                    indices_list[second][0] -= 1
                    indices_list[second][1] -= 1

                # D[i,j] came from D and Q.
                elif D[i,j] == Q[i,j] and D[i,j] == D[i-1,j-1] + fr.getSubsMatScore(s1[i-1], s2[j-1], subst_matrix_fn,beta):
                    trace_list.append(copy.deepcopy(trace_list[index])) # we need to split the trace_list sublist into 3 lists
                    indices_list.append(copy.deepcopy(indices_list[index])) #copy
                    # treat traceback[index] as the list where the traceback has come from P[i,j]
                    # j2 is the index after calling tracebackInQ: after traversing through the
                    # Q matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length and closing the gap.
                    j2 = self.tracebackInQ(Q,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_j gives the number of gaps we have to insert in sequence 1, and in
                    # trace_list[index][1] . It is also used to decrement the corresponding indices_list col index.
                    diff_j = j - j2
                    indices_list[index][1] -= diff_j
                    trace_list[index][0] += '-' * diff_j
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between j2 and j as we always add to
                    # trace_list in the reverse order
                    trace_list[index][1] += s2[j2:j][::-1]
                    trace_list[index][2] += ' ' * diff_j

                    # treat traceback[second] as the list where the traceback has come from D[i-1,j-1] (decrease i and j)
                    second = len(trace_list) - 1
                    trace_list[second][0] += s1[i-1]
                    trace_list[second][1] += s2[j-1]
                    if s1[i-1] == s2[j-1]:
                        trace_list[second][2] += '*'
                    else:
                        trace_list[second][2] += ':'
                    indices_list[second][0] -= 1
                    indices_list[second][1] -= 1

                # D[i,j] came from only P
                #indicates that a gap has been added in Sequence 1, so we have to decrement i

                elif D[i,j] == P[i,j]:
                    # i2 is the index after calling tracebackInP: after traversing through the
                    # P matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length, and closing the gap.
                    i2 = self.tracebackInP(P,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_i gives the number of the gap we have to insert in sequence 2, and in
                    # trace_list[index][1] . It is also used to decrement the corresponding indices_list row index.
                    diff_i = i - i2
                    indices_list[index][0] -= diff_i
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between i2 and i as we always add to
                    # trace_list in the reverse order
                    trace_list[index][0] += s1[i2:i][::-1]
                    trace_list[index][1] += '-' * diff_i
                    trace_list[index][2] += ' ' * diff_i

                # D[i,j] came from only Q.
                #indicates that a gap has been added in Sequence 2, so we have to decrement i
                elif D[i,j] == Q[i,j]:
                    # j2 is the index after calling tracebackInQ: after traversing through the
                    # Q matrix and returning to the D matrix. This corresponds to creating a gap, extending
                    # it to a certain length and closing the gap.
                    j2 = self.tracebackInQ(Q,D,indices_list[index][0],indices_list[index][1],beta,g)
                    # diff_j gives the number of gaps we have to insert in sequence 1, and in
                    # trace_list[index][1] . It is also used to decrement the corresponding indices_list col index.
                    diff_j = j - j2
                    indices_list[index][1] -= diff_j
                    trace_list[index][0] += '-' * diff_j
                    # Note that the indices for s2/s1 are one less than that for the D, P and Q matrices
                    # Also, while adding to the trace_list, we need to reverse the string between j2 and j as we always add to
                    # trace_list in the reverse order
                    trace_list[index][1] += s2[j2:j][::-1]
                    trace_list[index][2] += ' ' * diff_j

                # D[i,j] came from D[i-1,j-1]
                # Indicates a substitution
                elif D[i,j] == D[i-1,j-1] + fr.getSubsMatScore(s1[i-1], s2[j-1], subst_matrix_fn,beta):
                    trace_list[index][0] += s1[i-1]
                    trace_list[index][1] += s2[j-1]
                    if s1[i-1] == s2[j-1]:
                        trace_list[index][2] += '*'
                    else:
                        trace_list[index][2] += ':'
                    indices_list[index][0] -= 1
                    indices_list[index][1] -= 1

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
        print("\nSubstitution Matrix used: {}\n".format(subsMat))
        print("Number of optimal alignments between sequence 1 (with id {}) and sequence 2 (with id {}) is:{}".format(id1,id2,num_alignments))
        for index, alignment in enumerate(alignments):
            if complete_traceback == False:
                print ("One Random alignment (with alignment score",optimalScore,'):\n')
            else:
                print("Alignment",index+1,"(with alignment score",optimalScore,'):\n')

            # We need to print only 80 characters per line.
            # If number of characters in the ouput is exactly 80 (or a multiple of 80), the number of lines is the integer division
            # of length//linesize, for eg. for 80=> 80//80 = 1. If the number of characters in the ouptut is not exactly 80 or a
            # multiple of 80, number of lines is the integer division of length//linesize + 1. For eg, 159 => 159//80 + 1 = 2
            # This is implemented here using the ternary operator.
            lineSize = 80
            alignmentLength = len(alignment[0])
            numLines = ((len(alignment[0]) // lineSize) + 1) if alignmentLength % lineSize != 0 else (len(alignment[0]) // lineSize)
            for line in range(0,numLines):
                print(alignment[0][lineSize * line : lineSize * (line + 1)] + "\t<-  Sequence 1 ("+ id1+")")
                print(alignment[2][lineSize * line : lineSize * (line + 1)])
                print(alignment[1][lineSize * line : lineSize * (line + 1)] + "\t<- Sequence 2 ("+ id2+")")
                if line != numLines - 1:
                    print("Continued below... \n")

    def run(self,
            seq1_fasta_file,
            seq2_fasta_file,
            subst_matrix_fn,
            affine_cost_gap_open,
            affine_cost_gap_extend,
            complete_traceback):
            """This is the main function which parses fasta files,
            calls functions to create the D, P and Q
            matrices, and calls another function to print the final result"""

            fr = FriendClass()
            # Parse the fasta files. Get 2 sequences out of them
            # record1 and record2 are lists containing the sequences and ids in
            # fasta file 1 and 2 respectively.
            record1,record2 = fr.parseFastaFiles(seq1_fasta_file,seq2_fasta_file)
            # if there is a problem with the fasta file, list(SeqIO.parse) returns an empty list
            if len(record1) == 0 or len(record2) == 0:
                print("You have a problem with one of your FASTA files.  Hint: check if the first character is '>'")
                sys.exit(1) # error code 1
            id1 = record1[0].id
            s1 = str(record1[0].seq) # convert from Bio.Seq.Seq to str
            # Make sure s1 doesn't contain non-amino acid characters
            if fr.validateAminoSequence(s1) == 0:
                print("You have invalid character(s) in your 1st file")
                sys.exit(11)          # error code 11
            id2 = record2[0].id
            s2 = str(record2[0].seq) # convert from Bio.Seq.Seq to str
            #Make sure s2 doesn't contain any non-amino acid characters
            if fr.validateAminoSequence(s2) == 0:
                print("You have invalid character(s) in your 2nd file")
                sys.exit(12)  # error code 12

            # If gap costs are positive, take the additive inverse, return the
            # negative version of the same value.
            if affine_cost_gap_extend > 0:
                affine_cost_gap_extend = - affine_cost_gap_extend
                print ("Your gap extension cost is positive. I assume you want it to be negative, I have added a minus")
            if affine_cost_gap_open > 0:
                affine_cost_gap_open = - affine_cost_gap_open
                print ("Your gap opening cost is positive. I assume you want it to be negative, I have added a minus")

            # Affine gap function, g is actually g(1), which amounts to the cost of a new gap of 1.
            g = affine_cost_gap_open + affine_cost_gap_extend

            (D,P,Q,optimalScore) = self.buildMatrices(s1,s2,subst_matrix_fn,affine_cost_gap_open, affine_cost_gap_extend, g)
            alignment_strings = self.getAlignmentsFromTracebacks(s1,s2,subst_matrix_fn,D,P,Q,affine_cost_gap_open, affine_cost_gap_extend,g)
            num_alignments = len(alignment_strings)
            # If complete_traceback is False, assign only one alignment to the alignment_strings vector
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
    # run Gotoh with some parameters
    nw = Gotoh()
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
    parser.add_argument("gapExtendCost", type=int,
        help="Specify the cost of extending an existing gap")
    parser.add_argument("-c", action="store_true", default=False)
    parser.add_argument("-r", action="store_true", default=False)
    args=parser.parse_args()
    (id1, s1, id2, s2, optimalScore, alignment_strings, num_alignments) = nw.run(args.fasta1, args.fasta2,
     args.subsMatrixType, args.gapOpenCost, args.gapExtendCost, args.c)
    #nw.run('easy1.fasta','easy2.fasta','pam250',-11,-1,True)
