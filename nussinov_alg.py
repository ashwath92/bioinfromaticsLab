#-------------------------------------------------------------------------------
# Name:         Nussinov Concrete class
# Purpose:      This class has functions to implement the Nussinov algorithm
#               for maximizing the number of base pairs in an RNA structure.
#               It takes as input a FASTA file containing an RNA sequence, and
#               returns the base pairs in the dot-bracket notation.
#               It imports the 'FriendClass' class, has some common functions to
#               validate sequences, get substitution matrix scores and parse FASTA files.
#
# Author:      Ashwath Sampath
#
# Created:     16-1-2018
# Copyright:   (c) Ashwath Sampath 2018
#-------------------------------------------------------------------------------

from prakt.nussinov import NussinovBase
from friend import FriendClass
from Bio import SeqIO
import numpy as np
import argparse
import random
import copy
import sys

@NussinovBase.register
class Nussinov(NussinovBase):
    """This concrete class implements the abstract class NussinovBase
    and contains functions to execute the whole Nussinov algorithm."""

    def __init__(self):
        """Constructor"""
        self.indices_list= []
        # self.indices_list = [[]]
        self.trace_list = ["",""]


    def is_complementary(self,nucl1,nucl2):
        if (nucl1 == 'G' and nucl2 == 'C') or nucl1 == 'C' and nucl2 == 'G':
            return True
        elif (nucl1 == 'A' and nucl2 == 'U') or (nucl1 == 'U' and nucl2 == 'A'):
            return True
        # Weaker bonds between G and U are also formed in RNA
        #elif (nucl1 == 'G' and nucl2 == 'U') or (nucl1 == 'U' and nucl2 == 'G'):
            #return True
        else:
            return False


    def build_matrix(self,seq):
        """This function constructs the Nussinov matrix from the sequence passed
        in the arguments."""
        n = len(seq)
        # Create a Nussinov matrix where there are n rows and n+1 columns. The extra first column
        # is needed for calculation. This automatically initializes the principal diagonal elements
        # and all elements to the left of the principal diagonal to 0. N is the Nussinov matrix.
        N = np.zeros((n,n+1),dtype=int)
        # Start a loop on variable d. d is just the difference between j and i. This loop is needed
        # because the order of calculation is diagonal-wise (diagonals with d=j-i=2, d=3 and so on).
        # Max d = max difference between i and j: n (which corresponds to i=first char, j = last char)
        for d in range(2,n+1):
            i = 0 # i is reset to 0 after finishing each diagonal
            # E.g. If n=5, i goes upto 3 (n-d) in first iteration of the d for loop,
            # upto 2 in the second iteration, upto 1 in the third iteration, only
            # i=0 in the 4th and final iteration.
            while i < (n-d+1):
                j = i + d
                # k starts at i + 1 because of extra column at the beginning
                k = i + 1
                # Create a list of paired values: different values corresponding
                # to different i,j,k values will be present.
                paired = []
                unpaired = N[i,j-1]
                maxPaired = 0
                while k < j:
                    if self.is_complementary(seq[k-1], seq[j-1]):
                        # Formula below suitably adjusted because k starts at i+1
                        paired.append(N[i,k-1] + N[k,j-1] + 1)
                    k += 1
                # default maxPaired to 0 when there are no complementary bases,
                # we don't want a None.
                maxPaired = max(paired,default=0)
                N[i,j] = max(maxPaired,unpaired)
                i += 1
        #np.savetxt('out.txt',N,delimiter=',')

        # Return the Nussiov matrix
        return N

    def tracebackInN(self,N,seq,i,j):
        """ This recursive function gets the traceback from the Nussinov matrix
        that was created in the build_matrix function. The traceback starts
        from the N(1,n) term. This is a direct implementation of the algorithm
        given in the Uni-Freiburg slides. The 'trace_list' and 'indices_list'
        class variables are updated. These give the paired nucleotides and
        their indices respectively"""

        k = i + 1
        # Stopping condition of the algorithm: this indicates 2 characters cannot
        # be paired
        if j <= i:
            return
        # No pairing, go to the previous character in the sequence and try again
        elif N[i,j] == N[i,j-1]:
            return self.tracebackInN(N,seq,i,j-1)
        # seq[k-1] and seq[j-1] get paired if they are complementary
        else:
            # Loop through k, the second condition checks if seq[k-1] is
            # complementary to seq[j-1].

            while k < j:
                #for index,[k,j] in enumerate(self.indices_list):
                if self.is_complementary(seq[j-1],seq[k-1]):
                    if N[i,j] == N[i,k-1] + N[k,j-1] + 1:
                        # Corresponding indices and characters in [0] and [1] of
                        # indices_list and trace_list are paired.
                        #self.indices_list[self.index].append([k-1,j-1])
                        self.indices_list.append([k-1,j-1])
                        self.trace_list[0] += seq[k-1]
                        self.trace_list[1] += seq[j-1]
                        # Divide and conquer. 2 recursive calls: 1 for the 1st part
                        # of the sequence, 1 for the 2nd part.
                        # The indices simply correspond to the if condition above.
                        self.tracebackInN(N,seq,i,k-1)
                        self.tracebackInN(N,seq,k,j-1)
                        return
                        # Multiple tracebacks
                        #self.indices_list.append(copy.deepcopy(self.indices_list[self.index]))
                        #self.index +=1
                k += 1


    def printer(self,seq,id1):
        length = len(seq)
        # crete a list of '.'s as '.' refers to a non-paired nucleotide in the
        # dot-bracket notation. We need a list as strings are immutable.
        output = ['.' for i in range(length)]
        # indices_list contains all the nucleotide pairings. It has 2 sub-lists
        # which have the first element to pair  and the second element to pair resp.
        for element in self.indices_list:
            output[element[0]] = '('
            output[element[1]] = ')'
        # Convert the output into a string for display
        outputString = ''.join(output)
        # Count the number of ( in the output to get the optimal number of base pairs.
        # Note that we could've also counted the number of ) as they have to be equal.
        optimalNumberBp =  outputString.count('(')
        print("ID = ",id1)
        print("The optimal number of base pairs is: ", optimalNumberBp)
        print("The optimal base pair alignmnent in dot-bracket notation is:")
        # Convert the output list to a string
        print(outputString)
        print(seq)
        # Return the dot-bracket notation.
        return outputString

    def run(self,
            seq_fasta_fn):
        """
        Fold RNA with Nussinov algorithm.

        Args:
            seq_fasta_fn: path to fasta file containing sequence

        Returns:
            tuple of
            (id_seq: fasta id of sequence,
             seq: sequence,
             structure: dot-bracket string of optimal folding)
        """

        """This is the main function which parses fasta files,
            calls functions to create Needleman Wunsch and traceback
            matrices, and calls another function to print the final result"""

        fr = FriendClass()
        # Parse the fasta files. Get 2 sequences out of them
        # record is a list containing the sequences and ids in
        # the fasta file.
        record = fr.parseMultSequenceFastaFile(seq_fasta_fn)
        # if there is a problem with the fasta files, list(SeqIO.parse)
        # returns an empty list
        if len(record) == 0:
            print("You have a problem with your FASTA file.  H"
              "int: check if the first character is '>'")
            sys.exit(1) # error code 1
        id1 = record[0].id
        s1 = str(record[0].seq) # convert from Bio.Seq.Seq to str
        # Make sure s1 doesn't contain non-RNA characters
        if fr.validateRNASequence(s1) == 0:
            print("You have invalid character(s) in your file")
            sys.exit(11)         # error code 11

        # Build the Nussinov matrix
        N = self.build_matrix(s1)
        i = 0
        j = len(s1)
        # Get the traceback: it is stored in the class-variables trace_list and
        # indices_list
        self.tracebackInN(N,s1,i,j)
        # The
        dot_bracket = self.printer(s1,id1)
        return N, dot_bracket

if __name__ == '__main__':
    # run Nussinov with some parameters
    nuss = Nussinov()
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help='Specify the path to the sequence file')
    args=parser.parse_args()
    N,dot_bracket = nuss.run(args.fasta)
