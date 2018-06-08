#-------------------------------------------------------------------------------
# Name:        Xpgma Concrete class
# Purpose:     This class has functions to a build UPGMA and WPGMA phylogenetic
#              trees from a FASTA file containing multiple sequences. It imports
#              the 'NeedlemanWunsch' class to get pairwise alignments. The
#              'FriendClass' class has some common functions to validate
#               sequences, get substitution matrix scores and parse FASTA files.
#
# Author:      Ashwath Sampath
#
# Created:     23-01-2018
# Copyright:   (c) Ashwath Sampath 2018
#-------------------------------------------------------------------------------

from prakt.xpgma import XpgmaBase
from needleman_wunsch import NeedlemanWunsch
from friend import FriendClass
from Bio import SeqIO
from Bio import Phylo
import numpy as np
import argparse
import random
import math
import sys
from io import StringIO
from decimal import Decimal
import matplotlib as plt


@XpgmaBase.register
class Xpgma(XpgmaBase):
    """This concrete class implements the abstract class XpgmaBase
    and contains functions to execute the UPGMA and WPGMA algorithms."""

    def similarityToDistance(self,s_ab,a,b,nw,alignment,subsMat,gapOpenCost):
        """ This function converts a similarity score to a distance score."""

        # 1. Calculate S(a,b)_rand using the formula on this page, but with linear gap costs:
        # http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Feng-Doolittle
        # Find length of the sequence
        L = len(alignment[0]) # same length for alignment[0],[1] and [2]
        # Find number of gaps in alignment[0] and alignment[2]
        N_g = alignment[0].count('-') + alignment[2].count('-')
        fr = FriendClass()
        sum_xy = 0
        # Randomize a and b to calculate s_rand.
        list_a = list(a)
        list_b = list(b)
        random.shuffle(list_a)
        random.shuffle(list_b)
        rand_a = "".join(list_a)
        rand_b = "".join(list_b)
        for i,x in enumerate(a):
            for j,y in enumerate(b):
                    s_xy = fr.getSubsMatScore(rand_a[i],rand_b[j],subsMat,gapOpenCost)
                    Na_x = a.count(x)
                    Nb_y = b.count(y)
                    sum_xy += (Na_x * Nb_y * s_xy)
        s_ab_rand = (sum_xy/L) + (N_g * gapOpenCost)

        # 2. Calculate s_ab_max
        (traceback_aa, s_aa) = nw.buildMatrices(a,a,subsMat,gapOpenCost)
        (traceback_bb, s_bb) = nw.buildMatrices(b,b,subsMat,gapOpenCost)
        s_ab_max = (s_aa + s_bb)/2
        #s_ab_eff is the normalized similarity: between 0 and 1.
        s_ab_eff = (s_ab - s_ab_rand) / (s_ab_max - s_ab_rand)
        d = - math.log(s_ab_eff)
        return d

    def pairwiseAlignment(self,s1,s2,subsMat,gapOpenCost):
        # Get the similarity using the Needleman Wunsch algorithm
        nw = NeedlemanWunsch()
        (traceback, optimalScore) = nw.buildMatrices(s1,s2,subsMat,gapOpenCost)
        alignment_strings = nw.getAlignmentsFromTracebacks(s1,s2,traceback)
        num_alignments = len(alignment_strings)
        # Get only one random optimal alignment
        randomNum = random.randint(0,num_alignments-1)
        alignment = alignment_strings[randomNum] # alignment is a list: ['','','']
        distance = self.similarityToDistance(optimalScore,s1,s2,nw,alignment,subsMat,gapOpenCost)
        return distance


    def getMin(self,distanceMatrix):
        """ This function gets the minimum value in a matrix, along with the index."""
        minimum = Decimal('inf')
        for i in range(len(distanceMatrix)):
            for j in range(len(distanceMatrix)):
                # We need to consider only the upper triangle matrix as it is symmetric.
                if i < j and distanceMatrix[i,j] < minimum:
                    # Round to 3 digits
                    minimum = Decimal(distanceMatrix[i,j]).quantize(Decimal('.010'))
                    row = i
                    col = j

        return row,col,minimum


    def updateMatrix(self,distanceMatrix,clustering,clust1,clust2,row,col):
        """This function updates the distance matrix for UPGMA and WPGMA"""

        # For WPGMA, the value in the denominator in the distanceMatrix update
        # step is just 2 and the values in the numerator are 1.
        if clustering == "WPGMA":
            clust1 = clust2 = 1
        for j in range(0,len(distanceMatrix)):
                if j != col and j != row:
                    # Update the columns in row no. 'row'
                    distanceMatrix[row,j] = ((clust1 * distanceMatrix[row,j]) +
                        (clust2 * distanceMatrix[col,j]))/ (clust1 + clust2)
                    # Update the rows in column no. 'row'. Matrix is symmetric.
                    distanceMatrix[j,row] = distanceMatrix[row,j]
        # Delete the row at index = col
        distanceMatrix = np.delete(distanceMatrix,col,axis=0)
        # Delete the column at index = col
        distanceMatrix = np.delete(distanceMatrix,col,axis=1)
        return distanceMatrix

    def updateMatrixAndClusters(self,c,c2,distanceMatrix,row,col,dmin,next_cluster,clustering,mapping,idNewick):
        """This function updates the distance matrix after each time 2 clusters are merged (by
        calling the updateMatrix function.
        It also updates the cluster list c. """
        # If clustering is UPGMA, we need to use the clustering lengths in the
        # distance matrix update calculation. We create clust1 and clust2, which
        # tell us how many values are in the c[row] and c[col] clusters, i.e.
        # the 2 clusters to be merged (row and col are the indices of the min value
        # in the distanceMatrix).
        clust1 = len(c[row].split(','))
        clust2 = len(c[col].split(','))

        # Update the UPGMA/WPGMA distance matrix
        distanceMatrix = self.updateMatrix(distanceMatrix, clustering, clust1,clust2,row,col)

        # Update the c list, which consists of all the current cluster names
        # The clusters are combined into the cluster at index 'row' and the
        # cluster at index 'col' is deleted.
        # To do this, we first calculate a string called 'added', which will be
        # assinged to c[row]. It includes the distance between the clusters: dmin/2
        # if 2 singleton clusters are joined, but dmin/2 - previous cluster distance
        # (stored in the dictionary 'mapping') if one or both clusters are not singleton
        # addedId is the same as added, except that it uses the original ids from
        # the fasta file (needed for later printing).

        if clust1 == 1 and clust2 == 1:
            added = '({}:{},{}:{})'.format(c[row],dmin/2,c[col],dmin/2)
            addedId =  '({}:{},{}:{})'.format(idNewick[row],dmin/2,idNewick[col],dmin/2)
        elif clust1 == 1 and clust2 != 1:
            added = '({}:{},{}:{})'.format(c[row],dmin/2,c[col],dmin/2 - mapping[c[col]])
            addedId = '({}:{},{}:{})'.format(idNewick[row],dmin/2,idNewick[col],dmin/2 - mapping[c[col]])
        elif clust1 != 1 and clust2 == 1:
            added = '({}:{},{}:{})'.format(c[row],dmin/2 - mapping[c[row]],c[col],dmin/2)
            addedId = '({}:{},{}:{})'.format(idNewick[row],dmin/2 - mapping[c[row]],idNewick[col],dmin/2)
        elif clust1 != 1 and clust2 != 1:
            added = '({}:{},{}:{})'.format(c[row],dmin/2 - mapping[c[row]],c[col],dmin/2 - mapping[c[col]])
            addedId = '({}:{},{}:{})'.format(idNewick[row],dmin/2 - mapping[c[row]],idNewick[col],dmin/2 - mapping[c[col]])

        # added2 is similar to added, but is much less complex as it doesn't take
        # distances into account
        added2 = '({},{})'.format(c2[row],c2[col])

        # The mapping dictionary is required to help each cluster calculation, keys
        # are sets of clusters and values are the min distances between them
        mapping[added] = dmin/2
        c[row] = added
        del c[col]

        c2[row] = added2
        del c2[col]

        idNewick[row] = addedId
        del idNewick[col]
        # Give a name to the new cluster=> C+next_cluster
        new_cluster = 'C{}'.format(str(next_cluster))
        return c,c2,distanceMatrix,idNewick

    def UandWpgma(self,ids,s,n,subsMat,gapOpenCost,clustering):
        """ This function builds a UPGMA or a WPGMA phylogenetic tree. The same
        function is used for creating both types of trees as there is only a
        difference in the matrix update formulae. It takes as input the
        lists of ids and sequences from the fasta file along with the number of
        sequences 'n' (ids and s are lists), the substitution matrix, gap open cost
        and the type of tree. It returns the tree in Newick format."""

        # c is a list which holds all the clusters. It will also eventually
        # hold the distance between cluster.
        c = []
        # c2 is identical to c, but will NOT hold distances.
        c2 = []
        # Initialize c and c2 to have n singleton clusters 0,...,n-1. c will be used
        # to get the output in Newick format with distances, c2 will be used to
        # get the output in Newick format without distances.

        # idNewick contains the actual IDs in the fasta file
        idNewick = []
        # This dict will map from the generic node name (C+node number), to the
        # actual Id in the file.
        cToIdDict = dict()
        for index, id in enumerate(ids):
            cToIdDict['C'+str(index)] = id

        for i in range(0,n):
            c.append('C'+str(i))
            c2.append('C'+str(i))
            idNewick.append(cToIdDict.get(c[i]))

        # create distanceMatrix as a 2D numpy array, distance[i,i] is always zero.
        distanceMatrix = np.zeros((n,n))
        # next_cluster just contains the name of the next new cluster to be created.
        next_cluster = n
        for i in range(0,n-1):
            for j in range(i+1,n):
                # call the pairwiseAlignment function to calculate each cell
                # of the distance matrix. This internally calls the Needleman-Wunsch
                # function, and converts similarity scores into distance scores
                # using the Feng-Doolittle method.
                distanceMatrix[i,j] = self.pairwiseAlignment(s[i],s[j],subsMat,gapOpenCost)
                # Matrix is symmetric, so the following statement is valid
                distanceMatrix[j,i] = distanceMatrix[i,j]

        # Create an empty dictionary, which will be sent to the updateMatrixAndClusters
        # function
        mapping = dict()
        # Loop and combine clusters until there is only one cluster left
        while len(c) > 1:
            # stop when there is only 1 cluster left
            # Get the minimum element of the matrix with the corresponding indices
            row,col,dmin = self.getMin(distanceMatrix)
            # Round the matrix to 3 digits so that it is easier to display
            distanceMatrix = np.round(distanceMatrix,3)
            # Call the function to update the distance matrix and the cluster list c
            # (clusters and rows/columns in the matrix are combined in this function)
            # We also get back the output in almost-Newick format from the function
            # in lists c and c2.
            c,c2,distanceMatrix,idNewick = self.updateMatrixAndClusters(c, c2,
                distanceMatrix,row,col,dmin,next_cluster,clustering,mapping,idNewick)
            # Update the next cluster number
            next_cluster += 1

        # Add a semicolon as the Newick format requires it. As it gets added as
        # the second member of the c list, we need to use join to convert it into
        # a string 'newick', which is returned. The clusters without distances are
        # also returned, as are the clusters with the original ids from the file.
        c += ';'
        c2 += ';'
        idNewick += ';'
        newick = ''.join([cluster for cluster in c])
        newickNoDistance = ''.join([cluster for cluster in c2])
        newickIds = ''.join([cluster for cluster in idNewick])
        return(newick, newickNoDistance, distanceMatrix,newickIds)

    def printTree(self,newick,clustering):
        """This function takes a string in Newick format and prints it in the
        format of a Dendogram/Phylogenetic tree. For this, it uses the Phylo
        module of Biopython"""
        print("The {} tree in Newick format is {}".format(clustering, newick))
        print("The {} dendogram of the tree is displayed in the figure:".format(clustering))

        # String needs to be parsed using the io module's StringIO function
        handle = StringIO(newick)
        tree = Phylo.read(handle, "newick")
        Phylo.draw(tree, branch_labels=lambda c: c.branch_length)


    def run(self,
            seq_fasta_file,
            subst_matrix_fn,
            cost_gap_open,
            clustering):
            """This is the main function which parses the fasta file,
            calls functions to create the UPGMA and WPGMA trees, and
            calls another function to print the final result"""

            fr = FriendClass()
            # Parse the fasta file. Get 2 sequences out of them

            # record  is a list containing the sequences and ids in
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

            # Call the function which will call other functions to create a UPGMA/WPGMA
            # tree based on the 'clustering' argument that is sent.
            newick, newickNoDistance, distanceMatrix, newickIds = self.UandWpgma(ids,s,num_sequences,subst_matrix_fn,cost_gap_open,clustering)

            # Call a function which prints the 3 strings on the console
            self.printTree(newickIds,clustering)
            # Return both the Newick output with distances and the output with just
            # the amino acid names.
            return newick, newickNoDistance, newickIds

if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    gma = Xpgma()
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
    newick, newickNoDistance, newickIds = gma.run(args.fasta, args.subsMatrixType, args.gapOpenCost, args.clustering.upper())
