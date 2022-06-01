#step 0: training


#take in a file of R loop 
#take in file of control
#get kmer frequency for each
#store kmer frequency in a dictionary

#calculating transition probability:
#- take avg length of an r loop sequence
#- p of going from R loop to non-rloop: 1/length
#- p of going staying in R loop is 1-(1/length)

'''
***VITERBI TO RLOOP: instead of an ACGT frequency table, which shows the probabilities of each nucleotide being present, 
you have a table of all possible kmers, and the frequency at which they occur. ****
'''

#step 1: Initialization

need a matrix for the transition probabilities
need a matrix for the traceback itself

'''
first col represents what we expect to be the frequency of r loop/nonrloop 
literature indicates R loops occur across approx 5% of the genome
so Rloop = .05, nonrloop=.95
'''


'''

'''
           kmer 1 | kmer 2 | .... | kmer n
--------------------------------------------
| R loop | p that the kmer is an R loop region
--------------------------------------------
| Non R  | p that the kmer is not in R loop region

'''

#step 2: Fill
'''
for each fill, you need 3 components:
- p of the cell being in the n-1 state
- p of whatever state n is
- the transition probability from n-1 to n

calculate these probabilities for each possible route
store the maximum one in your matrix
'''

#step 3: traceback
'''
once the table is filled, start at the end and pick the max value of the column
move to the left, picking the max value of each column every time
store the path in a final variable
note the locations of ROW change and what direction it changes in
these represent switches from r loop and back
'''







