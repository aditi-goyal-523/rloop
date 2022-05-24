import argparse
import numpy as np

import seqlib

# Functions
def kmer_freq(seqs, k):
    k_count = {}
    total = 0
    for seq in seqs:
        seq = seq.upper()
        for i in range(len(seq) - k +1):
            kmer = seq[i:i+k]
            if "N" not in kmer:
                if kmer not in k_count: 
                    k_count[kmer]=0
                k_count[kmer]+=1
                total +=1
    
    freq = {}
	
    for kmer in k_count:
        freq[kmer] = k_count[kmer] / total
	
    return freq

def avg_length(seqs):
    length= 0
    counter= 0
   
    for seq in seqs:
        length+=len(seq)
        counter+=1
    
    return(length/counter)



#----------------------------------------------
#Argparse


parser = argparse.ArgumentParser(description='Predict R loop presence using HMM')
parser.add_argument('--plus', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--minus', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--nor', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--test', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--k', required=False, type=int, default=5,
 metavar="<int>", help='length of k-mer [%(default)i]')
# parser.add_argument('--x', required=False, type=int, default=4,
#  metavar='<int>', help='x-fold cross-validation [%(default)s]')
# parser.add_argument('--seed', required=False, type=int,
#  metavar='<int>', help='random seed')
arg = parser.parse_args()
#-----------------------------------------------
#################   setup   ####################

plus_fasta=seqlib.read_fasta(arg.plus)
minus_fasta=seqlib.read_fasta(arg.minus)
nor_fasta=seqlib.read_fasta(arg.nor)

k=arg.k

plus=[]
for name, seq in plus_fasta:
    plus.append(seq)

minus=[]
for name, seq in minus_fasta:
    minus.append(seq)

nor=[]
for name, seq in nor_fasta:
    nor.append(seq)

plus_freq=kmer_freq(plus, k)
minus_freq=kmer_freq(minus, k)
nor_freq=kmer_freq(nor, k)

#emission probabilities
plus_leave=1/avg_length(plus)
plus_stay=1-plus_leave

minus_leave=1/avg_length(minus)
minus_stay=1-minus_leave

nor_leave=(plus_leave+minus_leave)/2
nor_stay=(plus_stay+minus_stay)/2



#-----------------------------------------------
#############   Intialization   ################

test_fasta=seqlib.read_fasta(arg.test)

test=[]
for name, seq in test_fasta:
    test.append(seq)

test=''.join(test)
test.upper()

kmers=[]
for i in range(len(test) - k +1):
    kmer = test[i:i+k]
    if "N" not in kmer:
        kmers.append(kmer)

l=len(kmers)+1
plus_row=[None]*l
minus_row=[None]*l
nor_row=[None]*l

matrix=[plus_row, minus_row, nor_row]

matrix[0][0] = ("N", 0.05) #R loops occur across approximatley 5% of the genome
matrix[1][0] = ("N", 0.05) #R loops occur across approximatley 5% of the genome
matrix[2][0] = ("N", 0.90) #R loops occur across approximatley 5% of the genome

#-----------------------------------------------
#############   Scoring   ################

def plus_score(kmer, m, i):

    #make a dictionary of results
    result={'p2p':-1, 'p2m':-1, 'p2n':-1}
    
    #store the p for 
    result['p2p']=np.log(plus_freq[kmer]) + np.log(plus_stay) + np.log(m[0][i - 1][1])
    result['p2m']=np.log(minus_freq[kmer]) + np.log(plus_leave) + np.log(m[0][i - 1][1])
    result['p2n']=np.log(nor_freq[kmer]) + np.log(plus_leave) + np.log(m[0][i - 1][1])

    max_val=max(result.values())
    max_key=max(result, key=result.get)

    return((max_key, np.exp(max_val)))

def minus_score(kmer, m, i):

    #make a dictionary of results
    result={'m2m':-1, 'm2p':-1, 'm2n':-1}
    
    #store the p for 
    result['m2m']=np.log(minus_freq[kmer]) + np.log(minus_stay) + np.log(m[1][i - 1][1])
    result['m2p']=np.log(plus_freq[kmer]) + np.log(minus_leave) + np.log(m[1][i - 1][1])
    result['m2n']=np.log(nor_freq[kmer]) + np.log(minus_leave) + np.log(m[1][i - 1][1])

    max_val=max(result.values())
    max_key=max(result, key=result.get)

    return((max_key, np.exp(max_val)))

def nor_score(kmer, m, i):

    #make a dictionary of results
    result={'n2n':-1, 'n2p':-1, 'n2m':-1}
    
    #store the p for 
    result['n2n']=np.log(nor_freq[kmer]) + np.log(nor_stay) + np.log(m[2][i - 1][1])
    result['n2p']=np.log(plus_freq[kmer]) + np.log(nor_leave) + np.log(m[2][i - 1][1])
    result['n2m']=np.log(minus_freq[kmer]) + np.log(nor_leave) + np.log(m[2][i - 1][1])

    max_val=max(result.values())
    max_key=max(result, key=result.get)

    return((max_key, np.exp(max_val)))

#-----------------------------------------------
#############   Filling the Matrix   ################

for i in range(1, len(kmers)+1):
    section=kmers[i-1]
    matrix[0][i] = plus_score(section, matrix, i)
    matrix[1][i] = minus_score(section, matrix, i)
    matrix[2][i] = nor_score(section, matrix, i)

#-----------------------------------------------
#############   Traceback   ################

def trace_to_dir(tup):
    if '2p' in tup[0]: #aka the kmer ends in a plus R loop state
        return("p")
    elif '2m' in tup[0]:
        return("m")
    else:
        return("n")

traceback=[None]*len(kmers)

j=len(kmers)-1
while j>=0:
    plus_tup=matrix[0][j]
    minus_tup=matrix[1][j]
    nor_tup=matrix[2][j]

    max_tup=max(plus_tup, minus_tup, nor_tup)
    traceback[j]=trace_to_dir(max_tup)
    j-=1




    


    


