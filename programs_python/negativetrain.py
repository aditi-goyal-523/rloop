from more_itertools import random_permutation
import sys
import seqlib

#----------------------------------------------------
def complement(seq):
    reverse = seq[::-1]

    newseq = ""

    for nt in reverse:
        if nt == "A":
            newseq += "T"
        elif nt == "T":
            newseq += "A"
        elif nt == "C":
            newseq += "G"
        elif nt == "G":
            newseq += "C"

    return(newseq)


def randomizer(seq):
    return(''.join(random_permutation(seq)))
#----------------------------------------------------------------------------------

#load in data
filename=sys.argv[1]
myfasta = seqlib.read_fasta(filename)

#set empty holder variables
random_seqs=[]
complement_seqs=[]

#get random and complement for all s in myfasta
for name, s in myfasta:
    random_seqs.append(randomizer(s))
    complement_seqs.append(complement(s))

#set output name
parts=filename.split(".")
savename=parts[0]+"."+parts[1]+".negativecontrol"

#get length of both vectors
rs_length=len(random_seqs)
cs_length=len(complement_seqs)

#write output
with open(savename, 'w') as writer:
    #for each seq in randomized sequences
    for i in range(rs_length):
        #write fasta header
        writer.write(">seq %d of %d | randomized control | seq length: %d\n" % (i+1, rs_length, len(random_seqs[i])))
        #write the sequence
        writer.write("%s\n" % random_seqs[i])

    #for each seq in complement sequences
    for i  in range(cs_length):
        #write fasta header
        writer.write(">seq %d of %d | complement control | seq length: %d\n" % (i+1, cs_length, len(complement_seqs[i])))
        #write the sequence
        writer.write("%s\n" % complement_seqs[i])
writer.close()