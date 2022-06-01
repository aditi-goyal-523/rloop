import seqlib
import argparse
import random

#Functions
def kmer_freq(seqs, k):
    k_count = {}
    total = 0
    #seq = seq.upper()
    for seq in seqs:
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


def weight_by_freq(freq, seq, k):
    score = 1
    for i in range(len(seq) -k +1):
        kmer = seq[i:i+k]
        score *= freq[kmer]
    return score

def cross_validation(seqs, n):
    length = len(seqs)
    sector_size = round(length/n)
    set=[]

    #break the sequence up into sectors of size n each
    #n determined by user, levels of cross validation
    for i in range(n):
        end=(i+1)*sector_size
        sector = seqs[i:end]
        set.append(sector)

   #for the set of sectors, return one of them as the training, and rest as test sets     
    for i in range(len(set)):
        train=[]
        test=[]

        for j in range(len(set)):
            if i !=j:
                train.extend(set[j])
            else: 
                test=set[j]
        yield(train, test)
        


#----------------------------------------------
#Argparse
parser = argparse.ArgumentParser(description='Predict R loop presence using kmer analysis')
parser.add_argument('--pos', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--neg', required=True, type=str, metavar='<path>',
help='path to a fasta file, may be compressed')
parser.add_argument('--k', required=False, type=int, default=5,
 metavar="<int>", help='length of k-mer [%(default)i]')
parser.add_argument('--x', required=False, type=int, default=4,
 metavar='<int>', help='x-fold cross-validation [%(default)s]')
# parser.add_argument('--seed', required=False, type=int,
#  metavar='<int>', help='random seed')
arg = parser.parse_args()
#-----------------------------------------------
#Data Loading

seqs1 = [(1, seq) for name, seq in seqlib.read_fasta(arg.pos)] #each name becomes a "1"
seqs0 = [(0, seq) for name, seq in seqlib.read_fasta(arg.neg)]

seqs = seqs1 + seqs0 #get full set of sequences
	
random.shuffle(seqs)

#------------------------------------------------
#Cross Validation

accuracies=[]

pos_set_train=[]
neg_set_train=[]

for train, test in cross_validation(seqs, arg.x):
    #now there is a training set and testing set
    #each one needs to be divided into what they actually are
   
    for cat, seq in train:
        
        if cat == 1:
            pos_set_train.append(seq)
        else:
            neg_set_train.append(seq)

    pos_freq = kmer_freq(pos_set_train, arg.k) #gets frequencies of each kmer in positive set
    neg_freq = kmer_freq(neg_set_train, arg.k) #gets frequencies of each kmer in negative set

    true_pos=0
    true_neg=0
    false_pos=0
    false_neg=0

    for pair in test:
        cat, seq = pair

        #run this sequence using both frequencies, and determine which is the better model
        pos_score=weight_by_freq(pos_freq, seq, arg.k)
        neg_score=weight_by_freq(neg_freq, seq, arg.k)

        if cat == 1: #actually a positive
            if pos_score > neg_score: true_pos += 1 #positive model worked better
            else:                     false_neg += 1 #positive misclassified as negative
        else: #actually a negative
            if neg_score > pos_score: true_neg +=1 #negative model worked better
            else:                     false_pos += 1 #negative misclassified as positive

    accuracy = (true_pos + true_neg) / (true_pos + true_neg + false_pos + false_neg) 
    accuracies.append(accuracy)
	
    print(f"True Positives: {true_pos}\n True Negatives: {true_neg} \n False Positives {false_pos}\n False Negatives: {false_neg} \n Accuracy: {accuracy}")


