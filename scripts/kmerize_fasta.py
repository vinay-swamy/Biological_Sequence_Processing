
#%%
from Bio import SeqIO
from os import path 
from pickle import dump
import sys 



def kmerize_sequence(seq,k, n, label):
    '''
    k is the length of the kmer, and n is the distance between each kmer
    '''
    res = [seq[i:(i+k)] for i in range(0, len(seq)-k+1, n)]
    
    return((res, label))


def process_fasta(file, k, n, label):
    fpath = path.expanduser(file) if '~' in file else path.abspath(file)
    fasta = SeqIO.parse(fpath, format='fasta')
    kmer_list = [kmerize_sequence(
        str(record.seq), k, n,label) for record in fasta]
    return(kmer_list)





#%%

ref_fasta_file=sys.argv[1]
dummy_fasta_file=sys.argv[2]
kmer_size=int(sys.argv[3])# k
kmer_dist=int(sys.argv[4])# l
kmer_outfile=sys.argv[5]
lineSentence_outfile=sys.argv[6]

all_records_kmers = process_fasta(ref_fasta_file, kmer_size, kmer_dist, 1) + process_fasta(dummy_fasta_file, kmer_size, kmer_dist, 0)


# %%

with open(kmer_outfile, 'wb+') as ofl:
    dump(all_records_kmers, ofl)

with open(lineSentence_outfile, 'w+') as lso:
    [lso.write(' '.join(record[0]) + '\n') for record in all_records_kmers]

