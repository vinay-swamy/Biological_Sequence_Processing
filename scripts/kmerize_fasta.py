
#%%
from Bio import SeqIO
from os import path 
from pickle import dump
import sys 
def kmerize_sequence(seq,k, n=1):
    '''
    k is the length of the kmer, and n is the distance between each kmer
    '''
    res = [seq[i:(i+k)] for i in range(0, len(seq)-k+1, n)]
    return(res)

#%%

fasta=sys.argv[1]
kmer_size=int(sys.argv[2])
kmer_dist=int(sys.argv[3])
kmer_outfile=sys.argv[4]

fasta_path = path.expanduser(fasta) if '~' in fasta else path.abspath(fasta)
fasta=SeqIO.parse(fasta_path, format='fasta')
kmer_list=[kmerize_sequence(str(record.seq),kmer_size,kmer_dist) for record in fasta ]
# %%

with open(kmer_outfile, 'wb+') as ofl:
    dump(kmer_list, ofl)
# %%
