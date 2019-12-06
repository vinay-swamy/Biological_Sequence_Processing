#%%
from Bio import SeqIO
import sys 

#%%
infasta=sys.argv[1]
MIN_LENGTH=int(sys.argv[2])
outfasta=sys.argv[3]



#%%
fasta = SeqIO.parse(infasta, 'fasta')
# there are some rogue N's around, want to try and catch them 
good_tx = [record for record in fasta if 'N' not in record and len(record) >= MIN_LENGTH]



# %%
SeqIO.write(good_tx, outfasta, 'fasta')

