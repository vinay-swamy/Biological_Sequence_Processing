#%% imports
import pandas as pd
import pybedtools as pbt
import sys
MIN_TX_LENGTH=100
MAX_DUMMY_TX=200000
tx_length_file=sys.argv[1]
outfile=sys.argv[2]

#%%


tx_lengths = (pd.read_csv(tx_length_file,
                    sep='\t', names=['seqid', 'start', 'end'])
    .assign(length= lambda x: x['end'] - x['start'])['length'] )


lengths_to_shuffle=set(tx_lengths)
for i in range(1,11):
    col=pd.Series().append(tx_lengths + i).append(tx_lengths -i)
    lengths_to_shuffle.update(set(col))
l2s=pd.Series(list(lengths_to_shuffle))
l2s=l2s[l2s >= MIN_TX_LENGTH]
lengths_to_shuffle=l2s.sample(n=MAX_DUMMY_TX, replace=True, random_state=42)



#%%
bed_to_shuffle=pd.DataFrame(data={'seqid':'chrZ', 'start':0, 'end':lengths_to_shuffle})
bed_to_shuffle.to_csv(outfile, sep='\t', header=False, index=False)



