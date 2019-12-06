#%% imports
import pandas as pd
import sys
MIN_TX_LENGTH=100

tx_length_file=sys.argv[1]
MAX_DUMMY_TX=int(sys.argv[2])
outfile=sys.argv[3]

tx_lengths = (pd.read_csv(tx_length_file,
                    sep='\t', names=['seqid', 'start', 'end'])
    .assign(length= lambda x: x['end'] - x['start'])['length'] )
lengths_to_shuffle = tx_lengths.sample(n=MAX_DUMMY_TX, replace=True, random_state=42)
bed_to_shuffle=pd.DataFrame(data={'seqid':'chrZ', 'start':0, 'end':lengths_to_shuffle})
bed_to_shuffle.to_csv(outfile, sep='\t', header=False, index=False)



