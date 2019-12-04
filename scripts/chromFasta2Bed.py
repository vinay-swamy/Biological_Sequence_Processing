import sys
import gzip
with gzip.open(sys.argv[1], 'rt') as infasta, open(sys.argv[2], 'w+') as outfile:
    line=infasta.readline().strip('\n')
    header=line
    length=0
    First=False
    while line:
        if '>' in line and First :
            outfile.write( '\t'.join([header[1:].split(' ')[0],'0' ,str(length)]) +'\n')
            header=line
            length=0
        else:
            First=True
            length+=len(line.strip())
        line=infasta.readline().strip('\n')