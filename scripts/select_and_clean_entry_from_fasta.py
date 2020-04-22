#!/use/bin/env python3qq
import argparse
parser = argparse.ArgumentParser(description='Remove sequences from a fasta')
parser.add_argument('--infasta', type=str,help='fasta to remove from ')

parser.add_argument('--txToKeep', type=str, help='entries to keep')

parser.add_argument('--outfasta', type=str,help='file for filtered fasta ')
args=parser.parse_args()

with open(args.infasta) as infasta, open(args.txToKeep) as tx, open(args.outfasta,'w+') as outfasta:
    names=set()
    
    for line in tx:
        names.add(line.strip())
    oldline=infasta.readline().strip()#.split('.')[0]
    id_1 = oldline.split('|')[1]
    id_2 = oldline.split('|')[2].split(' ')[0]
    while oldline:
        found = id_1 in names or id_2 in names
        if not found and  oldline[0] == '>':
            write=False
        elif found and oldline[0] == '>':
            which_found = id_1 if id_1 in names else id_2
            write=True
            oldline = '>' + which_found
        if write:
            outfasta.write(oldline + '\n')
        oldline=infasta.readline().strip()#.split('.')[0]
        if  oldline  and  oldline[0]  == '>':
            id_1 = oldline.split('|')[1]
            id_2 = oldline.split('|')[2].split(' ')[0]
