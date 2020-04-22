#!/bin/bash 

wget -O ref/idmapping.dat.gz "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
mkdir -p new_fasta
cd new_fasta
wget -r "ftp://ftp.ebi.ac.uk/pub/databases/ena/coding/release/std/fasta/"
mv ftp.ebi.ac.uk/pub/databases/ena/coding/release/std/fasta/rel_* .
rm ftp.ebi.ac.uk/ -rf 
for i in * ; do gunzip $i ; done 
cat *.fasta >  ../ref/all_ena_protein_coding.fasta
