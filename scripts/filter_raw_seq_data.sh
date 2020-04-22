#!/bin/bash 

zgrep "EMBL" ref/idmapping.dat.gz  > ref/uniprot2allembl.tsv
sed 's/^/\^/' ref/uniprot_target_protein_uniprot_ids.txt > ref/uniprotids_re.txt
grep -Ef ref/uniprotids_re.txt ref/uniprot2allembl.tsv > ref/target_uniprot2target_embl.txt
cat ref/target_uniprot2target_embl.txt | cut -f3 >ref/target_embl_ids.txt 
python3 scripts/select_and_clean_entry_from_fasta.py --infasta ref/all_ena_protein_coding.fasta --txToKeep ref/target_embl_ids.txt --outfasta ref/target_protein_seqs.fa
grep -o -Ff ref/target_embl_ids.txt ref/target_protein_seqs.fa | grep -Ff - ref/target_uniprot2target_embl.txt > ref/target_uniprot_found_in_target_embl.txt
python3 scripts/clean_labels.py
cat data/annotation/label_mappings.csv | cut -d',' -f1 > data/annotation/final_embl_ids.txt 
python3 ~/scripts/select_entry_from_fasta.py --infasta ref/target_protein_seqs.fa --txToKeep  data/annotation/final_embl_ids.txt  --outfasta data/seqs/complete_target_protein_nt_seqs.fa