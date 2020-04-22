########## 
working_dir = config['working_dir']

ref_GTF=config['ref_GTF']
ref_genome=config['ref_genome']
bam_path = config['bam_path']
STARindex = bam_path + 'ref/STARindex'

########## 



rule all:
    input:
        expand('data/model_results/uniprot_rf_dm-{dm}_wc-{mincount}_kmer-{size}_dims-{dims}/model_results.csv',
                    dm=['0', '1'],
                    mincount=['3', '15', '30'],
                    size=[str(i) for i in ['8', '10', '12']],
                    dims=['100',  '300', '600', '1000']
        )


  
'''
NOTE: Need to first run scripts/get_raw_seq_data.sh on helix, and downlod uniprot_annotation from webportal, and upload back to biowulf before starting 

'''
rule prep_data:
    output:
        seqs = 'data/seqs/complete_target_protein_nt_seqs.fa',
        lab_df = 'data/annotation/label_mappings.csv', 
        all_anno  = 'data/annotation/all_seqmeta_info.csv'
    shell:
        '''
        zgrep "EMBL" ref/idmapping.dat.gz  > ref/uniprot2allembl.tsv
        sed 's/^/\^/' ref/uniprot_target_protein_uniprot_ids.txt > ref/uniprotids_re.txt
        grep -Ef ref/uniprotids_re.txt ref/uniprot2allembl.tsv > ref/target_uniprot2target_embl.txt
        cat ref/target_uniprot2target_embl.txt | cut -f3 >ref/target_embl_ids.txt 
        python3 scripts/select_and_clean_entry_from_fasta.py --infasta ref/all_ena_protein_coding.fasta --txToKeep ref/target_embl_ids.txt --outfasta ref/target_protein_seqs.fa
        grep -o -Ff ref/target_embl_ids.txt ref/target_protein_seqs.fa | grep -Ff - ref/target_uniprot2target_embl.txt > ref/target_uniprot_found_in_target_embl.txt
        python3 scripts/clean_labels.py
        cat data/annotation/label_mappings.csv | cut -d',' -f1 > data/annotation/final_embl_ids.txt 
        python3 ~/scripts/select_entry_from_fasta.py --infasta ref/target_protein_seqs.fa --txToKeep  data/annotation/final_embl_ids.txt  --outfasta data/seqs/complete_target_protein_nt_seqs.fa
        '''


rule kmerize_transcripts:
    input: 
        fasta='data/seqs/complete_target_protein_nt_seqs.fa',
    params:
        out_prefix = lambda wildcards: f'data/raw_model_data/uniprot_tx_{wildcards.size}'
    output:  
        full_lsf = 'data/raw_model_data/uniprot_tx_{size}_train.lsf', 
        full_txids = 'data/raw_model_data/uniprot_tx_{size}_train_txids.txt'
    shell:
        '''
        python3  scripts/kmerize_fasta_low_mem.py \
            --infasta {input.fasta} \
            --kmerSize {wildcards.size} \
            --outPfx {params.out_prefix}
           
        '''


'''
03/11/20 changed so we are training on only the target tx, maybe this will help things
'''

rule train_doc2vec_and_infer_vectors:
    input: 
        target_tx = 'data/raw_model_data/uniprot_tx_{size}_train.lsf',
        targ_txids = 'data/raw_model_data/uniprot_tx_{size}_train_txids.txt'
    params: 
        model = lambda wildcards: f'data/embedding_models/uniprot_dm-{wildcards.dm}_wc-{wildcards.mincount}_kmer-{wildcards.size}_dims-{wildcards.dims}.csv.gz',
        model_dir= f'data/embedding_models/'
    output: 
        out_matrix = 'data/embedded_model_data/uniprot_dm-{dm}_wc-{mincount}_kmer-{size}_dims-{dims}.csv.gz'
        
    shell:
        '''
        mkdir -p {params.model_dir}
        python3 scripts/train_doc2vec.py\
            --workingDir {working_dir} \
            --corpusFile {input.target_tx} \
            --corpusTxIDs {input.targ_txids} \
            --edim {wildcards.dims} \
            --dm {wildcards.dm} \
            --wc {wildcards.mincount} \
            --trainedModel {params.model}  \
            --outTrainMatrix {output.out_matrix}
        '''

rule run_random_forest:
    input:
        matrix = 'data/embedded_model_data/uniprot_dm-{dm}_wc-{mincount}_kmer-{size}_dims-{dims}.csv.gz', 
        lab_df = 'data/annotation/label_mappings.csv'
    params: 
        outdir=lambda wildcards: f'data/embedded_model_data/uniprot_dm-{wildcards.dm}_wc-{wildcards.mincount}_kmer-{wildcards.size}_dims-{wildcards.dims}/',
        csvname=lambda wildcards: f'uniprot_rf_dm-{wildcards.dm}_wc-{wildcards.mincount}_kmer-{wildcards.size}_dims-{wildcards.dims}.csv'
    output:
        model_res='data/model_results/uniprot_rf_dm-{dm}_wc-{mincount}_kmer-{size}_dims-{dims}/model_results.csv'
    shell:
        '''
        python3 scripts/run_skl_model.py \
            --workingDir {working_dir} \
            --inputFile {input.matrix} \
            --labFile {input.lab_df} \
            --nproc 32 \
            --outputDir {params.outdir} 
        '''
