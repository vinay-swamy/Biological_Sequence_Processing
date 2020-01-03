'''
Snakefile for project



'''

# def wc2fa(wc, k ,l):
#     if w  == 'PC

genome='ref/gencode_genome.fa'
transcripts='ref/gencode_transcript_seqs.fa'
ano='ref/gencode_ano.gtf.gz'




rule all:
    #input: expand('run_k-{k}_l-{l}/dummy_transcript_seqs.fa',k=[6, 10, 14, 20], l=[1, 2] )
    input: expand('run_k-{k}_l-{l}/data/PC_NC_Y_vec.pydata',k=[6, 10], l=[1, 2]),\
        expand('run_k-{k}_l-{l}/data/ref_dummy_Y_vec.pydata',k=[6, 10], l=[1, 2])

rule download_annotation:
    output:genome, transcripts, ano
    shell:
        '''
        mkdir -p ref/
        wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz | gunzip -c - > {genome}
        wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz | gunzip -c - > {ano}
        wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz | gunzip -c - > {transcripts}
        module load samtools 
        samtools faidx {genome}
       
        '''


'''
-remove transcripts that dont meet length req
-find and remove masked regions from genome 
-pad locations of known transcripts by pad_lengthon either side, and then remove the padd locatiosn from the genome
-get size distibution of konwn transcripts length
'''

rule prep_ref_data:
    input:raw_tx=transcripts, genome=genome, ano=ano, 
    params: min_transcript_length = config['min_transcript_length']
    output: \
        filt_tx = 'ref/ref_transcript_seqs.fa',\
        tbit = 'ref/gencode_genome.2bit',\
        chrom_masks = 'ref/chrom_masked_regions.bed',\
        chrom_bed = 'ref/gencode_chroms.bed',\
        bt_g = 'ref/chrom_lengths_bt.txt',\
        tloc = 'ref/transcript_locs.bed',\
        tloc_pad = 'ref/transcript_locs_padded.bed',\
        tx_len = 'ref/transcript_lengths.bed',\
        ntw = 'ref/nontranscribed_windows.bed'
    shell:
        '''
        python3 scripts/remove_tx_by_lengths.py {input.raw_tx} {params.min_transcript_length} {output.filt_tx}
        
        module load ucsc
        faToTwoBit {input.genome} {output.tbit} 
        twoBitInfo {output.tbit} -nBed stdout > {output.chrom_masks}
        python3 scripts/chromFasta2Bed.py {genome} {output.chrom_bed}
        cut {output.chrom_bed} -f1,3 > {output.bt_g}
        
        module load bedtools
        awk '$3 == "transcript"' {input.ano}  | cut -f1,4,5 > {output.tloc}
        bedtools slop -b 100 -i {output.tloc} -g {output.bt_g} > {output.tloc_pad}
        python3 scripts/chromFasta2Bed.py {output.filt_tx} {output.tx_len}
        
        module load bedtools
        bedtools subtract -a {output.chrom_bed} -b {output.chrom_masks} |
            bedtools subtract -a stdin -b {output.tloc_pad} > {output.ntw}
        
        '''


'''
randomly sample lengths annotated transcript length distribution, and then overlay those on to the non transcribd windows we found earlier (could not figure out how to get gnu shuf to take a random seed  )
'''
rule make_dummy_tx:
    input: tx_l = 'ref/transcript_lengths.bed', wins = 'ref/nontranscribed_windows.bed', bt_g = 'ref/chrom_lengths_bt.txt', tloc_pad = 'ref/transcript_locs_padded.bed', fa = genome
    params: NUM_DUM_TX = config['max_dummy_transcripts'], dummy_pf = 'ref', min_transcript_length=config['min_transcript_length']
    output: shufs = 'ref/shuf.bed', dtx = 'ref/dummy_tx.bed', dseq='ref/dummy_transcript_seqs.fa'
    shell:
        '''
        db=/tmp/{params.dummy_pf}.bed
        dt=/tmp/{params.dummy_pf}.txt
        rm -rf $db $dt
        module load bedtools        
        python3 scripts/make_dummy_tx.py {input.tx_l} {params.NUM_DUM_TX} {output.shufs}
        bedtools shuffle -seed 42 -incl {input.wins} -noOverlapping -i {output.shufs} -g {input.bt_g} |\
             bedtools intersect -v -a stdin -b {input.tloc_pad} > $db
        k=`wc -l  < $db`
        for i in $(seq 1 $k); do echo "dummy_${{i}}" >> $dt  ; done
        paste $db $dt > {output.dtx}
        bedtools getfasta -fi {input.fa} -bed {output.dtx}  > /tmp/tmp.fasta
        python3 scripts/remove_tx_by_lengths.py /tmp/tmp.fasta {params.min_transcript_length} {output.dseq}
        
        '''


rule split_tx_PC_NC:
    input: all_tx= 'ref/ref_transcript_seqs.fa', pctx='ref/gencode_protein_coding_transcripts.txt', nctx='ref/gencode_non_coding_transcripts.txt'
    params: min_transcript_length=config['min_transcript_length']
    output:PC='ref/PC_transcript_seqs.fa', NC='ref/NC_transcript_seqs.fa'
    shell:
        '''
        python3 scripts/select_entries_from_fasta.py {input.all_tx} {params.min_transcript_length} {input.pctx} {output.PC}
        python3 scripts/select_entries_from_fasta.py {input.all_tx} {params.min_transcript_length} {input.nctx} {output.NC}
        
        '''


rule kmerize_transcripts:
    input: dummy = 'ref/{second}_transcript_seqs.fa', ref = 'ref/{first}_transcript_seqs.fa'
    output: kmer_pydata = 'run_k-{k}_l-{l}/data/{first}_{second}_kmers.pydata', kmer_lineSentence = 'run_k-{k}_l-{l}/data/{first}_{second}_kmers.lsf'
    shell:
        '''
        which python3
        python3  scripts/kmerize_fasta_low_mem.py {input.ref}  {input.dummy} {wildcards.k} {wildcards.l} {output.kmer_pydata} {output.kmer_lineSentence}
        '''


rule train_doc2vec:
    input: corpus='run_k-{k}_l-{l}/data/ref_dummy_kmers.lsf'
    output: model = 'run_k-{k}_l-{l}/models/doc2vec_ep-15_PV-DBOW_M-300.pymodel'
    shell:
        '''
        which python3 
        python3 scripts/train_doc2vec.py {input.corpus} {output.model}
        '''

rule infer_vectors:
    input: model = 'run_k-{k}_l-{l}/models/doc2vec_ep-15_PV-DBOW_M-300.pymodel', data = 'run_k-{k}_l-{l}/data/{first}_{second}_kmers.pydata'
    output: X_data = 'run_k-{k}_l-{l}/data/{first}_{second}_X_mat.pydata', Y_data = 'run_k-{k}_l-{l}/data/{first}_{second}_Y_vec.pydata'
    shell:
        '''
        which python3 
        python3 scripts/infer_vectors.py {input.data} {input.model} {output.X_data} {output.Y_data}

        '''


