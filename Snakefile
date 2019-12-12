'''
Snakefile for project



'''

genome='ref/gencode_genome.fa'
transcripts='ref/gencode_transcript_seqs.fa'
ano='ref/gencode_ano.gtf.gz'

rule all:
    #input: expand('run_k-{k}_l-{l}/dummy_transcript_seqs.fa',k=[6, 10, 14, 20], l=[1, 2] )
    input: expand('run_k-{k}_l-{l}/data/X_mat.pydata',k=[6, 10, 14, 20], l=[1, 2])

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
- extract transcript locations, make a padded version, and make 2 bit
- get masked regions of genome(masked means unsequenced)

'''

rule remove_small_transcripts:
    input: transcripts
    params: min_transcript_length = config['min_transcript_length']
    output: 'run_k-{k}_l-{l}/gencode_transcript_seqs_valid_size.fa'
    shell:
        '''
        
        python3 scripts/remove_tx_by_lengths.py {input} {params.min_transcript_length} {output}

        '''



rule get_genome_info:
    input: genome
    output: tbit = 'ref/gencode_genome.2bit', chrom_masks = 'ref/chrom_masked_regions.bed',  bed = 'ref/gencode_chroms.bed', bt_g = 'ref/chrom_lengths_bt.txt'
    shell:
        '''
        module load ucsc
        faToTwoBit {genome} {output.tbit} 
        twoBitInfo {output.tbit} -nBed stdout > {output.chrom_masks}
        python3 scripts/chromFasta2Bed.py {genome} {output.bed}
        cut {output.bed} -f1,3 > {output.bt_g}
        '''


'''
pad annotated transcripts and get distribution of transcript lengths 
'''


rule get_transcript_info:
    input: bt_g = 'ref/chrom_lengths_bt.txt', a = ano, transcripts = 'run_k-{k}_l-{l}/gencode_transcript_seqs_valid_size.fa'
    output: tloc = 'run_k-{k}_l-{l}/transcript_locs.bed', tloc_pad = 'run_k-{k}_l-{l}/transcript_locs_padded.bed', tx_l = 'run_k-{k}_l-{l}/transcript_lengths.bed'
    shell:
        '''
        module load bedtools
        awk '$3 == "transcript"' {ano}  | cut -f1,4,5 > {output.tloc}
        bedtools slop -b 100 -i {output.tloc} -g {input.bt_g} > {output.tloc_pad}
        python3 scripts/chromFasta2Bed.py {input.transcripts} {output.tx_l}
        '''


rule make_valid_chrom_intervals:
    input: chrom_bed = 'ref/gencode_chroms.bed', masks = 'ref/chrom_masked_regions.bed', tloc_pad = 'run_k-{k}_l-{l}/transcript_locs_padded.bed'
    output: ntw = 'run_k-{k}_l-{l}/nontranscribed_windows.bed'
    shell:
        '''
        module load bedtools
        bedtools subtract -a {input.chrom_bed} -b {input.masks} | 
            bedtools subtract -a stdin -b {input.tloc_pad} > {output}
        '''

# https://stackoverflow.com/questions/5914513/shuffling-lines-of-a-file-with-a-fixed-seed

'''
randomly sample lengths annotated transcript length distribution, and then overlay those on to the non transcribd windows we found earlier (could not figure out how to get gnu shuf to take a random seed  )
'''
rule make_dummy_tx:
    input: tx_l = 'run_k-{k}_l-{l}/transcript_lengths.bed', wins = 'run_k-{k}_l-{l}/nontranscribed_windows.bed', bt_g = 'ref/chrom_lengths_bt.txt', tloc_pad = 'run_k-{k}_l-{l}/transcript_locs_padded.bed'
    params: NUM_DUM_TX = config['max_dummy_transcripts'], dummy_pf = lambda wildcards: f'run_k-{wildcards.k}_l-{wildcards.l}'
    output: shufs = 'run_k-{k}_l-{l}/shuf.bed', dtx = 'run_k-{k}_l-{l}/dummy_tx.bed',
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
        '''

rule make_dummy_tx_fasta:
    input: fa = genome, bed = 'run_k-{k}_l-{l}/dummy_tx.bed'
    params: min_transcript_length = config['min_transcript_length'], tfile = lambda wildcards: f'/tmp/run_k-{wildcards.k}_l-{wildcards.l}.bed'
    output: 'run_k-{k}_l-{l}/dummy_transcript_seqs.fa'
    shell:
        '''
        module load bedtools 
        bedtools getfasta -fi {input.fa} -bed {input.bed}  > {params.tfile}
        python3 scripts/remove_tx_by_lengths.py {params.tfile} {params.min_transcript_length} {output}

        '''

rule kmerize_transcripts:
    input: dummy = 'run_k-{k}_l-{l}/dummy_transcript_seqs.fa', ref = 'run_k-{k}_l-{l}/gencode_transcript_seqs_valid_size.fa'
    output: kmer_pydata = 'run_k-{k}_l-{l}/data/all_record_kmers.pydata', kmer_lineSentence = 'run_k-{k}_l-{l}/data/all_record_kmers.lsf'
    shell:
        '''
        python3  scripts/kmerize_fasta_low_mem.py {input.ref}  {input.dummy} {wildcards.k} {wildcards.l} {output.kmer_pydata} {output.kmer_lineSentence}
        '''


rule train_doc2vec:
    input: corpus='run_k-{k}_l-{l}/data/all_record_kmers.lsf'
    output: model = 'run_k-{k}_l-{l}/models/doc2vec_ep-15_PV-DBOW_M-300.pymodel'
    shell:
        '''
        python3 scripts/train_doc2vec.py {input.corpus} {output.model}
        '''

rule infer_vectors:
    input: model = 'run_k-{k}_l-{l}/models/doc2vec_ep-15_PV-DBOW_M-300.pymodel', data = 'run_k-{k}_l-{l}/data/all_record_kmers.pydata'
    output: X_data = 'run_k-{k}_l-{l}/data/X_mat.pydata', Y_data = 'run_k-{k}_l-{l}/data/Y_vec.pydata'
    shell:
        '''
        python3 scripts/infer_vectors.py {input.data} {input.model} {output.X_data} {output.Y_data}

        '''
