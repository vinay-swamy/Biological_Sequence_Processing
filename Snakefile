'''
Snakefile for project



'''

genome='ref/gencode_genome.fa.gz'
transcripts='ref/gencode_transcript_seqs.fa.gz'
ano='ref/gencode_ano.gtf.gz'

rule all:
    input: 'ref/nontranscribed_windows.bed',dtx='ref/dummy_tx.bed'

rule download_annotation:
    output:genome, transcripts, ano
    shell:
        '''
        mkdir -p ref/
        wget -O {genome} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
        wget -O {ano} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz
        wget -O {transcripts} ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.transcripts.fa.gz
        '''


'''
- extract transcript locations, make a padded version, and make 2 bit
- get masked regions of genome(masked means unsequenced)

'''

rule get_genome_masked_regions:
    input: genome, transcripts, ano
    output: tbit='ref/gencode_genome.2bit',chrom_masks='ref/chrom_masked_regions.bed' 
    shell:
        '''
        module load ucsc
        faToTwoBit {genome} {output.tbit} 
        twoBitInfo {output.tbit} -nBed stdout > {output.chrom_masks}
        '''


rule get_chrom_sizes:
    input: genome, transcripts, ano
    output: bed = 'ref/gencode_chroms.bed', bt_g = 'ref/chrom_lengths_bt.txt'
    shell:
        '''
        python3 scripts/chromFasta2Bed.py {genome} {output.bed}
        cut {output.bed} -f1,3 > {output.bt_g}
        '''

rule get_transcript_info:
    input: bt_g = 'ref/chrom_lengths_bt.txt', a=ano, t=transcripts
    output: tloc = 'ref/transcript_locs.bed', tloc_pad = 'ref/transcript_locs_padded.bed', tx_l='ref/transcript_lengths.bed'
    shell:
        '''
        module load bedtools
        zcat {ano} | awk '$3 == "transcript"' - | cut -f1,4,5 > {output.tloc}
        bedtools slop -b 100 -i {output.tloc} -g {input.bt_g} > {output.tloc_pad}
        python3 scripts/chromFasta2Bed.py {transcripts} {output.tx_l}
        '''


rule make_valid_chrom_intervals:
    input: chrom_bed = 'ref/gencode_chroms.bed', masks = 'ref/chrom_masked_regions.bed', tloc_pad = 'ref/transcript_locs_padded.bed'
    output: ntw='ref/nontranscribed_windows.bed'
    shell:
        '''
        module load bedtools
        bedtools subtract -a {input.chrom_bed} -b {input.masks} | 
            bedtools subtract -a stdin -b {input.tloc_pad} > {output}
        '''

# https://stackoverflow.com/questions/5914513/shuffling-lines-of-a-file-with-a-fixed-seed

rule make_dummy_tx:
    input: tx_l = 'ref/transcript_lengths.bed', wins = 'ref/nontranscribed_windows.bed', bt_g = 'ref/chrom_lengths_bt.txt'
    output: shufs='ref/wins_to_shuffle.bed',  dtx='ref/dummy_tx.bed',
    shell:
        '''
        module load bedtools        
        python3 scripts/make_dummy_tx.py {input.tx_l} {output.shufs}
        bedtools shuffle -seed 42 -incl {input.wins} -noOverlapping -i {output.shufs} -g {input.bt_g} > /tmp/dummy.bed
        k=`wc -l  < /tmp/dummy.bed`
        for i in $(seq 1 $k); do echo "dummy_${{i}}" >> /tmp/dummy.txt ; done
        paste /tmp/dummy.bed /tmp/dummy.txt > {output.dtx}
        '''
