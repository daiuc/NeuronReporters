# download ChIP-seq for ZEB1 Chip-Seq on biopolar spindle neuron Homosapiens
# https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8660358/SRR8660358
# GEO:GSE127650,FactorBook:ENCSR418KUS

rule Download_ChIP:
    message: 'Download SRA of ChIP-seq'
    output: 
        temp('resources/ChIP/{Gene}/SRR8660358')
    log: 'logs/Download_{Gene}_ChIP.log'
    params: 
        url = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8660358/SRR8660358'
    threads: 1 
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
        curl -L {params.url} -s -o {output}
        '''

rule fasterq_dump:
    message: 'fasterq-dump SRR'
    input: 
        rules.Download_ChIP.output
    output: 
        touch('resources/ChIP/{Gene}/fastq_dump.done')
    log: 'logs/fastq_dump_{Gene}.log'
    params: 
        out_dir = 'resources/ChIP/{Gene}'
    threads: 8 
    resources: cpu = 8 , mem_mb = 25000, time = 2100
    shell: 
        '''
        fasterq-dump \
            -e {threads} -O {params.out_dir} -3 \
            {input} 
        '''

def get_chip_fastq(wildcards):
    fastqs = glob.glob1('resources/ChIP/' + wildcards.Gene, '*fastq.gz')
    fastqs = ['resources/ChIP/' + wildcards.Gene + '/' + q for q in fastqs]
    # print(fastqs)
    return fastqs

rule FASTQC_CHIP:
    input: get_chip_fastq
    output: touch('results/ChIP/{Gene}/fastqc/fastqc.done')
    params:
        outdir = "results/ChIP/{Gene}/fastqc"
    threads: 12
    resources: time=1000, mem_mb=40000, cpu=12
    shell:
        '''
        fastqc="/nfs/sw/fastqc/fastqc-0.11.4/fastqc"
        
        touch {output}
        fastqc --threads {threads} \
            --outdir {params.outdir} {input}
        '''

def get_AlignChip_input(wildcards):
    fastqs = glob.glob1('resources/ChIP/' + wildcards.Gene, '*fastq.gz')
    fastqs = ['resources/ChIP/' + wildcards.Gene + '/' + q for q in fastqs]
    return fastqs

rule AlignChip:
    input: get_AlignChip_input
    output: temp('results/ChIP/{Gene}/bwa/chip.sam')
    wildcard_constraints: 
        Gene = '[A-Z0-9]+'
    params:
        bwa_refseq = config['FA_HS38_BWA']
    threads: 12
    resources: time=1800, mem_mb=36000, cpu=12
    shell:
        '''
        module load bwa/0.7.17

        bwa aln -t {threads} {params.bwa_refseq} {input} \
            > results/ChIP/{wildcards.Gene}/bwa/chip.sai
        bwa samse {params.bwa_refseq} \
            results/ChIP/{wildcards.Gene}/bwa/chip.sai {input} > {output}
                
        rm results/ChIP/{wildcards.Gene}/bwa/chip.sai

        '''
rule SortBamChip:
    input: rules.AlignChip.output
    output: temp('results/ChIP/{Gene}/bwa/chip_raw.bam')
    params:
        flagstat = 'results/ChIP/{Gene}/bwa/chip_raw_flagstat.txt'
    threads: 8
    resources: time=100, mem_mb=20000, cpu=8
    shell:
        '''
        samtools view -Shu {input} | \
        samtools sort \
            -o {output} \
            -@ {threads} \
            --write-index -
        samtools flagstat -@ {threads} {output} > {params.flagstat} 
        '''

# filter mapq >= 20, chr in chr1 - chrY
rule FilterBamChip:
    input: rules.SortBamChip.output
    output: 'results/ChIP/{Gene}/bwa/chip_filtered.bam'
    params:
        chroms = CHROMS,
        flagstat = 'results/ChIP/{Gene}/bwa/chip_filtered_flagstat.txt'
    threads: 4
    resources: time=100, mem_mb=20000, cpu=4
    shell:
        '''
        samtools view -b -h -q 20 -F 4 -@ {threads} -o {output} {input} {params.chroms}
        samtools index -@ {threads} {output}
        samtools flagstat -@ {threads} {output} > {params.flagstat} 
        '''

rule BigwigChip:
    input: rules.FilterBamChip.output
    output: "results/ChIP/{Gene}/deeptools/chip.bw"
    threads:8
    resources: time=600, mem_mb=36000, cpu=8
    shell:
        '''
         bamCoverage --bam {input} \
            -o {output} \
            --outFileFormat bigwig \
            --binSize 5 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize 3049315783 \
            --ignoreForNormalization chrX \
            --numberOfProcessors {threads} \
            --minMappingQuality 20
        '''


rule CallPeaksChip:
    input: rules.FilterBamChip.output
    output: "results/ChIP/{Gene}/macs2/chip_peaks.narrowPeak"
    params:
        out_dir = "results/ChIP/{Gene}/macs2",
        prefix = "chip"
    resources: time=600, mem_mb=30000, cpu=1
    shell:
        '''
        macs2 callpeak \
            -t {input} \
            -f BAM \
            --outdir {params.out_dir} \
            --name {params.prefix} \
            -g 3.1e9 \
            -q 0.05 \
            --keep-dup auto \
            --cutoff-analysis
        ls -lah {output}
        '''


rule IntersectChipPeaks:
    message: 'Intersect Chip peaks to find targets'
    input: 
        peaks = rules.CallPeaksChip.output,
        tf = 'resources/annotations/hs38/gencode_v31_tf_u2k_d1k.bed'
    output: 'results/ChIP/{Gene}/targets.bed'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
        intersectBed  -wo \
            -a {input.peaks} -b {input.tf} > {output}
        '''