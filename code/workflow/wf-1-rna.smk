



##------------- RNA-seq -----------------------

## Note: We used 3' enriched RNA-seq protocols. 
## Theoretically, all captured molecules should have a polyA tail, 
## 3'UTR and the last Exons only. When normalizing library size, 
## we should NOT normalize on gene length. A counts per million should suffice.

rule FastQC_RNA:
    input: FASTQC_RNA_INPUTS
    output: touch('results/RNAseq/fastqc/fasqc.done')
    params:
        outdir = "results/RNAseq/fastqc"
    threads: 12
    resources: time=1000, mem_mb=40000, cpu=12
    shell:
        '''
        module load fastqc/0.11.4
        
        touch {output}
        fastqc --threads {threads} \
            --outdir {params.outdir} {input}
        '''

rule MultiQC_RNA:
    input: rules.FastQC_RNA.output
    output: 'results/RNAseq/multiqc/multiqc_report.html'
    params:
        fastqc_dir = 'results/RNAseq/fastqc',
        out_dir = 'results/RNAseq/multiqc'
    threads: 1
    shell:
        '''
        module load multiqc
        ls {input}
        multiqc --outdir {params.out_dir} {params.fastqc_dir}
        '''


rule CutPolyA:
    input: "resources/RNAseq/fastq/{rna_fastq}.fastq.gz"
    output: temp("resources/RNAseq/fastq/PolyATrimmed/{rna_fastq}.fastq.gz")
    threads:4
    params: 
        adapter="AAAAAAA", 
        min_read_length=25, 
        min_overlap=5
    resources: time=480, mem_mb=30000, cpu=4
    shell:
        '''
        module load cutadapt/2.10
        cutadapt -a {params.adapter} \
            -m {params.min_read_length} -O {params.min_overlap} -j {threads}\
            -o {output} {input}
        '''




rule STAR_RNA:
    input: get_STAR_inputs
    output: 
        transBam = "results/RNAseq/STAR/{treatment}_{timepoint}_{rep}_{batch}.Aligned.toTranscriptome.out.bam",
        genomBam = temp("results/RNAseq/STAR/{treatment}_{timepoint}_{rep}_{batch}.Aligned.sortedByCoord.out.bam")
    params:
        genome_dir = config['STAR_INDEX'],
        out_prefix = "results/RNAseq/STAR/{treatment}_{timepoint}_{rep}_{batch}.",
        RG = '{treatment}_{timepoint}_{rep}_{batch}',
        read_files_comm = "zcat"
    threads: 8
    resources: time=2000, mem_mb=42000, cpu=8
    shell:
        '''
        module load star/2.7.1a 
        echo STAR {params.out_prefix}
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input} \
            --readFilesCommand {params.read_files_comm} \
            --outFilterType BySJout \
            --outSAMattributes NH HI AS NM MD RG \
            --outSAMunmapped None \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM \
            --outFileNamePrefix {params.out_prefix} \
            --limitBAMsortRAM 400000000000 \
            --outSAMattrRGline ID:{params.RG} SM:{params.RG} PL:ILLUMINA LB:{params.RG} 
        
        sleep 10
        ls {output.transBam}
        #samtools index -@ {threads} {output.transBam}
        samtools index -@ {threads} {output.genomBam}
        '''

# filter bam
rule filter_bam_rna:
    input: rules.STAR_RNA.output.genomBam
    output: temp('results/RNAseq/STAR/{treatment}_{timepoint}_{rep}_{batch}_filtered.bam')
    params:
        chroms = CHROMS
    threads: 4
    resources: time=100, mem_mb=20000, cpu=4
    shell:
        '''
        samtools view -b -h -q 20 -F 4 -@ {threads} -o {output} {input} {params.chroms}
        samtools index -@ {threads} {output} 
        '''

# mark duplicates, marking both pcr and optical duplicates, but removing sequencing duplicates only.
rule MarkDups_rna:
    input:
        bam = rules.filter_bam_rna.output,
    output:
        bam = "results/RNAseq/MarkDups/{treatment}_{timepoint}_{rep}_{batch}.bam",
        metrics = "results/RNAseq/MarkDups/{treatment}_{timepoint}_{rep}_{batch}_metrics.txt"
    params:
        ref = config['FA_HS38']
    threads: 1
    resources: time=600, mem_mb=40000, cpu=1
    shell:
        '''
        module load gatk/4.1.8.1
        gatk MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            -R {params.ref} \
            --REMOVE_SEQUENCING_DUPLICATES \
            --CREATE_INDEX
        '''

rule Bigwig_rna:
    input: rules.MarkDups_rna.output.bam
    output: "results/RNAseq/bigwig/{treatment}_{timepoint}_{rep}_{batch}.bw"
    threads:8
    resources: time=600, mem_mb=36000, cpu=8
    shell:
        '''
         module load deeptools/3.4.2


         bamCoverage --bam {input} \
            -o {output} \
            --outFileFormat bigwig \
            --binSize 10 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize 3049315783 \
            --ignoreForNormalization chrX \
            --numberOfProcessors {threads} \
            --minMappingQuality 20
        '''


# rule RSEM
rule RSEM:
    input: rules.STAR_RNA.output.transBam
    output: "results/RNAseq/RSEM/{treatment}_{timepoint}_{rep}_{batch}.genes.results"
    params:
        rsem_idx = config['RSEM_INDEX'],
        sample_name = "results/RNAseq/RSEM/{treatment}_{timepoint}_{rep}_{batch}"
    threads: 8
    resources: time=1000, mem_mb=40000, cpu=8
    shell:
        '''
        module load rsem/1.2.21
        rsem-calculate-expression \
            --bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345 --time \
            -p {threads} --ci-memory 30000 \
            {input} {params.rsem_idx} {params.sample_name}
        ls {output}
        '''



rule combineRsemCounts:
    input: getCombineRsemCountsInputs
    output: 'results/RNAseq/RSEM/combined_RSEM_counts.txt'
    threads: 1
    script: 'scripts/collectRsemCounts.R'

# # QC stats of RNAseq reads and alignments
# rule RNAseqQC:
#     input:
#         star_final_log = expand("results/RNAseq/STAR/{batch}_S{sample}.Log.final.out", batch=RNA_BATCHES, sample=RNA_SAMPLES)
#     output: "results/RNAseq/QC/RNAseqAlignedReadsQC.html"  # rmarkdown render
#     params:
#         output_dir="results/RNAseq/QC"
#     threads: 4
#     resources: time=600, mem_mb=36000, cpu=4
#     script:
#         "Scripts/RNAseqAlignedReadsQC.Rmd"

    