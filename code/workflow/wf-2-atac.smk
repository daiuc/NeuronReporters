#Note modify the rule to use bwa alignment instead. 
# see parameters in bam files in TFscreen/atac/bams_v3

# ID:bwa  PN:bwa  VN:0.7.17-r1194-dirty   CL:
# bwa sampe -P /gpfs/commons/groups/sanjana_lab/cdai/ref_genome/bwa_gencode_genome/gencode_GRCh38.primary_assembly.genome.fa \
# /gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/fastq/ATAC1.1.sai \
# /gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/fastq/ATAC1.2.sai \
# /gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/fastq/HG7FNBGX9_n01_ATAC1_11-14-18.fastq \
# /gpfs/commons/groups/sanjana_lab/cdai/TFscreen/atac/fastq/HG7FNBGX9_n02_ATAC1_11-14-18.fastq


# bwa aln -t $nThreadsBWA $REF $read1 > ${SAMPLE}.1.sai
# bwa aln -t $nThreadsBWA $REF $read2 > ${SAMPLE}.2.sai
# bwa sampe -P $REF ${SAMPLE}.1.sai ${SAMPLE}.2.sai $read1 $read2 > ${SAMPLE}.sam


# ###### Run Samtools to sort sam to sorted bam, and index

# echo "Sorting, indexing bam..."
# samtools sort -o ${SAMPLE}.sorted.bam -@ $nThreadsSAM -T ${SAMPLE} -@ $nThreadsSAM ${SAMPLE}.sam
# samtools index ${SAMPLE}.sorted.bam




##-------------- ATAC-seq -------------------

rule FASTQC_ATAC:
    input: FASTQC_ATAC_INPUTS
    output: touch('results/ATACseq/fastqc/fasqc.done')
    #output: "results/RNAseq/fastqc/{batch}_S{sample}/{batch}_S{sample}_fastqc.out"
    params:
        outdir = "results/ATACseq/fastqc"
    threads: 12
    resources: time=1000, mem_mb=40000, cpu=12
    shell:
        '''
        module load fastqc/0.11.4
        
        touch {output}
        fastqc --threads {threads} \
            --outdir {params.outdir} {input}
        '''

rule MultiQC_ATAC:
    input: rules.FASTQC_ATAC.output
    output: 'results/ATACseq/multiqc/multiqc_report.html'
    params:
        fastqc_dir = 'results/ATACseq/fastqc',
        out_dir = 'results/ATACseq/multiqc'
    threads: 1
    shell:
        '''
        module load multiqc
        ls {input}
        multiqc --outdir {params.out_dir} {params.fastqc_dir}
        '''



rule bwa:
    input: unpack(get_bwa_inputs)
    output: temp('results/ATACseq/bwa/{timepoint}_{rep}.sam')
    params:
        bwa_refseq = config['FA_HS38_BWA'],
    threads: 8
    resources: time=1800, mem_mb=40000, cpu=8
    shell:
        '''
        module load bwa/0.7.17

        bwa aln -t {threads} {params.bwa_refseq} {input.R1} \
            > results/ATACseq/bwa/{wildcards.timepoint}_{wildcards.rep}.1.sai
        bwa aln -t {threads} {params.bwa_refseq} {input.R2} \
            > results/ATACseq/bwa/{wildcards.timepoint}_{wildcards.rep}.2.sai 
        bwa sampe -P {params.bwa_refseq} \
            results/ATACseq/bwa/{wildcards.timepoint}_{wildcards.rep}.1.sai \
            results/ATACseq/bwa/{wildcards.timepoint}_{wildcards.rep}.2.sai \
            {input.R1} {input.R2} > {output}
        
        rm results/ATACseq/bwa/{wildcards.timepoint}_{wildcards.rep}.1.sai \
            results/ATACseq/bwa/{wildcards.timepoint}_{wildcards.rep}.2.sai  

        '''

rule sort_bam_atac:
    input: rules.bwa.output
    output: temp('results/ATACseq/bwa/{timepoint}_{rep}.bam')
    threads: 4
    resources: time=100, mem_mb=20000, cpu=4
    shell:
        '''
        samtools sort \
            -o {output} \
            -@ {threads} \
            -T {wildcards.timepoint}_{wildcards.rep}  \
            {input}
        
        samtools index -@ {threads} {output}
        '''

# filter mapq >= 20, chr in chr1 - chrY
rule filter_bam_atac:
    input: rules.sort_bam_atac.output
    output: temp('results/ATACseq/bwa/{timepoint}_{rep}_filtered.bam')
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
rule MarkDups_atac:
    input:
        bam = rules.filter_bam_atac.output,
    output:
        bam = "results/ATACseq/MarkDups/{timepoint}_{rep}.bam",
        metrics = "results/ATACseq/MarkDups/{timepoint}_{rep}_metrics.txt"
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

rule Bigwig_atac:
    input: rules.MarkDups_atac.output.bam
    output: "results/ATACseq/bigwig/{timepoint}_{rep}.bw"
    threads:8
    resources: time=600, mem_mb=36000, cpu=8
    shell:
        '''
         module load deeptools/3.4.2
         bamCoverage --bam {input} \
            -o {output} \
            --outFileFormat bigwig \
            --binSize 1 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize 3049315783 \
            --ignoreForNormalization chrX \
            --numberOfProcessors {threads} \
            --minMappingQuality 20
        '''


rule CallPeaks:
    input: rules.MarkDups_atac.output.bam # calling peak on each sample
    output: "results/ATACseq/macs2/{timepoint}_{rep}_peaks.narrowPeak"
    params:
        out_dir = "results/ATACseq/macs2/",
        prefix = "{timepoint}_{rep}"
    resources: time=600, mem_mb=30000, cpu=1
    shell:
        '''
        module load macs2/2.1.2
        macs2 callpeak \
            -t {input} \
            -f BAMPE \
            --outdir {params.out_dir} \
            --name {params.prefix} \
            -g 3.1e9 \
            -q 0.05 \
            --cutoff-analysis

        ls -lah {output}
        '''


# Construct consensus peaks
rule consensusPeaks:
    input: expand("results/ATACseq/macs2/{timepoint}_{rep}_peaks.narrowPeak", zip, \
                timepoint=atac.timepoint, rep=atac.rep)
    output: 
        bed = 'results/ATACseq/diffbind/consensusPeaks.bed',
        counts = 'results/ATACseq/diffbind/consensusPeaks_counts.txt' 
    params:
        MIN_OVERLAP = 2,
        SAMPLESHEET = config['DIFFBIND_SAMPLES']
    resources: time=60, mem_mb=30000, cpu=1
    script: "scripts/diffBind.R"



# rule ATACseqQC:
#     input:
#         bam = expand("results/ATACseq/MarkDups/atac_{batch}_{sample}_markDups.bam", batch=ATAC_BATCHES, sample=ATAC_SAMPLES)
#     output:
#         "results/ATACseq/QC/ATACseqLibraryQC.html"
#     params:
#         output_dir="results/ATACseq/QC"
#     threads: 2
#     resources: time=600, mem_mb=36000, cpu=4
#     script:
#         "Scripts/ATACseqLibraryQC.Rmd" # could also just use this script as input
        





