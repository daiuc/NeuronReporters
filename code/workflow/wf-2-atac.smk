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


rule PlotFragmentSizes:
    input: rules.MarkDups_atac.output.bam
    output: 'results/ATACseq/PlotFragmentSizes/{timepoint}_{rep}.pdf'
    threads: 1
    resources: cpu = 1, mem_mb = 30000, time = 2100
    script: 'scripts/ATACseqLibraryQC.R'


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


def getMultiBigwigSummaryParamsCmd(wildcards):
    if wildcards.FeatureType.upper() == 'GENOME':
        cmd = "multiBigwigSummary bins"
    else:
        cmd = "multiBigwigSummary BED-file"
    return cmd

def getMultiBigwigSummaryParamsBed(wildcards):
    if wildcards.FeatureType.upper() == 'GENOME':
        bed = ""
    else:
        bed = '--BED ' + config[wildcards.FeatureType]
    return bed

def getMultiBigwigSummaryBinSize(wildcards):
    if wildcards.FeatureType.upper() == 'GENOME':
        bs = 10000
    else:
        bs = 100
    return bs

rule MultiBigwigSummaryATAC:
    message: '### use multiBigwigSummary Genomewide'
    input: list(set(expand("results/ATACseq/bigwig/{timepoint}_{rep}.bw", zip, timepoint=atac.timepoint, rep=atac.rep)))
    output: 
        npz = 'results/ATACseq/BigwigSummary/{FeatureType}/results.npz',
        heatmap = 'results/ATACseq/BigwigSummary/{FeatureType}/heatmap.pdf',
        matrix = 'results/ATACseq/BigwigSummary/{FeatureType}/corr.matrix'
    params: 
       cmd = getMultiBigwigSummaryParamsCmd,
       bed = getMultiBigwigSummaryParamsBed,
       binSize = getMultiBigwigSummaryBinSize
    threads: 8
    resources: cpu = 8, mem_mb = 25000, time = 2100
    shell: 
        '''
        {params.cmd} -b {input} -o {output.npz} {params.bed} -p {threads}  -bs {params.binSize}
        
        plotCorrelation -in {output.npz} -c spearman -p heatmap \
            -o {output.heatmap} --plotNumbers \
            --outFileCorMatrix {output.matrix} \
            -T "Spearman Correlation - {wildcards.FeatureType}" --skipZeros --colorMap viridis
        '''


use rule MultiBigwigSummaryATAC as MultiBigwigSummaryRNA with:
    input: 
        list(set(expand('results/RNAseq/bigwig/{treatment}_{timepoint}_{rep}_{batch}.bw', zip, treatment=rna.treatment, timepoint=rna.timepoint, rep=rna.rep, batch=rna.batch)))
    output: 
        npz = 'results/RNAseq/BigwigSummary/{FeatureType}/results.npz',
        heatmap = 'results/RNAseq/BigwigSummary/{FeatureType}/heatmap.pdf',
        matrix = 'results/RNAseq/BigwigSummary/{FeatureType}/corr.matrix'




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



        




##################################################
####       below analysis are for revision     ###
##################################################



rule getIntergenicAnnotations:
    input: 'resources/annotations/hs38/gencode_v31_protein_u2k_d1k.bed'
    output: 'resources/annotations/hs38/gencode_v31_intergeneic_TSSup100kb-TSSup2kb.bed'
    params:
        CHROM_SIZES = config['HG38_GENOME_SIZES']
    shell: 
        '''
        slopBed -i {input} -g {params.CHROM_SIZES} -l 98000 -r 0 -s `# extend protein coding gene body by 100kb upstream` | \
            subtractBed -a stdin -b {input} -s -F 0.95 `# use 100kb extended gene body to substract gene body` | \
            sort -k1,1 -k2,2 -V > {output}
        '''

rule getIntergenicAndGeneBodyAnnotations:
    input: 'resources/annotations/hs38/gencode_v31_protein_u2k_d1k.bed'
    output: 'resources/annotations/hs38/gencode_v31_protein_u100k_d1k.bed'
    params:
        CHROM_SIZES = config['HG38_GENOME_SIZES']
    shell:
        '''
        slopBed -i {input} -g {params.CHROM_SIZES} -l 98000 -r 0 -s `# extend protein coding gene body by 100kb upstream` | \
            sort -k1,1 -k2,2 -V | \
            awk 'BEGIN  {{OFS="\t"; print "#chrom", "start", "end", "gene_name", "score", "strand"}}
                        {{print $0}}' > {output}
        '''

# def getComputeInteractionsFeatures(wildcards):
#     if wildcards.interact_feature == 'GeneBody':
#         feature = config['GeneBody']
#     elif wildcards.interact_feature == 'Intergenic':
#         feature = config['Intergenic']
#     elif wildcards.interact_feature == 'IntergenicAndGeneBody':
#         feature = config['IntergenicAndGeneBody']
#     return feature

def getComputeInteractionsFeatures(wildcards):
    return config[wildcards.interact_feature]


 
rule ComputeInteractions:
    input: 
        peaks = config['ATAC_CONSENSUS_PEAKS'],
        atac_counts = config['ATAC_CONSENSUS_PEAKS_READS'],
        jaspar2020 = config['JASPAR2020'],
        tflist = config['TF_LIST'],
        hitlist = config['HIT_LIST'],
        features = getComputeInteractionsFeatures,
        deseq2 = config['DESEQ_RESULT'],
        rna_sample_names = config['RNA_SAMPLE_NAMES']
    output: 
        expanded = temp('results/ATACseq/Interactions/{interact_feature}_expanded.tsv')
    wildcard_constraints: 
        interact_feature = 'GeneBody|IntergenicAndGeneBody'
    threads: 16
    resources: cpu = 16, mem_mb = 36000, time = 2100
    script: 'scripts/getInteractions.R'

rule ZipInteractionsFile:
    input: 
        expanded = rules.ComputeInteractions.output.expanded
    output: 
        expanded = 'results/ATACseq/Interactions/{interact_feature}_expanded.tsv.gz'
    shell: 
        '''
        bgzip -c {input.expanded} > {output.expanded}
        '''



plot_regions = {
    "NEUROG2": "chr4:112,512,907-112,516,272",
    "ZBTB18": "chr1:244,047,412-244,061,348",
    "FOXN2": "chr2:48,298,529-48,395,449",
    "POU3F2": "chr6:98,833,946-98,839,782",
    "VAX2": "chr2:70,896,924-70,934,684",
    "INSM1": "chr20:20,367,631-20,371,231"
}
rule plotPyGenomeTracks:
    '''manually run snakemake to plot regions in `plot_regions`
    '''
    input: "results/genometracks/{gene_name}/tracks.ini"
    output: "results/genometracks/{gene_name}/tracks.pdf"
    params:
        region = lambda wcs: plot_regions[wcs.gene_name]
    threads: 1
    shell:
        '''
        pyGenomeTracks --tracks {input} \
            --region {params.region} \
            --title {wildcards.gene_name} \
            -o {output}
        '''



