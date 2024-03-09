# QC mostly addressing reviewer comments

# read in previously aligned bam files including rna-seq, atac-seq
# this is in the TFscreen directory
previous_bams = pd.read_csv(config['PREVIOUS_BAMS'], 
                                sep='\t').set_index(['library', 'batch', 'sample'], 
                                                    drop=False)


def getFlagstatBamsInput(wildcards):
    '''
    given wildcards, return the absolute path of the bam file
    eg: /gpfs/commons/groups/sanjana_lab/cdai/TFscreen/data/H2HY7BGXC/star_results/bams/H2HY7BGXC_n01_S01.Aligned.toTranscriptome.out.bam
    '''
    path = previous_bams.query(
        'library == @wildcards.library & batch == @wildcards.batch & \
        sample == @wildcards.sample').paths[0]
    return(path)

rule getFlagstatBams:
    message: 'Get read library stats for QC'
    input: getFlagstatBamsInput
    output: 'results/QC/flagstat/{library}/{batch}_{sample}.txt'
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
        samtools flagstat -@ {threads} -O tsv {input} > {output}
        '''