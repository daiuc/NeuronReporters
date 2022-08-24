# this is for reanalysis of crispr-n screen using MAGeCK




rule WrangleCrisprRawCountTable:
    input: config['CRISPR_COUNTS']
    output: 
        raw = 'resources/crispr/counts/wrangled_raw_counts.tsv',
        norm = 'resources/crispr/counts/wrangled_ranknorm_counts.tsv'
    threads: 1
    resources: cpu = 1 , mem_mb = 15000, time = 2100
    script: 'scripts/prepCrisprRC.R'


def getMageckRRA_input(wildcards):
    return 'resources/crispr/counts/' + wildcards.crispr_contrast + '.counts.tsv'

def getMageckRRA_treatment(wildcards):
    return wildcards.crispr_contrast.split('-vs-')[0]

def getMageckRRA_control(wildcards):
    return wildcards.crispr_contrast.split('-vs-')[1]


rule MageckRRA:
    # input: getMageckRRA_input
    input: rules.WrangleCrisprRawCountTable.output.raw
    output: touch('results/crispr/mageckRRA/{crispr_contrast}/mageckRRA.done')
    wildcard_constraints: 
        crispr_contrast = 'Day[0-9]{1,2}_(Tdpos|Tdneg|iN)-vs-Day[0-9]{1,2}_(Tdpos|Tdneg|iN)'
    params: 
        treatment = getMageckRRA_treatment,
        control = getMageckRRA_control,
        day0_label = getMageckRRA_control,
        out_prefix = 'results/crispr/mageckRRA/{crispr_contrast}/RRA'
    shell: 
        '''
        mageck test -k {input} \
            -t {params.treatment} -c {params.control} \
            -n {params.out_prefix} \
            --remove-zero both --remove-zero-threshold 0 \
            --norm-method median \
            --sort-criteria neg
        '''

rule MageckMLE:
    input: 
        counts = rules.WrangleCrisprRawCountTable.output.raw,
        
    output: touch('results/crispr/mageckRRA/{crispr_contrast}/mageckRRA.done')
    params: 
        treatment = getMageckRRA_treatment,
        control = getMageckRRA_control,
        day0_label = getMageckRRA_control,
        out_prefix = 'results/crispr/mageckRRA/{crispr_contrast}/RRA'
    shell: 
        '''
        mageck mle -k {input} \
            -d ""
        '''