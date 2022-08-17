### This is a demo script running downloaded Mageck dataset



TimePoints = ['Day0', 'Day23']
Reps = ['Rep1', 'Rep2']


rule RunMageckCount:
    input: 
        sgRNA_lib = 'resources/mageck-demo/demo/fastq/library.csv',
        fastqs = expand('resources/mageck-demo/demo/fastq/GSC_0131_{timepoint}_{rep}.fastq.gz', timepoint = TimePoints, rep = Reps)
    output: touch('resources/mageck-demo/results/counts/count.done')
    params:
        sample_labels = ','.join(expand('{timepoint}_{rep}', timepoint=TimePoints, rep=Reps)),
        out_prefix = 'resources/mageck-demo/results/counts/GSC_0131'
    shell:
        '''
        mageck count -l {input.sgRNA_lib} \
            --output-prefix {params.out_prefix} \
            --sample-label {params.sample_labels} \
            --fastq {input.fastqs}
        '''


rule RunMageckRRA:
    input: 'resources/mageck-demo/results/counts/GSC_0131.count.txt'
    output: touch('resources/mageck-demo/results/RRA/RRA.done')
    shell:
        '''
        mageck test -k {input} \
            -t Day23_Rep1,Day23_Rep2 \
            -c Day0_Rep1,Day0_Rep2 \
            -n resources/mageck-demo/results/RRA/GSC_0131_RRA \
            --remove-zero both \
            --remove-zero-threshold 0
        '''
