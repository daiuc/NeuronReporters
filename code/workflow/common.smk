

rule MakeAnnotationBedFiles:
    message: 'Make Annotation BED files from GENCODE'
    input: 
        gencode = config['GENCODE_HS38'],
        tflist = config['TF_LIST']
    output: 
        protein = 'resources/annotations/hs38/gencode_v31_protein.bed',
        tf = 'resources/annotations/hs38/gencode_v31_tf.bed',
    threads: 1
    resources: cpu = 1, mem_mb = 25000, time = 2100
    script: 'scripts/makeBedFromGencodeAnnotation.R'


rule ExtendAnntationBedFile:
    message: 'Extend up 2kb and down 1kb of gene coordinates'
    input: 
        protein = rules.MakeAnnotationBedFiles.output.protein,
        tf = rules.MakeAnnotationBedFiles.output.tf
    output: 
        protein = 'resources/annotations/hs38/gencode_v31_protein_u2k_d1k.bed',
        tf = 'resources/annotations/hs38/gencode_v31_tf_u2k_d1k.bed'
    params:
        genome = config['HG38_GENOME_SIZES']
    threads: 1
    resources: cpu = 1, mem_mb = 15000, time = 2100
    shell: 
        '''
        slopBed -l 2000 -r 1000 -s -header \
            -i {input.protein} -g {params.genome} > {output.protein}
        slopBed -l 2000 -r 1000 -s -header \
            -i {input.tf} -g {params.genome} > {output.tf}
        '''