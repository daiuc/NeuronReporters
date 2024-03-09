# upload sra

def getUpload_to_SRA_params(dataset):
    if dataset == 'rna':
        return {
            'local_dir' : 'resources/RNAseq/fastq',
            'file_pattern' : '*.fastq.gz',
            'remote_dir' : 'rna-seq'
        }
    elif dataset == 'atac':
        return {
            'local_dir' : 'resources/ATACseq/fastq',
            'file_pattern' : '*.fastq.gz',
            'remote_dir': 'atac-seq'
        }
    elif dataset == 'crispr':
        return {
            'local_dir' : 'resources/crispr/counts',
            'file_pattern' : 'RC.csv',
            'remote_dir': 'crispr-screen'
        }


rule upload_to_sra:
    '''
    FTP server, username, password, ftp_folder must have been exported to env
    must mkdir within ftp_folder first
    '''
    output: touch("resources/upload_{dataset}_to_sra.done")
    wildcard_constraints: dataset = 'rna|atac|crispr'
    params:
        script = 'workflow/scripts/upload_to_sra.sh',
        local_dir = lambda wildcards: getUpload_to_SRA_params(wildcards.dataset)['local_dir'],
        remote_dir = lambda wildcards: getUpload_to_SRA_params(wildcards.dataset)['remote_dir'],
        file_pattern = lambda wildcards: getUpload_to_SRA_params(wildcards.dataset)['file_pattern']
    shell:
        '''
        LOCAL_FOLDER={params.local_dir}
        REMOTE_FOLDER={params.remote_dir}
        FILE_PATTERN={params.file_pattern}
        
        # Upload files to FTP server
        for file in "$LOCAL_FOLDER"/$FILE_PATTERN; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                echo "Uploading $file to ${{REMOTE_FOLDER}}/$filename on ftp..."
                ftp -n "$FTP_SERVER" <<EOF
                user "$FTP_USERNAME" "$FTP_PASSWORD"
                cd "$FTP_FOLDER"/$REMOTE_FOLDER
                put "$file" "$filename"
                quit
EOF
            echo "Done."
            fi
        done

        echo "All files uploaded successfully."


        '''
