# ./minimap2 -ax sr ref.fa read1.fa read2.fa > aln.sam 

# minimap2 -d ref.mmi ref.fa                     # indexing
# minimap2 -a ref.mmi reads.fq > alignment.sam   # alignment


rule minimap2_make_index:
    input:
        fasta = config['FILEPATHS']['REF']['fasta'],
    output:
        minimap2_index = config['FILEPATHS']['REF']['index_file']
    container:
        "https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0"
    threads: config['PARAMS']['thread_params']['threads']
    shell:
        """
        minimap2 -d {output.minimap2_index} {input.fasta} 
        """


rule minimap2_align:
    input:
        minimap2_index_file = rules.minimap2_make_index.output.minimap2_index,
        fastq_files_r1 = os.path.join(results_dir,"pheniqs/{sample}_split_fastq/{i}_r1.fastq.gz"),
        fastq_files_r2 = os.path.join(results_dir,"pheniqs/{sample}_split_fastq/{i}_r2.fastq.gz")
    output:
        # aligned_sam = directory(os.path.join(results_dir,"minimap2/{sample}_alignment/"))
        aligned_bam = os.path.join(results_dir,"minimap2/{sample}_alignment/{i}_alignment.bam")
    container:
        # "https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0"
        "https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:f83d48eb33ccc96dde5a2fb495d999e16ff0bc69-0"
    threads: config['PARAMS']['thread_params']['threads']
    # check how piping+threads works in snakemake
    #test if piping+threads gives the same result with test data
    shell:
        """
        minimap2 -ax sr {input.minimap2_index_file} {input.fastq_files_r1} {input.fastq_files_r2} -t {threads} | samtools sort -O bam -o {output.aligned_bam} -@ {threads}
        """

