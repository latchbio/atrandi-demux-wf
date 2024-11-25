import glob
import os


# def aggregate_input(wildcards):
#     '''
#     aggregate the file names of the random number of files
#     generated at the scatter step
#     '''

#     return expand(os.path.join(results_dir, "coverage/{{sample}}_coverage/{i}_coverage.txt"),
#         i=glob_wildcards(os.path.join(results_dir,"coverage/{sample}_coverage/{i}_coverage.txt")).i) #double {{}} for sample wildcards

def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.split_full.get(**wildcards).output[1] #directory of results, but not a pheniqs report
    return expand(os.path.join(results_dir, "coverage/{{sample}}_coverage/{i}_coverage.txt"),
        i=glob_wildcards(os.path.join(checkpoint_output, '{i}_r1.fastq.gz')).i) #double {{}} for sample wildcards


# pysam is shit, doesn't have normal working containers and it and all it's future descendants should be cursed upon forever. Also I hate my life
rule calculate_coverage:
    input:
        aligned_bam = rules.minimap2_align.output.aligned_bam
    output:
        coverage_tables = os.path.join(results_dir,"coverage/{sample}_coverage/{i}_coverage.txt"),
    threads: 192
    container:
        "docker://quay.io/biocontainers/mulled-v2-0fd299cadb7a80e2cc704b5d903ccc54893c512d:a1db556b69552bda96fc53f4abaaca85c8b8c212-0"
    shell:
        "samtools coverage {input.aligned_bam} > {output.coverage_tables}"

rule aggregate_tables:
    input:
        coverage_tables = aggregate_input
    output:
        single_coverage_table = os.path.join(results_dir,"coverage/{sample}_coverage.txt")
    script:
        "../scripts/calcualte_coverage_samtools.py"


# also adds barcode lookup
rule add_foldchange:
    input:
        coverage_table = rules.aggregate_tables.output.single_coverage_table,
        lookup_table = rules.parse_pheniqs_report.output.split_barcode_table
    output:
        fold_change_table = os.path.join(results_dir,"coverage/{sample}_foldchange.txt")
    params:
        filter_chromosomes = config['PARAMS']['coverage_params']['filter_chromosomes'],
    container:
        "docker://quay.io/biocontainers/mulled-v2-0fd299cadb7a80e2cc704b5d903ccc54893c512d:a1db556b69552bda96fc53f4abaaca85c8b8c212-0"
    script:
        "../scripts/calculate_foldchange.py"


