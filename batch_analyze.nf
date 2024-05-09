#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_dir = '' // Directory containing FASTQ files
params.input_reference = ''
params.output_dir = '' // Base directory for outputs
params.seq_sum = null  // Optional, null by default
params.seq_type = 'SE'  // Default to Single-End, change to 'PE' for Paired-End
params.start_time = 0
params.end_time = 100
params.output_prefix = 'sample'
params.threads = 1
params.minimum_read_length = 15
params.maximum_read_length = 0
params.trim_front_bp = 0
params.trim_tail_bp = 0
params.quality_threshold = 15
params.minimum_coverage = 1
params.minimap2_kmer = 15
params.force = false

process RunSequenoscopeAnalyze {
    tag "${sample_id}" // For easy identification of jobs
    publishDir "${params.output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fastqs), path(input_reference)

    output:
    path "${sample_id}/*"

    script:
    def opts = []
    if (params.start_time != 0) opts += "--start_time ${params.start_time}"
    if (params.end_time != 100) opts += "--end_time ${params.end_time}"
    if (params.output_prefix) opts += "-op ${params.output_prefix}"
    if (params.threads) opts += "-t ${params.threads}"
    if (params.minimum_read_length != 15) opts += "-min_len ${params.minimum_read_length}"
    if (params.maximum_read_length != 0) opts += "-max_len ${params.maximum_read_length}"
    if (params.trim_front_bp != 0) opts += "-trm_fr ${params.trim_front_bp}"
    if (params.trim_tail_bp != 0) opts += "-trm_tail ${params.trim_tail_bp}"
    if (params.quality_threshold != 15) opts += "-q ${params.quality_threshold}"
    if (params.minimum_coverage != 1) opts += "-min_cov ${params.minimum_coverage}"
    if (params.minimap2_kmer != 15) opts += "--minimap2_kmer ${params.minimap2_kmer}"
    if (params.force) opts += "--force"
    if (params.seq_sum != null && params.seq_type == 'SE') opts += "-seq_sum ${params.seq_sum}"
    """
    mkdir -p ${sample_id}
    sequenoscope analyze --input_fastq ${fastqs.join(' ')} --input_reference $input_reference -o ${sample_id} -seq_type ${params.seq_type} ${opts.join(' ')}
    """
}

workflow {
    def pattern = (params.seq_type == 'PE') ? '*_{1,2}.fastq' : '*.fastq'
    fastq_ch = Channel.fromFilePairs("$params.fastq_dir/$pattern", size: (params.seq_type == 'PE' ? 2 : 1))
    input_reference_ch = Channel.value(file(params.input_reference))

    fastq_ch
        .combine(input_reference_ch)
        .set{ combined_ch }

    RunSequenoscopeAnalyze(combined_ch)
}