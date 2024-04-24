#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_fastq = ''  // For single-end or the first file of paired-end
params.input_fastq2 = '' // Only for the second file of paired-end
params.input_reference = ''
params.output = ''
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

workflow {
    // Depending on the sequencing type, handle input creation
    def input_fastqs = (params.seq_type == 'PE') ?
                       [params.input_fastq, params.input_fastq2] :
                       [params.input_fastq]
    input_ch = Channel.of([input_fastqs, params.input_reference, params.output, params.seq_sum, params.seq_type])

    process RunSequenoscopeAnalyze {
        input:
        tuple val(input_fastqs), val(input_reference), val(output), val(seq_sum), val(seq_type)

        output:
        path("${output}")

        script:
        def opts = ""
        if (params.start_time != 0) opts += " --start_time ${params.start_time}"
        if (params.end_time != 100) opts += " --end_time ${params.end_time}"
        if (params.output_prefix) opts += " -op ${params.output_prefix}"
        if (params.threads) opts += " -t ${params.threads}"
        if (params.minimum_read_length != 15) opts += " -min_len ${params.minimum_read_length}"
        if (params.maximum_read_length != 0) opts += " -max_len ${params.maximum_read_length}"
        if (params.trim_front_bp != 0) opts += " -trm_fr ${params.trim_front_bp}"
        if (params.trim_tail_bp != 0) opts += " -trm_tail ${params.trim_tail_bp}"
        if (params.quality_threshold != 15) opts += " -q ${params.quality_threshold}"
        if (params.minimum_coverage != 1) opts += " -min_cov ${params.minimum_coverage}"
        if (params.minimap2_kmer != 15) opts += " --minimap2_kmer ${params.minimap2_kmer}"
        if (params.force) opts += " --force"

        def fastq_inputs = input_fastqs.join(' ')
        def seq_sum_option = params.seq_sum ? "-seq_sum $seq_sum" : ""

        // Build the command according to the sequencing type
        def command = ""
        if (params.seq_type == 'PE') {
            command = "sequenoscope analyze --input_fastq $fastq_inputs --input_reference $input_reference -o $output -seq_type $seq_type $opts"
        } else {
            command = "sequenoscope analyze --input_fastq $fastq_inputs --input_reference $input_reference -o $output -seq_type $seq_type $seq_sum_option $opts"
        }
        """
        $command
        """
    }

    RunSequenoscopeAnalyze(input_ch)
}
