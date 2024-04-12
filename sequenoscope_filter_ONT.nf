#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_fastq = ''
params.input_summary = ''
params.output = ''
params.output_prefix = 'sample'
params.classification = ''
params.min_ch = 1 
params.max_ch = 512
params.min_dur = 0
params.max_dur = 100
params.min_start = 0
params.max_start = 259200
params.min_q = 0
params.max_q = 100
params.min_len = 0
params.max_len = 50000
params.summarize = false
params.force = false

workflow {
    input_ch = Channel.of([params.input_fastq.split(',\\s*'), params.input_summary, params.output])

    process RunSequenoscopeAnalyze {
        input:
        tuple val(input_fastq), val(input_summary), val(output)

        output:
        path("${output}")

        script:
        def opts = " -op ${params.output_prefix}"
        if(params.classification) opts += " --classification ${params.classification}"
        if(params.min_ch != 1) opts += " -min_ch ${params.min_ch}"
        if(params.max_ch != 512) opts += " -max_ch ${params.max_ch}"
        if(params.min_dur > 0) opts += " -min_dur ${params.min_dur}"
        if(params.max_dur < 100) opts += " -max_dur ${params.max_dur}"
        if(params.min_start > 0) opts += " -min_start ${params.min_start}"
        if(params.max_start < 259200) opts += " -max_start ${params.max_start}"
        if(params.min_q > 0) opts += " -min_q ${params.min_q}"
        if(params.max_q < 100) opts += " -max_q ${params.max_q}"
        if(params.min_len > 0) opts += " -min_len ${params.min_len}"
        if(params.max_len < 50000) opts += " -max_len ${params.max_len}"
        if(params.summarize) opts += " --summarize"
        if(params.force) opts += " --force"

        def fastq_inputs = input_fastq.join(' ')
        """
        sequenoscope filter_ONT --input_fastq $fastq_inputs --input_summary $input_summary -o $output $opts
        """
    }

    RunSequenoscopeAnalyze(input_ch)
}