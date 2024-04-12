#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.test_dir = ''
params.control_dir = ''
params.output = ''
params.op = 'sample'
params.adaptive_sampling = false // Optional parameter
params.single = false
params.comparison_metric = 'taxon_%_covered_bases' // Corrected typo
params.violin_data_percent = 0.1
params.bin = 'minutes'
params.legend = false
params.force = false

workflow {
    input_ch = Channel.of([params.test_dir, params.control_dir, params.output])

    process RunSequenoscopePlot {
        input:
        tuple val(test_dir), val(control_dir), val(output)

        output:
        path("${output}")

        script:
        def opts = params.adaptive_sampling ? " -AS" : ""
        if(params.output_prefix) opts += " -op ${params.op}"
        if(params.single) opts += " -single"
        if(params.comparison_metric)opts += " --comparison_metric ${params.comparison_metric}"
        if(params.violin_data_percent)opts += " -VP ${params.violin_data_percent}"
        if(params.bin)opts += " -bin ${params.bin}"
        if(params.legend) opts += "--legend"
        if(params.force) opts += " --force

        """
        sequenoscope plot -T $test_dir -C $control_dir -o $output $opts
        """
    }

    RunSequenoscopePlot(input_ch)
}
