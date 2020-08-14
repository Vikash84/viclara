/******************************************
************** multiqc reports  ***********
******************************************/

// Has the run name been specified by the user? this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    
    when:
    !params.skipMultiQC
    
    input:
        file multiqc_config
        path ('fastqc/*')
        //path ('trimmomatic/stats*')
        file workflow_summary

    output:
        file "*multiqc_report.html" //into multiqc_report
        file "*_data"
        file "multiqc_plots"

    script:
        rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
        rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
        """
        multiqc . -f $rtitle $rfilename --config $multiqc_config -m custom_content -m fastqc
        """
}


