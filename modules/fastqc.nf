#!/usr/bin/env nextflow

/**************************************
*************** fastqc ****************
**************************************/

process fastqc {
    tag "$name"
    label "fastqc"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }

    when:
        !params.skipQC && !params.skipFastQC
    
    input:
        tuple val(name), file(reads) 

    output:
        path "*_fastqc.{zip, html}", emit: ch_fastqc_results

    script:
        name = reads.baseName
        """
        fastqc --quiet --nogroup --format fastq --threads $task.cpus $reads
        """
}