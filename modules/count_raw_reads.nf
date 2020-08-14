#!/usr/bin/env nextflow

/**************************************
*********** count reads ***************
**************************************/

process count_raw_reads {
    tag "$name"
    label "count_reads"
    publishDir "${params.outdir}/qc", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".tsv") > 0 ? "raw_counts/$filename" : "$filename" }
    
    input:
        tuple val(name), file(reads) 

    output:
        tuple val(name), file("*.counts.tsv"), emit: ch_read_counts
        
    script:        
        """
        count-reads.py --infile ${reads[0]} --output ${name}.counts.tsv
        """
}

process merge_raw_counts {
        label 'low_memory'
        tag 'merge_raw_counts'
        publishDir "${params.outdir}/", mode: 'copy',
            saveAs: { filename -> 
            if (filename.indexOf(".tsv") > 0 ) "readcounts/$filename"
            else null
            }

        input:
            file readcounts

        output:
            path "raw_readcounts.tsv", emit: ch_merged_counts

        script:
            
            """
            merge_read_counts.py --input ${readcounts} -o raw_readcounts.tsv
            """
    }