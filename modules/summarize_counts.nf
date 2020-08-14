#!/usr/bin/env nextflow

/**************************************
*********** count reads ***************
**************************************/

process merge_all_counts {
        label 'low_memory'
        tag 'merge_all_counts'
        publishDir "${params.outdir}/", mode: 'copy',
            saveAs: { filename -> 
            if (filename.indexOf(".tsv") > 0 ) "readcounts/$filename"
            else null
            }

        input:
            path readcounts

        output:
            path "readcounts.tsv", emit: ch_readcounts

        script:
            
            """
            merge_counts.py --input ${readcounts} -o readcounts.tsv
            """
}