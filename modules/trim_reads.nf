#!/usr/bin/env nextflow

/********************************************
**************** trimming *******************
********************************************/

if (params.trimmer == 'trimmomatic') {
    process trimmomatic {
        tag "$name"
        label 'process_medium'
        publishDir "${params.outdir}/trimmomatic", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("P.fastq.gz") > 0) "paired/$filename"
                else if (filename.indexOf(".U.fastq.gz") > 0) "unpaired/$filename"
                else if (filename.indexOf(".trimmed.fastq.gz") > 0 ) "single_end/$filename"
                else if (filename.indexOf("summary_stats.txt") > 0) "stats/$filename"
                else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
                else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
                else null
            }

        input:
            tuple val(name), file(reads)
            file wherearemyfiles

        output:
            tuple val(name), file("*.fastq.gz"), emit: ch_trimmed_reads_trimmomatic
            path "*summary_stats.txt", emit: ch_trimmed_trimmo_stats
            file "where_are_my_files.txt"

        script:
            lead = params.leading > 0 ? "LEADING:${params.leading}" : ''
            trail = params.trailing > 0 ? "TRAILING:${params.trailing}" : ''
            min_l = params.min_length > 0 ? "MINLEN:${params.min_length}" : ''
            avg_qual = params.average_qual > 0 ? "${params.average_qual}" : ''

            if (params.singleEnd && !params.adapters) {
                """
                trimmomatic SE \\
                    -threads $task.cpus \\
                    -${params.encoding} \\
                    -trimlog ${name}_trim_log.txt \\
                    -summary ${name}_summary_stats.txt \\
                    $reads \\
                    ${name}.trimmed.fastq.gz \\
                    $lead $trail SLIDINGWINDOW:4:$avg_qual $min_l
                """
            } else if (params.singleEnd && params.adapters) {
                """
                trimmomatic SE \\
                    -threads $task.cpus \\
                    -${params.encoding} \\
                    -trimlog ${name}_trim_log.txt \\
                    -summary ${name}_summary_stats.txt \\
                    $reads \\
                    ${name}.trimmed.fastq.gz \\
                    ILLUMINACLIP:${params.adapters}:2:30:10 \\
                    $lead $trail SLIDINGWINDOW:4:$avg_qual $min_l
                """
            } else if (!params.singleEnd && !params.adapters) {
                """
                trimmomatic PE \\
                    -threads $task.cpus \\
                    -${params.encoding} \\
                    -trimlog ${name}_trim_log.txt \\
                    -summary ${name}_summary_stats.txt \\
                    $reads \\
                    ${name}_R1.P.fastq.gz \\
                    ${name}_R1.U.fastq.gz \\
                    ${name}_R2.P.fastq.gz \\
                    ${name}_R2.U.fastq.gz \\
                    $lead $trail SLIDINGWINDOW:4:$avg_qual $min_l
                """
            } else {
                """
                trimmomatic PE \\
                    -threads $task.cpus \\
                    -${params.encoding} \\
                    -trimlog ${name}_trim_log.txt \\
                    -summary ${name}_summary_stats.txt \\
                    $reads \\
                    ${name}_R1.P.fastq.gz \\
                    ${name}_R1.U.fastq.gz \\
                    ${name}_R2.P.fastq.gz \\
                    ${name}_R2.U.fastq.gz \\
                    ILLUMINACLIP:${params.adapters}:2:30:10 \\
                    $lead $trail SLIDINGWINDOW:4:$avg_qual $min_l
                """
        }
    }
}