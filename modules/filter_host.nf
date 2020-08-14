#!/usr/bin/env nextflow

/*********************************************************************************
*************** host: alignment, sorting, indexing and extracting ****************
*********************************************************************************/

if (params.reference) {
    lastPath = params.reference.lastIndexOf(File.separator)
    ref_bwa_dir  = params.reference.substring(0,lastPath)
    ref_bwa_base = params.reference.substring(lastPath+1)
    ch_reference_fasta = Channel
        .fromPath(params.reference)
        .ifEmpty { exit 1, "Reference Fasta file not found: ${params.reference}" }
}


if (params.reference && params.aligner == 'bwa') {
    process build_bwa_index_ref {
        label 'high_memory'
        tag "$fasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
        saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
            file fasta

        output:
            file "*.{amb,ann,bwt,pac,sa}"

        script:
            
            """
            bwa index -a bwtsw ${fasta}
            """
    }


    process bwa_align_ref {
        label 'high_memory'
        tag "$name"
        publishDir "${params.outdir}/bwa", mode: 'copy', 
            saveAs: { filename ->
                if (filename.indexOf(".bam") > 0) "alignments/host/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
                }

        input:
            tuple val(name), file(reads)
            file fasta
            file index
            file wherearemyfiles
        
        output:
            tuple val(name), file("${name}.bam"),  emit: ch_ref_bams
            val "where_are_my_files.txt"

        script:
            rg="\'@RG\\tID:${name}\\tSM:${name.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\'"

            if (params.singleEnd) {
                """
                bwa mem -t ${task.cpus} -M -R $rg ${fasta} $reads > ${name}.bam
                """
            } else {
                """
                bwa mem -t ${task.cpus} -M -R $rg ${fasta} $reads \\
                | samtools view -@ ${task.cpus} -bS -O BAM -o ${name}.bam -
                """
            }
    }

    process sort_bam_ref {
        label 'high_memory'
        tag "$name"
        publishDir "${params.outdir}/bwa", mode: 'copy',
            saveAs: { filename ->
                if (filename.indexOf(".bam") > 0) "sorted/host/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
                }

        input:
            tuple val(name), ch_bam

        output:
            tuple val(name), path("${name}*.sorted.bam"), emit: ch_ref_sorted_bams

        script:
            """
            samtools sort \\
                $ch_bam \\
                -@ ${task.cpus} \\
                -o ${name}.sorted.bam
            samtools index ${name}.sorted.bam
            """
    }

    process extract_unmapped_ref {
        label 'high_memory'
        tag "${name}"
        publishDir "${params.outdir}/bwa", mode: 'copy',
            saveAs: { filename ->
                if (filename.indexOf(".fastq") > 0) "unmapped/host/$filename"
                else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
                else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
                else null
                }

        input:
            tuple val(name), file(ch_ref_sorted_bam)

        output:
            tuple val(name), path("${name}*.fastq"), emit: ch_ref_bam_unmapped_reads

        script:
            if (params.singleEnd) {
                """
                samtools fastq \\
                --threads ${task.cpus} \\
                -f 4 \\
                $ch_ref_sorted_bam > ${name}.fastq
                """
            } else {
                """
                samtools fastq \\
                --threads ${task.cpus} \\
                -1 ${name}_R1.fastq \\
                -2 ${name}_R2.fastq \\
                -0 /dev/null \\
                -s /dev/null \\
                -f 4 \\
                $ch_ref_sorted_bam
                """
            }
    }
}




