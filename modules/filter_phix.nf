#!/usr/bin/env nextflow

/******************************************************************************
******************** Preprocessing: build phix genome index *******************
******************************************************************************/

if (params.removePhiX) {
    phix_fasta = 'RTA'? params.genomes[ 'RTA' ].fasta ?: false : false
    phix_bwa_index = 'RTA' ? params.genomes[ 'RTA' ].bwa ?: false : false
    phix_bowtie2_index = 'RTA' ? params.genomes[ 'RTA' ].bowtie2 ?: false : false

    if ( phix_bwa_index && params.aligner == 'bwa' ){
        lastPath = phix_bwa_index.lastIndexOf(File.separator)
        bwa_dir =  phix_bwa_index.substring(0,lastPath)
        bwa_base = phix_bwa_index.substring(lastPath+1)
        ch_phix_bwa_index = Channel
            .fromPath(phix_bwa_index)
            .ifEmpty { exit 1, "BWA index not found: ${phix_bwa_index} you can download the files from (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz) and extract in the igenomes dir" }
    }
    if (phix_fasta && params.aligner == 'bwa') {
        ch_phix_fasta = Channel
            .fromPath(phix_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${phix_fasta}, you can download the files from (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz) and extract in the igenomes dir" }
    }
}


process build_bwa_index_phix {
    label 'build_index'
    tag "$fasta"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
        saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
        file fasta

    output:
        file 'phix_bwa_index'

    script:

        """
        bwa index -a bwtsw ${fasta}
        mkdir phix_bwa_index && mv ${fasta}* phix_bwa_index
        """
}

/*********************************************************************************
*************** phix: alignment, sorting, indexing and extracting ****************
*********************************************************************************/


process bwa_align_phix {
    label 'bwa_align'
    tag "$name"
    publishDir "${params.outdir}/bwa", mode: 'copy', saveAs: {filename ->
        if (filename.indexOf(".bam") > 0) "alignments/phix/$filename"
        else if (!params.saveAlignedIntermediates && filename == "where_are_my_files.txt") filename
        else if (params.saveAlignedIntermediates && filename != "where_are_my_files.txt") filename
        else null
        }
    
    input:
        tuple val(name), file(reads)            
        file index
        file wherearemyfiles

    output:
        tuple val(name), path("*.bam"), emit: ch_phix_bams
        file "where_are_my_files.txt"


    script:
        prefix="${name}"
        rg="\'@RG\\tID:${name}\\tSM:${name.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${name}\\tPU:1\'"

        if (params.singleEnd) {
        """
        bwa mem -t ${task.cpus} -M -R $rg ${bwa_dir}/${index} $reads > ${prefix}.bam
        """
        } else {
        """
        bwa mem -t ${task.cpus} -M -R $rg ${bwa_dir}/${index} $reads \\
        | samtools view -@ ${task.cpus} -bS -O BAM -o ${prefix}.bam -
        """
        }
}

process sort_bam_phix {
    label 'high_memory'
    tag "${name}"
    publishDir "${params.outdir}/bwa", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") > 0) "sorted/phix/$filename"
            else null
        }

    input:
        tuple val(name), file(ch_phix_bwa_bam)

    output:
        tuple val(name), path("${name}*.sorted.bam"), emit: ch_phix_sorted_bams

    script:
        """
        samtools sort \\
            $ch_phix_bwa_bam \\
            -@ ${task.cpus} \\
            -o ${name}.sorted.bam
            samtools index ${name}.sorted.bam
        """
}

process extract_unmapped_phix {
    label 'extract_unmapped'
    tag "${name}"
    publishDir "${params.outdir}/bwa", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0) "unmapped/phix/$filename"
            else null
        }

    input:
        tuple val(name), file(ch_phix_sorted_bam)

    output:
        tuple val(name), path("${name}*.fastq"), emit: ch_phix_bam_unmapped_reads

    script:
        if (params.singleEnd) {
        """
        samtools fastq \\
            --threads ${task.cpus} \\
            -f 4 \\
            $ch_phix_sorted_bam > ${name}.fastq
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
            $ch_phix_sorted_bam
        """
    }
}

