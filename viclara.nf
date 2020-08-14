#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* ViClaRA -- Virus Classification and Reference-based Assembly
* Author: ajodeh.juma@gmail.com
*/

/******************************************************************
*          Help messages, inputs and validation checks            *
******************************************************************/

if (!nextflow.version.matches('20.+')) {
    println "This workflow requires nextflow version 20.XX.XX and above -- Your nextflow version $nextflow.version"
    exit 1
}

if (params.help) {
    exit 0, helpMessage()
}

if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (!params.reads) {
    exit 1, "input missing, use [--reads]"}
if (!params.reference) {
    exit 1, "input missing, use [--reference]"}

// Has the run name been specified by the user? this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}


/**********************************************************
***************** channels for config files ***************
**********************************************************/

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
//ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/************************************************
************* header log info *******************
************************************************/

def summary = [:]
summary['Help Message']             = params.help
summary['Pipeline Name']            = "V I C L A R A"
summary['Pipeline Version']         = params.version
summary['Run Name']                 = custom_runName ?: workflow.runName
summary['Data Type']                = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Output dir']               = params.outdir
summary['Launch dir']               = workflow.launchDir
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['User']                     = workflow.userName
summary['Config Profile']           = workflow.profile
summary['Reference']                = params.reference
if (params.adapters) {
    summary['Adapter sequences']    = params.adapters
}
if (params.trimmer == 'trimmomatic') {
    summary['Trimmer']              = 'trimmomatic'
    summary['Trimming']             = "leading: $params.leading / trailing: $params.trailing / encoding: $params.encoding / minimum_length: $params.min_length / average quality: $params.average_qual"
}
if (params.trimmer == 'bbduk') {
    summary['Trimmer']              = 'bbduk'
}
if (params.aligner == 'bwa') {
    summary['Aligner']              = 'bwa'
}
summary['Save prefs']               = "Ref Genome: "+(params.saveReference ? 'Yes' : 'No')+" / Trimmed FastQ: "+(params.saveTrimmed ? 'Yes' : 'No')+" / Alignment intermediates: "+(params.saveAlignedIntermediates ? 'Yes' : 'No')
summary['Remove Phix174']           = params.removePhiX
summary['Human genome']             = params.hg
if (params.classify && params.krn2_db && params.krn2_task && params.krn2_library){
  summary['Classify Reads']         = params.classify
  summary['Kraken2 Task']           = params.krn2_task
  summary['Kraken2 database']       = params.krn2_db
  summary['Kraken2 Library']        = params.krn2_library
}

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "---------------------------------------------------------------"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'virclara-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'viclara Workflow Summary'
    section_href: 'https://github.com/ajodeh-juma/viclara'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

// define text colours
c_green = "\033[0;32m";
c_yellow = "\033[0;33m";
c_red = "\u001B[31m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";
c_dim = "\033[2m";

/*****************************************************************
********** create channel for input fastq read files *************
******************************************************************/

if (params.readPaths) {
    if (params.singleEnd) {
        ch_raw_reads = Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0], checkIfExists: true)]]}
            .dump()
            .ifEmpty {exit 1, "$params.readPaths was empty - no input files supplied"}
    } else {
        ch_raw_reads = Channel
            .from(params.readPaths)
            .map {row -> [row[0], [file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true)]]}
            .dump()
            .ifEmpty {exit 1, "$params.readPaths was empty - no input files supplied"}
    }
} else {
    if (params.singleEnd) {
        ch_raw_reads = Channel
            .fromPath(params.reads)
            .map { file -> tuple(file.baseName, file) }
    } else {
        ch_raw_reads = Channel
            .fromFilePairs(params.reads, size: -1 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
            }
}


/*********************************************************
*********** create channels for PhiX fasta file **********
*********************************************************/
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
            .ifEmpty { exit 1, "BWA index not found: ${phix_bwa_index}" }
    }
    if (phix_fasta && params.aligner == 'bwa') {
        ch_phix_fasta = Channel
            .fromPath(phix_fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Genome fasta file not found: ${phix_fasta}, you can download the files from (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz) and extract in the igenomes dir" }
    }
}

/*************************************************************************
**************** create channel fo the human genome ******************
*************************************************************************/

if (params.hg) {
    // lastPath = params.hg.lastIndexOf(File.separator)
    // ref_bwa_dir  = params.hg.substring(0,lastPath)
    // ref_bwa_base = params.hg.substring(lastPath+1)
    ch_hg_fasta = Channel
        .fromPath(params.hg)
        .ifEmpty { exit 1, "Reference Fasta file not found: ${params.hg}" }
}

/*************************************************************************
**************** create channel fo the reference genome ******************
*************************************************************************/

if (params.reference) {
    ch_reference_fasta = Channel
        .fromPath(params.reference)
        .ifEmpty { exit 1, "Reference Fasta file not found: ${params.reference}" }
    // lastPath = params.reference.lastIndexOf(File.separator)
    // ref_bwa_dir  = params.reference.substring(0,lastPath)
    // ref_bwa_base = params.reference.substring(lastPath+1)
    // ch_ref_bwa_index = Channel
    //     .fromPath(ref_bwa_dir, checkIfExists: true)
}


/******************************************************
************ validate kraken2 input options ***********
******************************************************/

if (params.classify && params.krn2_task != 'download-taxonomy' && params.krn2_task != 'download-library' && params.krn2_task != 'build' && params.krn2_task != 'standard') {
    exit 1, "Invalid kraken2 task option: ${params.krn2_task}. Valid options: 'download-taxonomy', 'download-library', 'build', 'standard'"
}

// stage channel for kraken2 input
if (params.classify && !params.krn2_db && !params.krn2_task && !params.krn2_library){
  exit 1, "Missing options, database path: ${params.kraken_db}, task: ${params.krn2_task} and library: ${params.krn2_library}"
} else if (params.classify && params.krn2_db && !params.krn2_task && !params.krn2_library) {
  exit 1, "Missing options, task: ${params.krn2_task} and library: ${params.krn2_library}"
} else if (params.classify && params.krn2_db && params.krn2_task && !params.krn2_library) {
  exit 1, "Missing options, library: ${params.krn2_library}"
} else if (params.classify && params.krn2_db && params.krn2_task && params.krn2_library) {
  ch_krn2_db = Channel.fromPath(params.krn2_db, type: 'dir')
} else if (params.classify && params.krn2_db && params.krn2_task == 'standard'){
  ch_krn2_db = Channel.fromPath(params.krn2_db, type: 'dir')
}

/**********************************************************************************
***** channels to check where files are if saving options are set true ************
**********************************************************************************/
ch_where_trimmomatic = Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
ch_where_bwa_bams = Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
ch_where_phix_bwa = Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
ch_where_human_bwa = Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)
ch_where_ref_bwa = Channel.fromPath("$baseDir/assets/where_are_my_files.txt", checkIfExists: true)


/*****************************************
***************** modules ****************
*****************************************/

include {get_software_versions}                                                                     from './modules/software_versions'
include {fastqc}                                                                                    from './modules/fastqc'
include {trimmomatic}                                                                               from './modules/trim_reads'
include {build_bwa_index_phix; bwa_align_phix; sort_bam_phix; extract_unmapped_phix}                from './modules/filter_phix'
include {build_bwa_index_human; bwa_align_human; sort_bam_human; extract_unmapped_human}            from './modules/filter_human'
include {build_bwa_index_ref; bwa_align_ref; sort_bam_ref; extract_unmapped_ref}                    from './modules/filter_host'
include {build_krakendb; classify_reads; merge_kraken_reports; visualize_species_summary}           from './modules/build_krakendb'
include {multiqc}                                                                                   from './modules/multiqc'
include {count_raw_reads; merge_raw_counts}                                                         from './modules/count_raw_reads'
include {count_trimmed_reads; merge_trimmed_counts}                                                 from './modules/count_trimmed_reads'
include {count_filtered_phix_reads;merge_filtered_phix_counts}                                      from './modules/count_filtered_phix_reads'
include {merge_all_counts}                                                                          from './modules/summarize_counts'

/*****************************************
***************** workflow ***************
*****************************************/


workflow {

// parse software version numbers
get_software_versions()

// fastqc
fastqc(ch_raw_reads)

// trimming 
if (!params.skipTrimming) {
    ch_trimmed_reads_trimmo = trimmomatic(ch_raw_reads, ch_where_trimmomatic.collect())
    if (!params.singleEnd) {
        ch_trimmed_reads_trimmomatic = trimmomatic.out.ch_trimmed_reads_trimmomatic.map { row -> [row[0], [ file(row[1][0], checkIfExists: true), file(row[1][2], checkIfExists: true) ] ] } 
    } else {
        ch_trimmed_reads_trimmomatic = trimmomatic.out.ch_trimmed_reads_trimmomatic
    }
} else {
    ch_trimmed_reads_trimmomatic = ch_raw_reads
}

// remove phix sequences to get unmapped reads
if (params.removePhiX && params.aligner == 'bwa') {
    build_bwa_index_phix(ch_phix_fasta)
    bwa_align_phix(ch_trimmed_reads_trimmomatic, ch_phix_bwa_index.collect(), ch_where_phix_bwa.collect())
    sort_bam_phix(bwa_align_phix.out.ch_phix_bams)
    extract_unmapped_phix(sort_bam_phix.out.ch_phix_sorted_bams)    
}

/****************************************************************************
****************** remove human sequences to get unmapped reads **************
****************************************************************************/

if (params.hg && params.aligner == 'bwa') {
    build_bwa_index_human(ch_hg_fasta)
    bwa_align_human(extract_unmapped_phix.out.ch_phix_bam_unmapped_reads, ch_hg_fasta.collect(), build_bwa_index_human.out.collect(), ch_where_human_bwa.collect())
    sort_bam_human(bwa_align_human.out.ch_human_bams)
    extract_unmapped_human(sort_bam_human.out.ch_human_sorted_bams)
}

/****************************************************************************
****************** remove host sequences to get unmapped reads **************
****************************************************************************/

if (params.reference && params.aligner == 'bwa') {
    build_bwa_index_ref(ch_reference_fasta)
    bwa_align_ref(extract_unmapped_human.out.ch_human_bam_unmapped_reads, ch_reference_fasta.collect(), build_bwa_index_ref.out.collect(), ch_where_ref_bwa.collect())
    sort_bam_ref(bwa_align_ref.out.ch_ref_bams)
    extract_unmapped_ref(sort_bam_ref.out.ch_ref_sorted_bams)
}

/****************************************************************
*********** build kraken2 database(s) and classify reads ********
****************************************************************/

build_krakendb(ch_krn2_db)
ch_unmapped_reads = extract_unmapped_human.out.ch_human_bam_unmapped_reads
classify_reads(ch_unmapped_reads, build_krakendb.out.collect())


/******************************************
********* summarize read counts ***********
******************************************/
count_raw_reads(ch_raw_reads)
ch_raw_counts = count_raw_reads.out.ch_read_counts.map { tuple -> [file(tuple[1], checkIfExists: true)] }.collect()
ch_raw_counts = merge_raw_counts(ch_raw_counts)

count_trimmed_reads(ch_trimmed_reads_trimmomatic)
ch_trimmed_counts = count_trimmed_reads.out.ch_read_counts.map { tuple -> [file(tuple[1], checkIfExists: true)] }.collect()
ch_trimmed_counts = merge_trimmed_counts(ch_trimmed_counts)

count_filtered_phix_reads(extract_unmapped_phix.out.ch_phix_bam_unmapped_reads)
ch_filtered_phix_counts = count_filtered_phix_reads.out.ch_read_counts.map { tuple -> [file(tuple[1], checkIfExists: true)] }.collect()
ch_filtered_phix = merge_filtered_phix_counts(ch_filtered_phix_counts)

ch_counts = merge_raw_counts.out.ch_merged_counts.mix(merge_trimmed_counts.out.ch_merged_counts, merge_filtered_phix_counts.out.ch_merged_counts).collect()

merge_all_counts(ch_counts)

/**************************************************************************
************ parse the report outputs and merge in mpa version ************
**************************************************************************/
ch_reports_text_files = classify_reads.out.ch_kraken_report.map { tuple -> [file(tuple[1], checkIfExists: true)] }.collect()
merge_kraken_reports(ch_reports_text_files)
visualize_species_summary(merge_all_counts.out.ch_readcounts, merge_kraken_reports.out.ch_merged_species_summary)



// /**********************************
// ********** multiqc report ********
// **********************************/

multiqc(ch_multiqc_config, 
    fastqc.out.ch_fastqc_results.collect().ifEmpty([]), 
    //trimmomatic.out.ch_trimmed_trimmo_stats.collect(),
    create_workflow_summary(summary))


}


/***************** 
***** output *****
*****************/

c_green = "\033[0;32m";
c_reset = "\033[0m";
c_grey = "\u001B[1;30m";
workflow.onComplete {
log.info """
  ______________________________________________________________

  Execution status: ${ workflow.success ? 'OK' : 'failed' }
  ${c_green}Results are reported here: $params.outdir${c_reset}
  Please cite:
  ______________________________________________________________
""".stripIndent()
}

/**************************
**** check hostname(s) ****
**************************/

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}



/*************
*    help    *
*************/

def helpMessage() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_light_blue = "\u001B[1:34m";
    c_dim = "\033[2m";
    log.info """
    -----------------------------------------------------------------------------------------------------------------------------
    ${c_light_blue}ViClaRA${c_reset} -- ${c_light_blue}Vi${c_reset}rus ${c_light_blue}Cl${c_reset}assification and ${c_light_blue}R${c_reset}eference-based ${c_light_blue}A${c_reset}ssembly

    ${c_yellow}Usage example: ${c_reset}
    nextflow run viclara.nf --reads "testdataset/*_R{1,2}_001.fastq.gz" \\
        -profile conda \\ 
        --removePhiX \\
        --aligner bwa \\
        --adapters adapters.fa \\
        --hg GCF_000001405.39_GRCh38.p13_genomic.fasta \\
        --reference GCF_000002315.6_GRCg6a_genomic.fna \\
        --outdir viclara_output \\
        --work-dir /tmp
    -----------------------------------------------------------------------------------------------------------------------------

    ${c_yellow}Mandatory arguments:${c_reset}
    ${c_green} --reads${c_reset}                        Path to input data (reads) in FASTQ format (must be surrounded with quotes)
    ${c_green} -profile${c_reset}                       Configuration profile to use. Available: conda, docker, singularity, slurm, test

    ${c_yellow}Generic:${c_reset}
    --singleEnd                     Specifies that the input is single-end reads

    ${c_yellow}Filter off-target reads:${c_reset}
    --removePhiX                    Remove PhiX sequences (uses the illumina PhiX reference by default)
    --hg                            Path to the human reference genome file in FASTA format (optional)
    --reference                     Path to the host reference genome file in FASTA format
    --saveReference                 Save the reference file in the results directory [default: $params.saveReference]
    --saveGenomeIndex               If generated by the pipeline save the BWA index in the results directory

    ${c_yellow}Trimming:${c_reset}
    --trimmer                       Trimming tool to use, available options (trimmomatic, bbduk) [default: $params.trimmer]
    --adapters                      Path to the adapters fasta
    --encoding                      Output quality offset: (phred33 (Sanger), phred64) [default: $params.encoding]
    --leading                       Instructs Trimmomatic to cut bases off the start of a read, if below a threshold quality [default: $params.leading]
    --trailing                      Instructs Trimmomatic to cut bases off the end of a read, if below a threshold quality [default: $params.trailing]
    --average_qual                  Instructs Trimmomatic the average quality required in the sliding window [default: $params.average_qual]
    --min_length                    Instructs Trimmomatic to drop the read if it is below a specified length [default: $params.min_length]
    --saveTrimmed                   Save the trimmed FASTQ files from the trimming step [default: $params.saveTrimmed]

    ${c_yellow}Alignment:${c_reset}
    --aligner                       Alignment/mapping tool to use, available options (bwa, bowtie2) [default: $params.aligner]
    --saveAlignedIntermediates      Save the BAM files from the aligment step [default: $params.saveAlignedIntermediates]

    ${c_yellow}Classification:${c_reset}
    --classify                      Perform classification using Kraken2 [default: $params.classify]
    --krn2_task                     Operation to be performed by kraken2, options ['download-taxonomy', 'download-library', 'build', 'standard']
    --krn2_db                       Path to the kraken2 database name
    --krn2_library                  Library to download sequences from, options: ['archaea', 'bacteria', 'plasmid','viral', 'human', 'fungi', 'plant', 'protozoa','nr', 'nt', 'env_nr', 'env_nt', 'UniVec','UniVec_Core']
    
    ${c_yellow}QC:${c_reset}
    --skipQC                        Skip all QC steps apart from MultiQC [default: $params.skipQC]
    --skipMultiQC                   Skip MultiQC [default: $params.skipMultiQC]
    --skipFastQC                    Skip FastQC [default: $params.skipFastQC]
    --skipTrimming                  Skip the trimming step [default: $params.skipTrimming]

    ${c_yellow}General options:${c_reset}
    --cores                         cores per process for local use [default: $params.cores]
    --max_cores                     max cores used on the machine for local use [default: $params.max_cores]
    --memory                        memory limit (GB) for alignment [default: $params.memory]
    --outdir                        The output directory where the results will be saved [default: $params.outdir]
    -w/--work-dir                   The temporary directory where intermediate data will be saved
    -name                           Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic


      

    """.stripIndent()
}