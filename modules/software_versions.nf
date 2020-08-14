/*****************************************************************
****************** parse software version numbers ****************
*****************************************************************/

process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename -> 
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
        file 'software_versions_mqc.yaml'
        file "software_versions.csv"

    script:
        """
        echo $workflow.manifest.version &> v_ngi_viclara.txt
        echo $workflow.nextflow.version &> v_nextflow.txt
        fastqc --version &> v_fastqc.txt
        trimmomatic -version &> v_trimmomatic.txt
        samtools --version &> v_samtools.txt
        multiqc --version &> v_multiqc.txt
        echo `bwa 2>&1 | grep Version | cut -d " " -f2 | cut -f1 -d "-"` &> v_bwa.txt
        scrape_software_versions.py &> software_versions_mqc.yaml
        """
}