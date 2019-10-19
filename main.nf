#!/usr/bin/env nextflow
/*
========================================================================================
                         trishulagenetics/genocan
========================================================================================
 trishulagenomics/genocan Analysis Pipeline. Started 19/07/2019.
 #### Homepage / Documentation
 https://github.com/trishulagenomics/genocan
 #### Authors
 Andries van Tonder <andries@trishulagenomics.com> - https://github.com/trishulagenomics>
----------------------------------------------------------------------------------------
*/

def helpMessage() {
	log.info"""
	=============================================
     trishulagenomics/genocan v${params.version}
	=============================================
	Usage:

	The typical command for running the pipeline is as follows:

	nextflow run trishulagenomics/genocan --reads '*_R{1,2}.fastq.gz' -profile docker

	Mandatory arguments:

	/*
	 * NEED TO DECIDE WHETHER INPUT IS WHOLE LANE OR PRE-DEMULTIPLEXED!
	 */

		--reads			                Path to input data (must be surrounded with quotes)
		-profile		                Configuration profile to use [docker/awsbatch]

	Options:
		--singleEnd		                Specifies that the input is single end reads
		--outdir		                Specify output directory where the results will be saved
		-name			                Name for the pipeline run. If not specified, Nextflow will generate a random name

	References			          If not specified in the configuration file or if you wish to overwrite any of the references
	--fasta			              Path to fasta reference
    --minimap2_index              Path to minimap2 index file
    --fasta_index                 Path to fasta index
    --saveReference               Saves reference genome indices for later reusage
    --snpcapture                  Runs in SNPCapture mode (specify a BED file if you do this!)

	""".stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message

params.help = false
if (params.help){
	helpMessage()
	exit 0
}

// Configurable variables

params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.minimap2_index = false
params.fasta_index = false
params.saveReference = false
params.snpcapture = false

multiqc_config = file(params.multiqc_config)
wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")

ref_fasta = file(params.fasta)

// Validate inputs

if (params.fasta) {
  fasta = file(params.fasta)
  if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}

// Use user-specified run name

custom_runName = params.name

if ( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ) {
	custom_runName = workflow.runName
}

// Check that workDir/outdir paths are s3 buckets if using AWSbatch

//if(workflow.profile == 'awsbatch'){
  //  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    //if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
//}

// Validate inputs

/*
 * Create a channel for input read files
 */

params.singleEnd = false

Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trim }

// Header log info

log.info """
============================================
 trishulagenomics/genocan v${params.version}
============================================
"""

def summary = [:]
summary['Pipeline name'] =  'trishulagenomics/genocan'
summary['Pipeline version'] = params.version
summary['Run name'] = custom_runName ?: workflow.runName
summary['Reads'] = params.reads
summary['Fasta reference'] = params.fasta
summary['minimap2 index'] = params.minimap2_index
summary['Data type'] = params.singleEnd ? 'Single-end' : 'Paired-end'
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'trishulagenomics-genocan-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'trishulagenomics/genocan Workflow Summary'
    section_href: 'https://github.com/trishulagenomics/genocan'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

// Check that Nextflow version is up to date enough

try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version out of date')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue but may not work like you wish.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * Parse software version numbers
 */

process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version &> v_fastqc.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    fastp --version &> v_fastp.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    minimap2 --version &> v_minimap2.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    python s3://trishulagenomics/scripts/scrape_software_versions.py > software_versions_mqc.yaml
    """
}

process build_minimap2_index {
    tag {fasta}

    publishDir path: "${params.outdir}/minimap2_index", mode: 'copy', saveAs: { filename ->
            if (params.saveReference) filename
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    when: !params.minimap2_index && params.fasta

    input:
        
    file ref_fasta
    file wherearemyfiles

    output:
        
    file "*.mmi" into minimap2_index_minimap2
    file "where_are_my_files.txt"

    """
    minimap2 -d "${ref_fasta}.mmi" "${ref_fasta}"
    """
}

process build_fasta_index {
    tag {fasta}

    publishDir path: "${params.outdir}/fasta_index", mode: 'copy', saveAs: { filename ->
            if (params.saveReference) filename
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    when: !params.fasta_index && params.fasta

    input:
    file ref_fasta
    file wherearemyfiles

    output:
    file "${ref_fasta}.fai" into fasta_index
    file "where_are_my_files.txt"

    """
    samtools faidx "${ref_fasta}"
    """
}

process fastqc {
    tag "$name"
    
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

process fastp {
    tag "$name"
    
    publishDir "${params.outdir}/FastP", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_trim

    output:
    set val(name), file("*trim.fastq.gz") into trimmed_fastq
    file "*.json" into fastp_for_multiqc

    script:
    
    if(params.singleEnd){
    
    """
    fastp --in1 ${reads[0]} --out1 "${reads[0].baseName}_trim.fastq.gz" -w ${task.cpus} --json "${reads[0].baseName}"_fastp.json
    """
    
    } else {
    
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${reads[0].baseName}_trim.fastq.gz" --out2 "${reads[1].baseName}_trim.fastq.gz" -w ${task.cpus} --json "${reads[0].baseName}"_fastp.json
    """
    
    }
}

process minimap2_align {
    tag "$name"
    
    publishDir "${params.outdir}/mapping/minimap2", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_fastq
    file minimap2_index from minimap2_index_minimap2
            
    output:
    file "*_sorted.bam" into minimap2_sorted_bam_idxstats, minimap2_sorted_bam_filter
    file "*.bai"
            
    script:

    if(params.singleEnd){
    
    """ 
    minimap2 -ax sr $minimap2_index ${reads[0]} -t ${task.cpus} > ${name}.sam
    samtools view -S -b ${name}.sam > ${name}.bam
    samtools sort ${name}.bam -o ${name}_sorted.bam
    samtools index ${name}_sorted.bam
    """ 
    
    } else {
    
    """ 
    minimap2 -ax sr $minimap2_index ${reads[0]} ${reads[1]} -t ${task.cpus} > ${name}.sam
    samtools view -S -b ${name}.sam > ${name}.bam
    samtools sort ${name}.bam -o ${name}_sorted.bam
    samtools index ${name}_sorted.bam
    """ 
    
    }
}

process samtools_idxstats {
    tag "$prefix"
    
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'

    input:
    file bam from minimap2_sorted_bam_idxstats

    output:
    file "*.stats" into minimap2_idxstats_for_multiqc

    script:
    
    prefix = "$bam" - ~/(\.bam)?$/

    """
    samtools flagstat $bam > ${prefix}.stats
    """
}

process samtools_filter {
    tag "$prefix"
    
    publishDir "${params.outdir}/samtools/filter", mode: 'copy'

    input:
    file bam from minimap2_sorted_bam_filter
    
    output:
    file "*filtered.bam" into bam_filtered_qualimap, bam_filtered_call_variants
    file "*.bai"

    script:
    
    prefix = "$bam" - ~/(\.bam)?$/

    """
    samtools view -h $bam -f4 -q 15 -b -o ${prefix}.filtered.bam
    samtools index ${prefix}.filtered.bam
    """
}

process qualimap {
    tag "${bam.baseName}"
    
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    when:
    !params.skip_qualimap

    input:
    file bam from bam_filtered_qualimap
    file ref_fasta

    output:
    file "*" into qualimap_results

    script:
    
    snpcap = ''
    
    if(params.snpcapture) snpcap = "-gff ${params.bedfile}"
    
    """
    qualimap bamqc -bam $bam -nt ${task.cpus} -outdir . -outformat "HTML" ${snpcap}
    """
}

process variant_call {
    tag "$prefix"
    
    publishDir "${params.outdir}/samtools/variants", mode: 'copy'

    input:
    file bam from bam_filtered_call_variants
    file ref_fasta
    file fasta_index

    output:
    file "${prefix}_variants.vcf" into intial_vcf
    file "${prefix}_variants_filtered.vcf" into filtered_vcf

    script:
    
    prefix = "$bam" - ~/(\_sorted.filtered.bam)?$/

    """
    bcftools mpileup -Oz -f "${ref_fasta}" ${bam} -o ${prefix}.vcf
    bcftools call --multiallelic-caller --variants-only --no-version --threads ${task.cpus} -Oz ${prefix}.vcf -o ${prefix}_variants.vcf
    bcftools filter -i 'DP>10' ${prefix}_variants.vcf > ${prefix}_variants_filtered.vcf
    """
}

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('samtools/*') from minimap2_idxstats_for_multiqc.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect().ifEmpty([])
    file ('qualimap/*') from qualimap_results.collect().ifEmpty([])
    file ('fastp/*') from fastp_for_multiqc.collect().ifEmpty([])

    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}