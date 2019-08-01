#!/usr/bin/env nextflow
/*
========================================================================================
                         trishulagenetics/genocan
========================================================================================
 trishulagenetics/genocan Analysis Pipeline. Started 19/07/2019.
 #### Homepage / Documentation
 https://github.com/trishulagenetics/genocan
 #### Authors
 Andries van Tonder <trishulagenetics@gmail.com> - https://github.com/trishulagenetics>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
	log.info"""
	=============================================
     trishulagenetics/genocan v${params.version}
	=============================================
	Usage:

	The typical command for running the pipeline is as follows:

	nextflow run trishulagenetics/genocan --reads '*_R{1,2}.fastq.gz' -profile docker

	Mandatory arguments:

	/*
	 * NEED TO DECIDE WHETHER INPUT IS WHOLE LANE OR PRE-DEMULTIPLEXED!
	 */

		--reads			Path to input data (must be surrounded with quotes)
		-profile		Configuration profile to use [docker/awsbatch]

	Options:
		--singleEnd		Specifies that the input is single end reads
		--outdir		Specify output directory where the results will be saved
		--email			Specify email address to send run details to
		-name			Name for the pipeline run. If not specified, Nextflow will generate a random name

	References			If not specified in the configuration file or if you wish to overwrite any of the references
		--fasta			Path to fasta reference

	""".stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message

if (params.help){
	helpMessage()
	exit 0
}

// Configurable variables

params.name = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)

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

if (workflow.profile == 'awsbatch') {
	if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

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
 trishulagenetics/genocan v${params.version}
============================================
"""

def summary = [:]
summary['Pipeline name'] =  'trishulagenetics/genocan'
summary['Pipeline version'] = params.version
summary['Run name'] = custom_runName ?: workflow.runName
summary['Reads'] = params.reads
summary['Fasta reference'] = params.fasta
summary['Data type'] = params.singleEnd ? 'Single-end' : 'Paired-end'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
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
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

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
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    trimmomatic -version > v_trimmomatic.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """

}

process trimmomatic {
    tag "$name"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_trim

    output:
    set val(name), file("trim_*.fastq.gz") into trimmed_fastq
    file "*_trim.out" into trim_logs

    script:
    trunc_string = "MINLEN:30 CROP:30"
    if(params.trim_truncate > 30){
      trunc_string = "MINLEN:${params.trim_truncate} CROP:${params.trim_truncate}"
    }
    
    if (params.singleEnd) {
      """
      trimmomatic SE \
      -threads ${task.cpus} \
      -trimlog ${name}_trim.log \
      -phred33 \
      """
    } else {
      """
      trimmomatic PE \
      -threads ${task.cpus} \
      -trimlog ${name}_trim.log \
      -phred33 \
      ${reads} trim_${reads[0]} U_${reads[0]} \
      ILLUMINACLIP:${params.trim_adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 \
      ${trunc_string} 2> ${name}_trim.out
      """
      }
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

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('trimmomatic/*') from trim_logs.collect()
    file ('software_versions/*') from software_versions_yaml

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

/*
 * Create output results description
 */

 process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[trishulagenetics/genocan] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[trishulagenetics/genocan] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[trishulagenetics/genocan] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[trishulagenetics/genocan] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[trishulagenetics/genocan] Pipeline Complete"

}