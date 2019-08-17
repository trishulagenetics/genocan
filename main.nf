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

		--reads			                  Path to input data (must be surrounded with quotes)
		-profile		                  Configuration profile to use [docker/awsbatch]

	Options:
		--singleEnd		                Specifies that the input is single end reads
		--outdir		                  Specify output directory where the results will be saved
		--email		 	                  Specify email address to send run details to
		-name			                    Name for the pipeline run. If not specified, Nextflow will generate a random name

	References			                If not specified in the configuration file or if you wish to overwrite any of the references
		--fasta			                  Path to fasta reference
    --bwa_index                   Path to bwa index files
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
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.bwa_index = './bwa_index'
params.fasta_index = false
params.saveReference = false
params.snpcapture = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")

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

if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Validate inputs

Channel.fromPath("${params.fasta}")
    .ifEmpty { exit 1, "No genome specified! Please specify one with --fasta or --bwa_index"}
    .into {ch_fasta_for_bwa_indexing; ch_fasta_for_faidx_indexing; ch_fasta_for_variant_call; ch_fasta_for_bwamem_mapping; ch_fasta_for_qualimap}

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
summary['bwa index'] = params.bwa_index
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

def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-eager-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/eager Workflow Summary'
    section_href: 'https://github.com/nf-core/eager'
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
    bwa &> v_bwa.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    python ${baseDir}/bin/scrape_software_versions.py > software_versions_mqc.yaml
    """
}

process build_bwa_index {
    tag {fasta}

    publishDir path: "${params.outdir}/bwa_index", mode: 'copy', saveAs: { filename ->
            if (params.saveReference) filename
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    when: !params.bwa_index && params.fasta

    input:
        
    file fasta from ch_fasta_for_bwa_indexing
    file wherearemyfiles

    output:
        
    file "*.{amb,ann,bwt,pac,sa,fasta,fa}" into bwa_index
    file "where_are_my_files.txt"

    """
    bwa index $fasta
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
    file fasta from ch_fasta_for_faidx_indexing
    file wherearemyfiles

    output:
    file "${fasta}.fai" into fasta_index
    file "${fasta}"
    file "where_are_my_files.txt"

    """
    samtools faidx $fasta
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

process bwa_align {
    tag "$name"
    
    publishDir "${params.outdir}/mapping/bwamem", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_fastq
    file fasta from ch_fasta_for_bwamem_mapping
    file index from bwa_index
            
    output:
    file "*_sorted.bam" into bwa_sorted_bam_idxstats, bwa_sorted_bam_filter
    file "*.bai"
            
    script:

    if(params.singleEnd){
    """ 
    bwa mem bwa_index ${reads[0]} -t ${task.cpus} | samtools sort -@ ${task.cpus} -o ${name}_sorted.bam
    samtools index -@ ${task.cpus} ${name}_sorted.bam
    """ 
    } else {
    """ 
    bwa mem bwa_index ${reads[0]} ${reads[1]} -t ${task.cpus} | samtools sort -@ ${task.cpus} -o ${name}_sorted.bam
    samtools index -@ ${task.cpus} ${name}_sorted.bam
    """ 
    }

}

process samtools_idxstats {
    tag "$prefix"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'

    input:
    file(bam) from bwa_sorted_bam_idxstats

    output:
    file "*.stats" into bwamem_idxstats_for_multiqc

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
    file bam from bwa_sorted_bam_filter
    
    output:
    file "*filtered.bam" into bam_filtered_qualimap, bam_filtered_call_variants
    file "*.fastq.gz"
    file "*.bai"

    script:
    
    prefix = "$bam" - ~/(\.bam)?$/

    """
    samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q 15 -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -F4 -q 15 -o ${prefix}.filtered.bam)
    samtools fastq -tn "${prefix}.unmapped.bam" | gzip > "${prefix}.unmapped.fastq.gz"
    samtools index -@ ${task.cpus} ${prefix}.filtered.bam
    """
}

process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    when:
    !params.skip_qualimap

    input:
    file bam from bam_filtered_qualimap
    file fasta from ch_fasta_for_qualimap

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
    file fasta from ch_fasta_for_variant_call
    file fasta_index

    output:
    file "${prefix}_variants.vcf" into intial_vcf
    file "${prefix}_variants_filtered.vcf" into filtered_vcf

    script:
    
    prefix = "$bam" - ~/(\_sorted.filtered.bam)?$/

    """
    bcftools mpileup -Ou -f ${fasta} ${bam} | bcftools call --multiallelic-caller --variants-only --no-version --threads ${task.cpus} -Oz -o ${prefix}_variants.vcf
    bcftools filter -i 'DP>10' ${prefix}_variants.vcf > ${prefix}_variants_filtered.vcf
    """
}

process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('samtools/*') from bwamem_idxstats_for_multiqc.collect().ifEmpty([])
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
    Rscript ${baseDir}/bin/markdown_to_html.R $output_docs results_description.html
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