#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Developed by the Genome Institute of Singapore for
the National Precision Medicine Programme

Copyright: 2022 Genome Institute of Singapore
License: The MIT License (MIT)

See LICENSE for more copyright information
*/

/*
----------------------------------------------------------------------
FUNCTIONS
----------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    Usage: nextflow run main.nf -config nextflow.config -params-file params.yaml 
                                -work-dir ./ --outdir ./
                                [-resume] [--keep_workdir] [--help]

    Options:
    -config           Generic workflow settings
    -params-file      Sample-specific settings
    -profile          Job launch settings
    -resume           Re-use existing results (optional, omit to re-run from scratch)
    --keep_workdir    Keep work directory (optional, omit for auto-deletion)
    --help            Print this help message
    """.stripIndent()
}

def nextflowMessage() {
    log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def version_message() {
    println("SG10K-lcWGS  ~  version ${workflow.manifest.version}")
}

def minimalInformationMessage() {
    log.info "User name    : " + workflow.userName
    log.info "Command Line : " + workflow.commandLine
    log.info "Project Dir  : " + workflow.projectDir
    log.info "Launch Dir   : " + workflow.launchDir
    log.info "Work Dir     : " + workflow.workDir
    log.info "Results Dir  : " + params.outdir
    log.info "Info Dir     : " + params.info_dir
    log.info "Profile      : " + workflow.profile
}

def startMessage() {
    this.nextflowMessage()
    this.version_message()
    this.minimalInformationMessage()
}


/*
----------------------------------------------------------------------
USAGE
----------------------------------------------------------------------
*/

if (params.help) exit 0, helpMessage()

/*
----------------------------------------------------------------------
LAUNCH INFO
----------------------------------------------------------------------
*/

startMessage()

/*
----------------------------------------------------------------------
PROCESSES
---------------------------------------------------------------------
*/

include { af_dist } from './modules/afdist'

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/

// main

workflow {


    Channel
        .fromList( params.samples )
        .ifEmpty { ['chr_no': params.chr_no, 'vcf_sg10k': params.vcf_sg10k, 'vcf_1kg': params.vcf_1kg] }
        .set { samples }


    Channel
        samples.branch { rec ->
            def vcf_sg10k_file = rec.vcf_sg10k ? file( rec.vcf_sg10k ) : null

            output: rec.chr_no && vcf_sg10k_file
                def vcf_sg10k_idx = file( "${rec.vcf_sg10k}.tbi" )
                def vcf_1kg_file = rec.vcf_1kg ? file( rec.vcf_1kg ) : null
                def vcf_1kg_idx = file( "${rec.vcf_1kg}.tbi" )

                return tuple( rec.chr_no, vcf_sg10k_file, vcf_sg10k_idx, vcf_1kg_file, vcf_1kg_idx )
        }
        .set { vcf_inputs }

    af_dist( vcf_inputs )


}

/*
----------------------------------------------------------------------
COMPLETION INFO
----------------------------------------------------------------------
*/
workflow.onComplete {
    log.info "Started     : " + workflow.start
    log.info "Completed   : " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Status      : " + workflow.success
    log.info "Publish dir : " + params.outdir
}

workflow.onError {
    log.info "Workflow execution stopped with the following message:"
    log.info "Exit status   : " + workflow.exitStatus
    log.info "Error message : " + workflow.errorMessage
    log.info "Error report  : " + (workflow.errorReport ?: '-')
}
