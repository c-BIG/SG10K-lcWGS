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

include { phasingcommon } from './modules/phasingcommon'
include { ligatecommon } from './modules/ligatecommon'

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/

// main

workflow {

    vcf = file( params.vcf )
    vcf_idx = file( params.vcf + ".csi" )
    gmap = file( params.gmap )

    Channel
        .fromPath(params.phase_common_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[1], row[0], row[2]) }
        .set { phase_common_run_ch }

    phasingcommon( phase_common_run_ch, vcf, vcf_idx, gmap )

    Channel
        //.fromList( chr1..chr22 )
        .of('ch22')
        .set {chromosomes_list}

    chromosomes_list.combine(phasingcommon.out,by:0)
        .map { chr_no, bcf, idx -> [chr_no, bcf, idx] }
        .set { ligate }

    //ligatecommon(chromosomes_list.combine(phasingcommon.out,by:0))
    //chromosomes_list.combine(phasingcommon,by:0).set{ligate}
    //chromosomes_list.combine(phasingcommon,by:0).view()
    ligatecommon(ligate)
    //ligatecommon( phasingcommon.out.collect() )

}

//         .map { row-> tuple(row.chr_no, row.chunk_region1, row.chunk_region2, file(row.vcf), file(row.index), file(row.gmap)) }

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
