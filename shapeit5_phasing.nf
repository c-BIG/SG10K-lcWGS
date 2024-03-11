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

include { testcommon } from './modules/testcommon'
include { testligate } from './modules/testligate'
include { testrare } from './modules/testrare'
include { testconcat } from './modules/testconcat'

/*
----------------------------------------------------------------------
WORKFLOW
---------------------------------------------------------------------
*/

// main

workflow {

    PREFIX = params.prefix

    // get the count of common chunks per chromosomes, to pass groupkey to initiate the ligate/concat process immediate to the no. of chunks processed by phase common
    common_chunk_frequency = Channel.fromPath(params.phase_common_list)
        .map { csv -> csv.splitCsv(header: false , sep: '\t') // Applies splitCsv to a file // Note: this is not the channel operator, but has same args
            .groupBy{ it[1] } // Group by chr // Groovy groupBy Collection function
        }
        .flatMap { map_grouped_by_chr -> map_grouped_by_chr
            .collect{ key, value -> value // key = chr, value = list of entries for chr // Groovy collect Collection function
                .collect { [ value.size() ] + it } // Appends number of entries in chr to start of list
            } // Makes a list of lists
        } // Emits each outer list
        .flatMap() // Emits each inner list from each outer list
        .map { row -> tuple(row[2], row[0]) }
        .unique()//.view()


    // get the count of rare chunks per chromosomes, to pass groupkey to initiate the ligate/concat process immediate to the no. of chunks processed by phase rare
    rare_chunk_frequency = Channel.fromPath(params.phase_rare_list)
        .map { csv -> csv.splitCsv(header: false , sep: '\t') // Applies splitCsv to a file // Note: this is not the channel operator, but has same args
            .groupBy{ it[1] } // Group by chr // Groovy groupBy Collection function
        }
        .flatMap { map_grouped_by_chr -> map_grouped_by_chr
            .collect{ key, value -> value // key = chr, value = list of entries for chr // Groovy collect Collection function
                .collect { [ value.size() ] + it } // Appends number of entries in chr to start of list
            } // Makes a list of lists
        } // Emits each outer list
        .flatMap() // Emits each inner list from each outer list
        .map { row -> tuple(row[2], row[0]) }
        .unique()//.view()

    // channel to prepare a list from common chunks 20cM "chr, chunk order, chunk regions"
    Channel
        .fromPath(params.phase_common_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[1], row[0], row[2]) }
        .set { phase_common_chunks_ch }

    // channel to prepare a list from rare chunks 4cM "chr, chunk order, chunk regions"
    Channel
        .fromPath(params.phase_rare_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[1], row[0], row[2], row[3]) }
        .set { phase_rare_chunks_ch }

    // make a tuple of unrelated samples chr and its bcf, bcf index files path and gmap
    Channel
        .fromPath(params.input_files_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), file(row[1]) + '.csi', file(row[2])) }
        .set { input_files_ch }

    //combine the common chunk regions and samples bcf and gmap files list
    phase_common_chunks_ch
        .combine(input_files_ch,by:0)
        .set{commoncombinedset}

    //run phase common varinats per chunk region
    testcommon(commoncombinedset)

    //collect phase common output and combine the count list of common chunk regions per chr to ligate
    common_chunk_frequency
        .combine(testcommon.out.common_bcf_files,by:0)
        .map{chr, count, bcf -> tuple( groupKey(chr, count), bcf ) }
        .groupTuple( by: 0 )
        .map {chr, bcf -> tuple(chr, bcf.flatten() ) }
        .set{ligate}

    //run ligate of phase common outputs
    testligate(ligate)

    // combine the rare 4cM chunks region, samples bcf and its index and phase common ligated scaffold and its index per chr
    phase_rare_chunks_ch
        .combine(input_files_ch,by:0)
        .combine(testligate.out.common_bcf_files,by:0)
        .map { chr, chunk_no, scaffold_reg, input_reg, bcf, gmap, scaffold, scaffold_idx -> tuple(chr, chunk_no, scaffold_reg, input_reg, bcf, gmap, scaffold, scaffold_idx) }
        .view().set{phaserare}

    //run phase rare varinats per chunk region
    testrare(phaserare)

    // collect phase rare output and combine the count list of chunk regions per chr to concat
    rare_chunk_frequency
        .combine(testrare.out.rare_scaffold_files,by:0)
        .map{chr, count, bcf -> tuple( groupKey(chr, count), bcf ) }
        .groupTuple( by: 0 )
        .map {chr, bcf -> tuple(chr, bcf.flatten() ) }
        .view().set{concat}

    //run concate process of all the phase rare outputs
    testconcat(concat)

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
