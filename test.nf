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

    common_chunk_frequency = Channel.fromPath(params.phase_common_list)
        .map { csv -> csv.splitCsv(header: false , sep: '\t') // Applies splitCsv to a file // Note: this is not the channel operator, but has same args
            .groupBy{ it[1] } // Group by chr // Groovy groupBy Collection function
        }
        .flatMap { map_grouped_by_chr -> map_grouped_by_chr
            .collect{ key, value -> value // key = chr, value = list of entries for chr // Groovy collect Collection function
                //.collect { it + [ value.size() ] } // Appends number of entries in chr to end of list
                .collect { [ value.size() ] + it } // Appends number of entries in chr to start of list
            } // Makes a list of lists
        } // Emits each outer list
        .flatMap() // Emits each inner list from each outer list
        //.view()
        .map { row -> tuple(row[2], row[0]) }
        .unique()//.view()
        //.toList().view()
        //.flatten().view()
        //.collate(2).collectEntries().view()
        //.flatten().view().toSpreadMap().view()
        //common_chunk_frequency.view()

    //common_chunk_frequency.collate(2).collectEntries().view()
    //common_chunk_frequency.toList().view().flatten().view().toSpreadMap().view
    //common_chunk_frequency.flatten().toSpreadMap()


    rare_chunk_frequency = Channel.fromPath(params.phase_rare_list)
        .map { csv -> csv.splitCsv(header: false , sep: '\t') // Applies splitCsv to a file // Note: this is not the channel operator, but has same args
            .groupBy{ it[1] } // Group by chr // Groovy groupBy Collection function
        }
        .flatMap { map_grouped_by_chr -> map_grouped_by_chr
            .collect{ key, value -> value // key = chr, value = list of entries for chr // Groovy collect Collection function
                //.collect { it + [ value.size() ] } // Appends number of entries in chr to end of list
                .collect { [ value.size() ] + it } // Appends number of entries in chr to start of list
            } // Makes a list of lists
        } // Emits each outer list
        .flatMap() // Emits each inner list from each outer list
        //.view()
        .map { row -> tuple(row[2], row[0]) }
        .unique()//.view()
        //rare_chunk_frequency.view()


    Channel
        .fromPath(params.phase_common_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[1], row[0], row[2]) }
        //.view()
        .set { phase_common_chunks_ch }

    Channel
        .fromPath(params.phase_rare_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[1], row[0], row[2], row[3]) }
        //.view()
        .set { phase_rare_chunks_ch }


    Channel
        .fromPath(params.input_files_list)
        .splitCsv(header: false, sep: '\t')
        .map { row -> tuple(row[0], file(row[1]), file(row[1]) + '.csi', file(row[2])) }
        //.view()
        .set { input_files_ch }

    phase_common_chunks_ch.combine(input_files_ch,by:0).set{commoncombinedset}
    //commoncombinedset.view()


    testcommon(commoncombinedset)
    //testcommon.out.view()


    //chromosomes_list.view()
    //chromosomes_list.combine(testcommon.out,by:0).map { chr_no, bcf, idx -> [chr_no, bcf, idx] }.view().set{ligate}

    //common_chunk_frequency_1 = [ "chr21": 3, "chr22": 3 ]


    //common_chunks_size = phase_common_chunks_chr_grpkey.map{ chr, freq -> [ chr ] }.toList()
    //common_chunks_size = phase_common_chunks_chr_grpkey.map{ chr, freq -> [ chr ] }.toList()
            // Adding number of chunks as elements
    //         .map{ chr -> [ chr, chr.size() ] }
    //         .transpose()
    //         .view()

    //common_chunks_size1 = input_files_ch.map{ chr, file1, file2, file3, freq -> [ chr, freq ] }.toList()
    //         // Adding number of chunks as elements
    //         .map{ it -> [ it, it.size() ] }
    //         .transpose()
    //         .view()

    common_chunk_frequency.combine(testcommon.out,by:0).map{chr, count, bcf -> tuple( groupKey(chr, count), bcf ) }.groupTuple( by: 0 ).map {chr, bcf -> tuple(chr, bcf.flatten() ) }.set{ligate}
    //common_chunk_frequency.combine(testcommon.out,by:0).map{chr, count, bcf -> tuple( groupKey(chr, count), bcf.flatten() ) }.groupTuple( by: 0 ).map {chr, bcf -> tuple(chr, bcf.flatten() ) }.view().set{ligate}
    //common_chunk_frequency.combine(testcommon.out,by:0).view()

    ////testcommon.out.map { chr, bcf -> tuple(chr, bcf.flatten() ) }.groupTuple( by: 0 ).map {chr, bcf -> tuple(chr, bcf.flatten() ) }.view().set{ligate}
    /////////testcommon.out.groupTuple( by: 0 ).map {chr, bcf -> tuple(chr, bcf.flatten() ) }.view().set{ligate}
    ///////// testcommon.out.map { chr, bcf -> tuple( groupKey(chr, common_chunk_frequency_1[chr]), bcf.flatten() )  }.groupTuple( by: 0 ).map {chr, bcf -> tuple(chr, bcf.flatten() ) }.view().set{ligate}


    //testcommon.out.map { chr, bcf -> [ groupKey(chr, common_chunk_frequency[chr]), bcf.flatten() ]  }.groupTuple( by: 0 ).view().set{ligate}
    //testcommon.out.map { chr, bcf -> tuple( groupKey(chr, chr.size()), bcf.flatten() )  }.view().groupTuple( by: 0 ).set{ligate}

    //testcommon.out.map { chr, bcf -> [ groupKey(chr, common_chunk_frequency[chr]), bcf.flatten() ]  }.groupTuple( by: 0 ).view().set{ligate}
    //testcommon.out.map { chr, bcf -> tuple( groupKey(chr, common_chunk_frequency[chr]), bcf.flatten() ) }.groupTuple( by: 0 ).view().set{ligate}
    //testcommon.out.map{ chr, bcf -> [chr, bcf.flatten() ] }.groupTuple( by: 0 ).view().set{ligate}


    //testcommon.out.groupTuple( by: 0 ).map { chr, bcf -> tuple( groupKey(chr, phase_common_chunks_chr_grpkey[chr]), bcf.flatten() ) }.set{ligate}
    //testcommon.out.map { chr, bcf -> tuple( groupKey(chr, phase_common_chunks_chr_grpkey[chr]), bcf.flatten() ) }.groupTuple( by: 0 ).set{ligate}

    testligate(ligate)

    //phase_rare_chunks_ch.combine(input_files_ch,by:0).set{rarecombinedset}
    //rarecombinedset.combine(testligate.out,by:0).flatten().view()
    phase_rare_chunks_ch.combine(input_files_ch,by:0).combine(testligate.out,by:0).map { chr, chunk_no, scaffold_reg, input_reg, bcf, gmap, scaffold, scaffold_idx -> tuple(chr, chunk_no, scaffold_reg, input_reg, bcf, gmap, scaffold, scaffold_idx) }.view().set{phaserare}

    testrare(phaserare)

    rare_chunk_frequency.combine(testrare.out,by:0).map{chr, count, bcf -> tuple( groupKey(chr, count), bcf ) }.groupTuple( by: 0 ).map {chr, bcf -> tuple(chr, bcf.flatten() ) }.view().set{concat}

    testconcat(concat)

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
