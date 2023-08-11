nextflow.enable.dsl=2
version = "2.0.0"

def help_message(String version) {
	println(
		"""
		==============================================
		Imputation Analysis (GLIMPSE)  ~  ${version}
		==============================================
		Mandatory arguments:
		--input_list	Path to file with list of sample name, sample depth and vcf files (one sample per line).
		--ref		Path to reference file with chunk & bin file (one chromosome per line).
		--outdir	Path of the output directory.
		
		Additional arguments:
		--help          Show this help message
		--version       Show the version of the workflow
		""".stripIndent()
	)
}

def version_message(String version) {
	println("Imputation Analysis (GLIMPSE)  ~  ${version}")
}

def help_or_version(Map params, String version){
	if (params.help) {
		help_message(version)
		System.exit(0)
	}
	if (params.version) {
		version_message(version)
		System.exit(0)
	}
}

def default_params(){
	def params = [:]
	params.help = false
	params.version = false
	params.input_list = false
	params.ref = false
	params.outdir = false
	return params
}

def check_params(Map params, nextflow.script.WorkflowMetadata workflow) {
	// merge defaults with user params
	final_params = default_params() + params

	// print help or version messages if requested
	help_or_version(final_params, version)

	// param checks
	check_mandatory_parameter(final_params, "input_list")
	check_mandatory_parameter(final_params, "ref")
	check_mandatory_parameter(final_params, "outdir")

	return final_params
}

def check_mandatory_parameter(Map params, String parameter_name){
	if ( !params[parameter_name] ){
		println "ERROR: Missing mandatory argument; specify '--help' for usage instructions"
		System.exit(1)
	} else {
		return params[parameter_name]
	}
}

final_params = check_params(params, workflow)


include { impute_phase } from './modules/impute_phase'
include { impute_ligate } from './modules/impute_ligate'
include { impute_concat } from './modules/impute_concat'

//Inputs
samples_ch = Channel
	.fromPath(final_params.input_list)
	.splitCsv(header:true)
	.map{ row -> tuple(row.run_id, row.sample_id, file(row.sample_cram), file(row.sample_cram + ".crai")) }

detail_ch = Channel
	.fromPath(final_params.ref)
	.splitCsv(header:true)
	.map{ row -> tuple(row.chr_no, file(row.chunk_txt), file(row.ref_bin)) }

ref_fasta = Channel
	.fromPath("s3://nalagenetics-npm-grids-collaboration/SG10K_imputation/reference/Homo_sapiens_assembly38.fasta")
fasta_index = Channel
	.fromPath("s3://nalagenetics-npm-grids-collaboration/SG10K_imputation/reference/Homo_sapiens_assembly38.fasta.fai")

phase_inputs = detail_ch.combine(samples_ch)


//Workflow
workflow{
phase_outputs = impute_phase(phase_inputs.combine(ref_fasta.combine(fasta_index)))
ligate_outputs = impute_ligate(phase_outputs)
impute_concat(ligate_outputs[0].groupTuple(by:[0,1]))
}

workflow.onComplete {

    println ( workflow.success ? """
        ==============================================
        Imputation Analysis (GLIMPSE) ~ EXECUTION SUMMARY
        ==============================================
        Workflow    : ${workflow.scriptName} ${version}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}
