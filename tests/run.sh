#!/bin/bash

# run the workflow on SH10K phased data with below command.

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../prepare_ref_panel.nf \
	-config ../nextflow.config \
	-params-file params.yml \
	-work-dir ./work \
	--publish_dir ./output \
	-profile local 
