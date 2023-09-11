#!/bin/bash

# run the workflow on SH10K phased data with below command.

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../count_variants.nf \
	-config ../nextflow.config \
	-params-file params.count.yml \
	-work-dir ./work \
	--outdir ./output \
	-profile local \
	-resume
