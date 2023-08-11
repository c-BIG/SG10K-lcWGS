#!/bin/bash

# run the workflow on SH10K phased data with below command.

export CAPSULE_LOG=none
NXF_ANSI_LOG=false

nextflow run ../main.nf \
	-config ../nextflow.config \
	--input_list samplesheet.csv \
	--ref reference.csv \
	--outdir output \
	-profile local
