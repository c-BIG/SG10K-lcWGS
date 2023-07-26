#Load python modules
import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import logging
from utils import basic_utils

#logger config
logging.basicConfig(format='%(asctime)s::%(levelname)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level = logging.DEBUG, handlers=[ logging.StreamHandler(sys.stdout)])
    

#load parameter from command line
def load_parameters():
    parser = argparse.ArgumentParser(
        description="Phasing & Imputation for 313 PRS Markers with GLIMPSE2"
    )
    parser.add_argument(
        "-S",
        "--sample",
        help="The given sample name",
        required=True,
    )
    parser.add_argument(
        "-C",
        "--chromosome",
        help="The specified chromosome",
        required=True,
    )
    parser.add_argument(
        "-I",
        "--input",
        help="The original input CRAM/BAM file containing genotype likelihoods",
        required=True,
    )
    parser.add_argument(
        "-R",
        "--reference",
        help="The path to folder reference panel of haplotypes in GLIMPSE2 bin format",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--chunk",
        help="The pre-processed chunk file",
        required=True,
    )
    parser.add_argument(
        "-F",
        "--fasta",
        help="Reference genome in Fasta format",
        required=True,
    )
    args = parser.parse_args()
    return args

def run_command(cmd, cmd_type="normal"):
    # run the command
    if cmd_type == "normal":
        cmds = cmd.split(" ")
        res = subprocess.run(cmds,capture_output=True)
        outp= res.stdout.decode('utf-8')
        error_msg = res.stderr.decode("utf-8").lower()
    elif cmd_type == "pipe":
        res = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
        outp = res.stdout.read().decode("utf-8")
        error_msg = res.stderr.read().decode("utf-8").lower()
        
    #error handling
    error_status = ["error", "no ", "unrecognized", "cannot", "not "]
    if (error_msg!='' and any(x in error_msg for x in error_status)):
        logging.error("CMD ERROR!")
        raise ValueError(f"Command: {cmd}\n Error Message: {error_msg}")
        
    return outp

#GLIMPSE phase function
def impute_phase(sample_id, chr_id, sample_cram, ref_bin, chunk_txt, fasta):
    
    #prepare input
    chunk_files = pd.read_table(chunk_txt, header=None)
    chunks = chunk_files[[0, 2, 3]].to_numpy()
    ref_bin_list =os.listdir(ref_bin)
    
    for i in chunks:
        
        #import binary reference
        region_id = f"{i[1]}".replace(":", "_").replace("-", "_")
        bin_target = [i for i in ref_bin_list if region_id in i][0]
        
        #phasing + imputation
        phase = f"GLIMPSE2_phase --bam-file {sample_cram} --reference {ref_bin}/{bin_target} --output {sample_id}_{chr_id}.imputed.{i[0]}.bcf --fasta {fasta}"
        index = f"bcftools index -f {sample_id}_{chr_id}.imputed.{i[0]}.bcf"
        
        #execute command
        try:
            logging.info(f"Phasing + Imputing {sample_id} on region {i[1]} ...")
            run_command(phase)
            logging.info(f"Creating Index file {sample_id} on region {i[1]} ...")
            run_command(index)
        except ValueError:
            logging.warning(f"Failed to Impute {sample_id} on region {i[1]} ...")

#Function to run the entire thing
def main():
    parameters = load_parameters()
    sample_id = parameters.sample
    chr_id = parameters.chromosome
    sample_cram = parameters.input
    ref_bin = parameters.reference
    chunk_txt = parameters.chunk
    fasta_file = parameters.fasta
    
    #Run function
    impute_phase(sample_id, chr_id, sample_cram, ref_bin, chunk_txt, fasta_file)

#run all
if __name__ == "__main__":
    main()
