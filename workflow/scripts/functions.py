############################################################################
#                                 Scripts                                  #
#                          MSc. Matheus Cosentino                          #
############################################################################
#    ______     ___   ____    ____  __    ______   ______   .______        #
#   /      |   /   \  \   \  /   / |  |  /      | /  __  \  |   _  \       #
#  |  ,----'  /  ^  \  \   \/   /  |  | |  ,----'|  |  |  | |  |_)  |      #
#  |  |      /  /_\  \  \      /   |  | |  |     |  |  |  | |      /       #
#  |  `----./  _____  \  \    /    |  | |  `----.|  `--'  | |  |\  \----.  #
#   \______/__/     \__\  \__/     |__|  \______| \______/  | _| `._____|  #
#                                                                          #
############################################################################


###########################################
# --- 1. Import libraryes to be used --- #
##########################################

import os, re, glob, time, sys, subprocess, platform, yaml
import urllib.request
from snakemake.io import expand
from collections import defaultdict
from Bio import SeqIO
import gzip

################################
# --- 2. Global Variables --- #
###############################
OUT_DIR = ""
BLAST_IDENTITIES = []
MODULES = {}
DATA = ""
SAMPLE = []


##########################################################
# --- 3. Function to obtain final outputs per module --- #
##########################################################

def get_final_outputs():
  final_outputs = []
  if MODULES["12s_diversity"]:
   # final_outputs.extend(expand("{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt", sample=SAMPLE, pident=BLAST_IDENTITIES, out_dir=OUT_DIR))    
    #final_outputs.extend(expand("{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt", out_dir=OUT_DIR, pident=BLAST_IDENTITIES, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/{sample}/Abundance/{sample}_{pident}_Abundance_table.tsv", out_dir=OUT_DIR, pident=BLAST_IDENTITIES, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/{sample}/Fasta_by_Genus_{pident}/", out_dir=OUT_DIR, sample=SAMPLE, pident=BLAST_IDENTITIES))
  if MODULES["quality_control"]:
    # Fixes the NameError and logically collects fastp reports for the QC module.
    # Assuming fastp creates .html reports in this path.
    final_outputs.extend(expand("{out_dir}/{sample}/Fastp/{sample}_filtered.html", out_dir=OUT_DIR, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/multiqc_all/multiqc_report.html", out_dir=OUT_DIR, sample=SAMPLE))
  #if MODULES["rarefaction"]:
  # final_outputs.append(expand("{out_dir}/Samtools/{sample}_Brute_Abundancy.txt", out_dir=OUT_DIR, sample=SAMPLE))

  
  return final_outputs


###################################################
# --- 4. Get Taxonomy reports of all samples --- #
##################################################
def get_all_basta_read_outputs(wildcards):
  paths = []
  for s in SAMPLE:
    meta = SAMPLE_META.get(s)
    paths.append(expand("{out_dir}/{sample}/Basta/{sample}_{pident}_LCA_Taxonomy.txt", sample=SAMPLE, pident=BLAST_IDENTITIES, out_dir=OUT_DIR))
  return paths

#############################################
# --- 5. Get inputs for MultiQC report --- #
############################################
def get_multiqc_inputs():
  """
  Collects all files for all samples to be included in the aggregate MultiQC report.
  """
  inputs = []
  if MODULES.get("12s_diversity"):
    inputs.extend(expand("{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt", out_dir=OUT_DIR, sample=SAMPLE, pident=BLAST_IDENTITIES))
  if MODULES.get("quality_control"):
    inputs.extend(expand("{out_dir}/{sample}/Fastp/{sample}_filtered.json", out_dir=OUT_DIR, sample=SAMPLE))
  return inputs