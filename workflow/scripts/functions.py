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
    final_outputs.extend(expand("{out_dir}/{sample}/Basta/{sample}_{pident}_LCA_Taxonomy.txt", sample=SAMPLE, pident=BLAST_IDENTITIES, out_dir=OUT_DIR))
    #final_outputs.extend(expand("{out_dir}/{sample}/Krona/{sample}_Basta_Krona.html", out_dir=OUT_DIR, sample=SAMPLE))
    
  if MODULES["rarefaction"]:
    final_outputs.append(expand("{out_dir}/Rarefaction_Curve/Basta_Rarefaction_Curve.pdf", out_dir=OUT_DIR)),
    final_outputs.append(expand("{out_dir}/{sample}/Rarefaction_Curve/{sample}_Basta_Rarefaction_Curve.pdf", out_dir=OUT_DIR, sample=SAMPLE))
  
  if MODULES["abundancy"]:
    final_outputs.append(expand("{out_dir}/Samtools/{sample}_Brute_Abundancy.txt", out_dir=OUT_DIR, sample=SAMPLE))
  
  return final_outputs


###################################################
# --- 3. Get Taxonomy reports of all samples --- #
##################################################
def get_all_basta_read_outputs(wildcards):
  paths = []
  for s in SAMPLE:
    meta = SAMPLE_META.get(s)
    paths.append(expand("{out_dir}/{sample}/Basta/{sample}_{pident}_LCA_Taxonomy.txt", sample=SAMPLE, pident=BLAST_IDENTITIES, out_dir=OUT_DIR))
  return paths