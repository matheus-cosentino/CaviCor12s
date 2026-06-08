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
# --- 1. Import libraries to be used --- #
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
PHYLO_TARGET = []
MODULES = {}
DATA = ""
SAMPLE = []
GENE = []

##########################################################
# --- 3. Function to get final outputs per module --- #
##########################################################

def get_final_outputs():
  final_outputs = []
  if MODULES.get("12s_diversity"):
    rank_list = getattr(sys.modules[__name__], 'LCA_RANKS', None)
    final_outputs.extend(expand("{out_dir}/{sample}/Abundance/{sample}_{pident}_Abundance_{rank}.tsv", out_dir=OUT_DIR, pident=BLAST_IDENTITIES, sample=SAMPLE, rank=rank_list))
    final_outputs.extend(expand("{out_dir}/{sample}/Fasta_by_LCA_{pident}_{rank}/", out_dir=OUT_DIR, sample=SAMPLE, pident=BLAST_IDENTITIES, rank=rank_list))
  if MODULES.get("quality_control"):
    final_outputs.extend(expand("{out_dir}/{sample}/Fastp/{sample}_filtered.html", out_dir=OUT_DIR, sample=SAMPLE))    
    final_outputs.extend(expand("{out_dir}/multiqc_all/{pident}_multiqc_report.html", out_dir=OUT_DIR, pident=BLAST_IDENTITIES))
  if MODULES.get("phylogeny"):
    final_outputs.extend(expand("{out_dir}/{sample}/Phylo/{sample}_{taxid}_{pident}_{gene}_Aligned.fasta", out_dir=OUT_DIR, sample=SAMPLE, taxid=PHYLO_TARGET, pident=BLAST_IDENTITIES, gene=GENE))
    final_outputs.extend(expand("{out_dir}/{sample}/Phylo/{sample}_{taxid}_{pident}_{gene}_aligned_tree.nwk", out_dir=OUT_DIR, sample=SAMPLE, taxid=PHYLO_TARGET, pident=BLAST_IDENTITIES, gene=GENE))
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
def get_multiqc_inputs(wildcards):
  inputs = []
    # 3. Quality Data (Fastp)
  if MODULES.get("quality_control"):
    inputs.extend(expand("{out_dir}/{sample}/Fastp/{sample}_filtered.json", 
                         out_dir=OUT_DIR, 
                         sample=SAMPLE))
    inputs.extend(expand("{out_dir}/{sample}/Vsearch/{sample}_vsearch_mqc.tsv", 
                       out_dir=OUT_DIR, 
                       sample=SAMPLE))

  # 2. Diversity Data (LCA/Blast)
  if MODULES.get("12s_diversity"):
    rank_list = getattr(sys.modules[__name__], 'LCA_RANKS', None)
    if not rank_list:
      rank_list = ["genus"]
    if isinstance(rank_list, str):
      rank_list = [rank_list]

    inputs.extend(expand("{out_dir}/{sample}/LCA/{sample}_{pident}_{rank}_mqc.tsv", 
                         out_dir=OUT_DIR, 
                         sample=SAMPLE, 
                         pident=wildcards.pident,
                         rank=rank_list))
    # BLAST Metrics
    inputs.extend(expand("{out_dir}/{sample}/Blast/{sample}_{pident}_blast_summary_mqc.tsv", 
                         out_dir=OUT_DIR, 
                         sample=SAMPLE, 
                         pident=wildcards.pident))


  return inputs

####################################
# --- 6. Get NCBI API Key Info --- #
####################################
def get_ncbi_api_key(key_file_path):
  """
  Read the NCBI API key from a file and return a single-line key.

  Behaviour:
  - If a line containing 'api_key=' exists, returns the value after '='.
  - Otherwise returns the first non-empty, non-comment line.
  - Returns an empty string if no plausible key is found.
  """
  if not key_file_path or not os.path.exists(key_file_path):
    return ""

  with open(key_file_path, "r") as f:
    for raw in f:
      line = raw.strip()
      if not line:
        continue
      # ignore comment lines
      if line.startswith('#'):
        continue
      # handle key=value format
      if 'api_key=' in line:
        # split on first '=' and strip
        parts = line.split('=', 1)
        key = parts[1].strip() if len(parts) > 1 else ''
        if key:
          return key
        else:
          continue
      # if line looks like a raw key (no spaces, not too long), accept it
      if ' ' not in line and 5 < len(line) < 200:
        return line

  return ""

###########################################################
# --- 7. Extract Gene Sequences from BLAST Database --- #
###########################################################
def extract_genes(gene_query, blast_db, taxid, output_fasta, blast_results_file, threads=4, evalue=1e-10):
  """
  Run BLAST and extract gene sequences from mitochondrial database.
  
  Args:
    gene_query: Path to query gene FASTA file
    blast_db: Path to BLAST database
    taxid: Taxonomy ID to filter
    output_fasta: Output FASTA file with extracted gene sequences
    blast_results_file: TSV file with BLAST results (intermediate)
    threads: Number of threads for BLAST
    evalue: E-value threshold for BLAST
  """
  
  # Step 1: Run BLASTN to find gene sequences
  print(f"Running BLAST for {gene_query} against {blast_db}...", file=sys.stderr)
  
  blast_cmd = [
    "blastn",
    "-query", gene_query,
    "-db", blast_db,
    "-taxids", str(taxid),
    "-outfmt", "6 sseqid sstart send sstrand",
    "-evalue", str(evalue),
    "-max_target_seqs", "1",
    "-num_threads", str(threads),
    "-out", blast_results_file
  ]
  
  result = subprocess.run(blast_cmd, capture_output=True, text=True)
  if result.returncode != 0:
    print(f"BLAST error: {result.stderr}", file=sys.stderr)
    sys.exit(1)
  
  print(f"BLAST results saved to {blast_results_file}", file=sys.stderr)
  
  # Step 2: Extract sequences using blastdbcmd
  print(f"Extracting sequences from BLAST hits...", file=sys.stderr)
  
  # Check if blast_results file has content
  if not os.path.exists(blast_results_file) or os.path.getsize(blast_results_file) == 0:
    print(f"Warning: No BLAST results found. Creating empty output file.", file=sys.stderr)
    open(output_fasta, 'w').close()
    return
  
  with open(blast_results_file) as f, open(output_fasta, 'w') as out_f:
    for line in f:
      parts = line.strip().split('\t')
      if len(parts) < 4:
        continue
      
      seqid, start, end, strand = parts[0], parts[1], parts[2], parts[3]
      start, end = int(start), int(end)
      
      # Build blastdbcmd command
      strand_flag = "-" if strand == "minus" else "+"
      
      blastdbcmd_cmd = [
        "blastdbcmd",
        "-db", blast_db,
        "-entry", seqid,
        "-range", f"{start}-{end}",
        "-strand", strand_flag
      ]
      
      result = subprocess.run(blastdbcmd_cmd, capture_output=True, text=True)
      if result.returncode == 0:
        out_f.write(result.stdout)
      else:
        print(f"Warning: Failed to extract {seqid} {start}-{end}: {result.stderr}", 
              file=sys.stderr)
  
  print(f"Gene sequences saved to {output_fasta}", file=sys.stderr)

