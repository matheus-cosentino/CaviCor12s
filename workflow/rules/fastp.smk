############################################################################
#                               Rule Fastp                                 #
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

rule clean_reads:
  input:
    fastq = os.path.join(DATA, "{sample}.fastq.gz")
  output:
    fastq_clean = "{out_dir}/{sample}/Fastp/{sample}_filtered_fastq.gz",
    html= "{out_dir}/{sample}/Fastp/{sample}_filtered.html",
    json= "{out_dir}/{sample}/Fastp/{sample}_filtered.json"
  params:  
    qual = config["fastp"]["min_quality"][0],
    min_len = config["fastp"]["min_length"][0],
    max_len = config["fastp"]["max_length"][0]
  threads: 
    1
  conda:
    FASTP
  log:
    "{out_dir}/{sample}/Fastp/{sample}.log"
  shell:
    """
    fastp --in1 {input.fastq} --out1 {output.fastq_clean} --html {output.html} --json {output.json} --thread {threads} --qualified_quality_phred {params.qual} --max_len1 {params.max_len} --length_required {params.min_len}  --dedup --trim_poly_g --trim_poly_x  2> {log}        
    """

