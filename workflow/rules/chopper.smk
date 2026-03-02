############################################################################
#                               Rule Chopper                               #
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
    fastq_clean = "{out_dir}/{sample}/Chopper/{sample}_filtered_fastq.gz"
  params:  
    qual = config["chopper"]["min_quality"][0],
    min_len = config["chopper"]["min_length"][0],
    max_len = config["chopper"]["max_length"][0]
  threads: 
    1
  conda:
    CHOPPER
  log:
    "{out_dir}/{sample}/Chopper/{sample}.log"
  shell:
    """
    gunzip -c {input.fastq} | \
    chopper --quality {params.qual} --minlength {params.min_len} --maxlength {params.max_len} --threads {threads} 2> {log} | gzip > {output.fastq_clean} 2> {log}
    """