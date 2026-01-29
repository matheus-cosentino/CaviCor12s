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
    fastq = 
  output:
    fastq_clean = "results/{sample}/Chopper/{sample}_filtered_fastq.gz"
  params:
    qual = config["chopper"]["min_quality"],
    min_len = config["chopper"]["min_length"],
    max_len = config["chopper"]["max_length"]
  threads: 
    4
  conda:
    CHOPPER
  log:
    
  shell:
    """
    gunzip -c {input.fastq} | \
    chopper --quality {params.qual} --minlength {params.min_len} --maxlength {params.max_len} --threads {threads} 2> {log} | gzip > {output.fastq_clean} >> {log} 2>&1
    """