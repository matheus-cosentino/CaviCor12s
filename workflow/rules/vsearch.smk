############################################################################
#                               Rules Vsearch.                             #
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

rule vsearch:
  message:
    """
    > VSearch >> Clustering Amplicon Reads <<
    > Input >> {input.chopper} <<
    > Output >> {output} <<
    """
  input:
    chopper= "results/{sample}/Chopper/{sample}_filtered.fastq"
  output:
    centroid= "results/{sample}/Vsearch/{sample}_centroids.fasta",
    consenso= "results/{sample}/Vsearch/{sample}_consenso.fasta"
  params:
    identity = config["vsearch"]["identity"][0]
  conda:
    VSEARCH
  threads: 
    4
  log:
    "results/{sample}/Vsearch/{sample}_Vsearch.log"
  shell:
    """ 
    vsearch --cluster_fast {input.chopper} --id {params.identity} --threads {threads} --centroids {output.centroid} --consout {output.consenso} --sizeout  >> {log} 2>&1
    """