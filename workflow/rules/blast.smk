###########################################################################
#                              Rules Blast                                #
#                         MSc. Matheus Cosentino                          #
###########################################################################
#   ______     ___   ____    ____  __    ______   ______   .______        #
#  /      |   /   \  \   \  /   / |  |  /      | /  __  \  |   _  \       #
# |  ,----'  /  ^  \  \   \/   /  |  | |  ,----'|  |  |  | |  |_)  |      #
# |  |      /  /_\  \  \      /   |  | |  |     |  |  |  | |      /       #
# |  `----./  _____  \  \    /    |  | |  `----.|  `--'  | |  |\  \----.  #
#  \______/__/     \__\  \__/     |__|  \______| \______/  | _| `._____|  #
#                                                                         #
###########################################################################


### get blast database
rule get_mito_db:
  output:
      directory("resources/blast_db/mito")
  log:
      "results/log/get_mito_db.log"
  shell:
      """
      mkdir -p {output}
      cd {output}
      # Baixa o banco pré-formatado (muito mais rápido que fasta bruto)
      wget https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz   > {log} 2>&1
      tar -xzvf mito.tar.gz
      rm mito.tar.gz 
      """

### blast mito
rule blast_mito:
  input:
    query="results/{sample}_consensus_cluster.fasta", 
    db_dir="resources/blast_db/mito"
  output:
    "results/{sample}/{sample}_{pident}_blast_mito.txt"
  threads: 
    8
  params:
    max_target_seqs=config["blast"]["max_target_seqs"],
    evalue=config["blast"]["evalue"],
    perc_identity=config["blast"]["perc_identity"]
    
  shell:
    """
    DB_PATH="{input.db_dir}/mito"    
    blastn -query {input.query} -db $DB_PATH -num_threads {threads} -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} -perc_identity {params.perc_identity} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" > {output}
    """