############################################################################
#                                Rule Blast                                #
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


rule get_mito_db:
  message:
    """
    > Blast >> Get Mito Blast DB <<
    > Output >> {output} <<
    """
  output:
    directory("resources/blast_db/mito")
  conda:
    BLAST
  log:
    "results/Databases_log/get_mito_db.log"
  shell:
    """
    mkdir -p {output}
    echo "Downloading DB..." > {log}
    wget "https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz" -O {output}/mito.tar.gz >> {log} 2>&1
    echo "Extracting Files..." >> {log}
    tar -xzvf {output}/mito.tar.gz -C {output} >> {log} 2>&1
    rm {output}/mito.tar.gz
    echo "Done!" >> {log}
    """

rule blast_mito:
  message:
    """
    > Blastn >> Mito Blastn <<
    > Input >> {input.query} <<
    > Output >> {output} <<
    > Identity >> {wildcards.pident} <<
    """
  input:
    query="results/{sample}/Vsearch/{sample}_consenso.fasta", 
    db_dir="resources/blast_db/mito"
  output:
    "results/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
  threads: 
    4
  params:
    max_target_seqs=config["blast"]["max_target_seqs"][0],
    evalue=config["blast"]["evalue"][0]
  conda:
    BLAST
  log:
    "results/{sample}/Blast/{sample}_{pident}_Blastn_12s.log"
  shell:
    """
    DB_PATH="{input.db_dir}/mito"    
    blastn -query {input.query} -db $DB_PATH -num_threads {threads} -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} -perc_identity {wildcards.pident} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" > {output} >> {log} 2>&1
    """