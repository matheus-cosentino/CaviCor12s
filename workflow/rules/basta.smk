############################################################################
#                               Rules Basta                                #
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


rule setup_basta_taxonomy:
  message:
    """
    > Basta >> Get Basta Taxonomy Map <<
    > Output >> {output} <<
    """
  output:
    "resources/basta_db/taxonomy.sqlite"
  params:
    db_dir="resources/basta_db"
  conda:
    BASTA
  log:
    "logs/Databases_log/setup_basta_taxonomy.log"    
  shell:
    """
    basta download gb -d {params.db_dir} >> {log} 2>&1
    """

rule basta_search:
  message:
    """
    > Basta >> LCA Mito12S  <<
    > Input >> {input.blast} <<
    > Output >> {output} <<
    > Identity >> {pident} <<
    """
  input:
    blast="{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt",
    tax_db="resources/basta_db/taxonomy.sqlite"
  output:
    "{out_dir}/{sample}/Basta/{sample}_{pident}_LCA_Taxonomy.txt"
  conda:
    BASTA
  threads:
    2
  log:
    "{out_dir}/{sample}/Basta/{sample}_{pident}_LCA_Taxonomy.log"
  params:
    lca_depth = config["basta"]["lca_depth"],
    best_hit = config["basta"]["best_hit"],
    percentage = config["basta"]["percentage"] 
  shell:
    """   
    basta sequence {input.blast} {output}  {input.tax_db} --best_hit {params.best_hit} --majority {params.lca_depth} --percentage {params.percentage} >> {log} 2>&1
    """

