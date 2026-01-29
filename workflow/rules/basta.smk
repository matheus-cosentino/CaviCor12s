###########################################################################
#                              Rules Basta                                #
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
  shell:
    """
    basta download gb -d {params.db_dir}
    """

rule basta_search:
  message:
    """
    > Blastn >> Mito Blastn <<
    > Input >> {input.blast} <<
    > Output >> {output} <<
    > Identity >> {wildcards.pident} <<
    """
  input:
    blast="results/{sample}/Blast{sample}_{pident}_Blastn_12s.txt",
    tax_db="resources/basta_db/taxonomy.sqlite"
  output:
    "{out_dir}/{sample}/Basta/{sample}_{pident}_LCA_Taxonomy.txt"
  shell:
    """   
    basta sequence {input.blast} {output}  {input.tax_db} --best_hit 1 --majority 90 --percentage 97
    """

