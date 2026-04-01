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
    > Blast >> Get Mito Blast DataBase <<
    > Output >> {output} <<
    """
  output:
    directory("resources/blast_db/mito")
  conda:
    BLAST
  log:
    "logs/Databases_log/get_mito_db.log"
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
    query="{out_dir}/{sample}/Vsearch/{sample}_consenso.fasta", 
    db_dir="resources/blast_db/mito"
  output:
    "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
  threads: 
    4
  params:
    max_target_seqs=config["blast"]["max_target_seqs"][0],
    evalue=config["blast"]["evalue"][0]
  conda:
    BLAST
  log:
    "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.log"
  shell:
    """
    export BLASTDB=$(pwd)/{input.db_dir}
    blastn -query {input.query} -db mito -num_threads {threads} -evalue {params.evalue} -max_target_seqs {params.max_target_seqs} -perc_identity {wildcards.pident} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames slineage" -out {output} 2>&1>> {log} 
    """

rule get_taxdump:
  message:
    """
    > TamDump DB >> Get Database <<
    > Output >> {output} <<
    """
  conda:
    TAXONKIT
  threads: 
    1
  output:
    directory("resources/taxonomy/taxdump")
  log:
    "logs/Databases_log/get_taxdump.log"
  shell:
    """
    mkdir -p {output}
    wget -c https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O {output}/taxdump.tar.gz >> {log} 2>&1
    tar -zxvf {output}/taxdump.tar.gz -C {output} >> {log} 2>&1
    """

rule get_lineages:
  message:
    """
    > LCA TamDump >> LCA Algo <<
    > Input >> {input.blast} & {input.db} <<
    > Output >> {output} <<
    """
  input:
    blast = "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt",
    db = rules.get_taxdump.output
  output:
    "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
  conda:
    TAXONKIT
  threads: 
    1
  shell:
    """
    taxonkit lca -i 13 --data-dir {input.db} {input.blast} | \
    taxonkit reformat -i 14 --data-dir {input.db} -f "{{p}};{{c}};{{o}};{{f}};{{g}}" > {output}
    
    """

rule summarize_blast_lca:
  message:
    """
    > LCA Summarize >> LCA Summary <<
    > Input >> {input.blast_out} <<
    > Output >> {output.abundance_table} <<
    """
  input:
    blast_out= "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
  output:
    abundance_table="{out_dir}/{sample}/Abundance/{sample}_{pident}_Abundance_table.tsv"
  run:
    from collections import defaultdict
    import re
    genus_abundance = defaultdict(int)
    centroids_processados = set()

    with open(input.blast_out, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            
            centroid_id = parts[0]
            if centroid_id in centroids_processados:
                continue
            
            centroids_processados.add(centroid_id)

            # Extrair o 'size' (abundância de reads no cluster)
            size_match = re.search(r"size=(\d+)", centroid_id)
            count = int(size_match.group(1)) if size_match else 1
            
            # Pegar a linhagem (LCA já calculado pelo TaxonKit)
            lineage = parts[-1] 
            taxa = lineage.split(';')
            
            # Extrair Gênero (ajuste o índice se necessário, aqui usamos o 5º nível)
            if len(taxa) >= 5:
                genus = taxa[4]
            elif len(taxa) > 0 and taxa[0] != "":
                genus = taxa[-1]
            else:
                genus = "Unclassified"

            genus_abundance[genus] += count

    # Salvar a tabela final
    with open(output.abundance_table, 'w') as out:
        out.write("Genus\tAbundance\n")
        for genus, total in sorted(genus_abundance.items(), key=lambda x: x[1], reverse=True):
            out.write(f"{genus}\t{total}\n")

rule extract_fastas_by_lca:
    message:
        """
        > Extracting Fasta by LCA taxon
        > Input: {input.lca}
        > Output Directory: {output.fasta_dir}
        """
    input:
        lca = "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt",
        query = "{out_dir}/{sample}/Vsearch/{sample}_consenso.fasta"
    output:
        fasta_dir = directory("{out_dir}/{sample}/Fasta_by_Genus_{pident}/")
    run:
        import os
        from Bio import SeqIO

        # 1. Mapear cada ID de sequência para o táxon atribuído pelo LCA
        seq_to_taxon = {}
        with open(input.lca, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2: continue
                
                # O ID da query é o primeiro campo [cite: 19, 20]
                query_id = parts[0]
                # A linhagem formatada é o último campo (após taxonkit reformat) [cite: 21]
                lineage = parts[-1]
                
                # Extração do táxon seguindo sua lógica: 
                # Tenta o gênero (5º nível), se não houver, pega o último nível disponível [cite: 22]
                taxa = [t for t in lineage.split(';') if t]
                if len(taxa) >= 5:
                    taxon = taxa[4]
                elif taxa:
                    taxon = taxa[-1]
                else:
                    taxon = "Unclassified"
                
                # Sanitização simples para evitar problemas com nomes de arquivos
                taxon_clean = taxon.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
                seq_to_taxon[query_id] = taxon_clean

        # 2. Garantir que o diretório de saída existe
        os.makedirs(output.fasta_dir, exist_ok=True)

        # 3. Ler o FASTA de entrada e distribuir as sequências
        # Usamos um dicionário para manter os arquivos abertos apenas durante a escrita
        handles = {}
        try:
            for record in SeqIO.parse(input.query, "fasta"):
                # O ID no arquivo LCA deve bater com o ID no FASTA [cite: 19, 25]
                if record.id in seq_to_taxon:
                    t_name = seq_to_taxon[record.id]
                    file_path = os.path.join(output.fasta_dir, f"{t_name}.fasta")
                    
                    if t_name not in handles:
                        handles[t_name] = open(file_path, "w")
                    
                    SeqIO.write(record, handles[t_name], "fasta")
        finally:
            for h in handles.values():
                h.close()