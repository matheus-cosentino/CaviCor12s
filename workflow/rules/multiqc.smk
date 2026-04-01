############################################################################
#                               Rules MultiQC                              #
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

rule blast_summary_mqc:
  input:
    blast_out = "{out_dir}/{sample}/Blast/{sample}_{pident}_Blastn_12s.txt"
  output:
    mqc_file = "{out_dir}/{sample}/Blast/{sample}_{pident}_blast_summary_mqc.tsv"
  run:
    import pandas as pd
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",  "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    try:
      df = pd.read_csv(input.blast_out, sep='\t', names=cols, header=None)
      total_hits = len(df)
      avg_ident = df['pident'].mean() if total_hits > 0 else 0
      max_bit = df['bitscore'].max() if total_hits > 0 else 0
    except:
      total_hits, avg_ident, max_bit = 0, 0, 0
    with open(output.mqc_file, 'w') as f:
      f.write("# id: blast_metrics\n")
      f.write("# section_name: 'BLAST Top Hits Metrics'\n")
      f.write("# plot_type: 'table'\n")
      f.write("Sample\tTotal Hits\tAvg Identity\tMax Bitscore\n")
      f.write(f"{wildcards.sample}\t{total_hits}\t{avg_ident:.2f}\t{max_bit}\n")

rule mqc_genus_abundance:
    input:
        lca_out = "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt"
    output:
        mqc_file = "{out_dir}/{sample}/LCA/{sample}_{pident}_genus_mqc.tsv"
    run:
        from collections import defaultdict
        import re
        
        genus_counts = defaultdict(int)
        processed_ids = set()

        with open(input.lca_out, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2: continue
                
                query_id = parts[0]
                if query_id in processed_ids: continue
                processed_ids.add(query_id)

                # Extrai abundância do cluster (size=N) 
                size_match = re.search(r"size=(\d+)", query_id)
                count = int(size_match.group(1)) if size_match else 1
                
                # Extrai o Gênero (5º nível da linhagem) 
                lineage = parts[-1]
                taxa = [t for t in lineage.split(';') if t]
                genus = taxa[4] if len(taxa) >= 5 else (taxa[-1] if taxa else "Unclassified")
                
                genus_counts[genus] += count

        # Escreve no formato de Custom Content do MultiQC
        with open(output.mqc_file, 'w') as f:
            # Cabeçalhos de configuração do MultiQC
            f.write("# id: 'genus_abundance_plot'\n")
            f.write("# section_name: 'Abundância Taxonômica (Gênero)'\n")
            f.write("# plot_type: 'bargraph'\n")
            f.write("# pconfig:\n")
            f.write("#    title: 'Leituras por Gênero'\n")
            f.write("#    ylab: 'Número de Reads'\n")
            
            # Cabeçalho da tabela (Amostra + Gêneros encontrados)
            genera = sorted(genus_counts.keys())
            f.write("Sample\t" + "\t".join(genera) + "\n")
            
            # Dados da amostra
            counts_str = "\t".join(str(genus_counts[g]) for g in genera)
            f.write(f"{wildcards.sample}\t{counts_str}\n")


rule multiqc_aggregate:
  message:
    """
    > Generate a multiqc report in HTML format for all samples
    > Input: {input.files}
    """ 
  conda:
    MULTIQC
  input:
    files = get_multiqc_inputs
  output:
    report = os.path.join(OUT_DIR, "multiqc_all", "{pident}_multiqc_report.html"),
    # MultiQC appends '_data' to the filename, not the wildcard.
    data_dir = directory(os.path.join(OUT_DIR, "multiqc_all", "{pident}_multiqc_report_data"))
  params:
    extra = "--title 'CaviCor12s Aggregate Report'"
  log:
    os.path.join(OUT_DIR, "logs", "{pident}_multiqc_aggregate.log")
  shell:
    """
    multiqc \
      --quiet \
      --export \
      --force \
      --outdir $(dirname {output.report}) \
      --filename $(basename {output.report}) \
      {params.extra} \
      {input.files} > {log} 2>&1 
    """