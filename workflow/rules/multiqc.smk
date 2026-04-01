############################################################################
#                               MultiQC Rules                              #
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
      avg_ident = df['pident'].max() if total_hits > 0 else 0
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

                # Extracts cluster abundance (size=N) 
                size_match = re.search(r"size=(\d+)", query_id)
                count = int(size_match.group(1)) if size_match else 1
                
                # Extracts Genus (5th level of the lineage) 
                lineage = parts[-1]
                taxa = [t for t in lineage.split(';') if t]
                genus = taxa[4] if len(taxa) >= 5 else (taxa[-1] if taxa else "Unclassified")
                
                genus_counts[genus] += count

        # Write in MultiQC Custom Content format
        with open(output.mqc_file, 'w') as f:
            # MultiQC configuration headers
            f.write("# id: genus_abundance_plot\n")
            f.write("# section_name: 'Taxonomic Abundance (Genus)'\n")
            f.write("# plot_type: 'bargraph'\n")
            f.write("# pconfig:\n")
            f.write("#    title: 'Reads per Genus'\n")
            f.write("#    ylab: 'Number of Reads'\n")
            
            # Table header (Sample + Found Genera)
            genera = sorted(genus_counts.keys())
            f.write("Sample\t" + "\t".join(genera) + "\n")
            
            # Sample data
            counts_str = "\t".join(str(genus_counts[g]) for g in genera)
            f.write(f"{wildcards.sample}\t{counts_str}\n")

rule vsearch_summary_mqc:
    input:
        log = os.path.join(OUT_DIR, "{sample}/Vsearch/{sample}_Vsearch.log")
    output:
        mqc_file = os.path.join(OUT_DIR, "{sample}/Vsearch/{sample}_vsearch_mqc.tsv")
    run:
        import re
        import os

        # Reading the log content
        with open(input.log, 'r') as f:
            content = f.read()

        # Regex to capture data from the log
        # E.g.: "50587495 nt in 89317 seqs"
        reads_match = re.search(r"in (\d+) seqs", content)
        # E.g.: "Clusters: 12945"
        clusters_match = re.search(r"Clusters: (\d+)", content)
        # E.g.: "Singletons: 11506, 12.9% of seqs"
        singletons_match = re.search(r"Singletons: (\d+), ([\d.]+)% of seqs", content)

        reads = reads_match.group(1) if reads_match else "0"
        clusters = clusters_match.group(1) if clusters_match else "0"
        singletons = singletons_match.group(1) if singletons_match else "0"
        perc_sing = singletons_match.group(2) if singletons_match else "0"

        # Writing the file formatted for MultiQC
        with open(output.mqc_file, 'w') as f:
            f.write("# id: vsearch_stats\n")
            f.write("# section_name: 'VSEARCH Clustering Summary'\n")
            f.write("# plot_type: 'table'\n")
            f.write("Sample\tTotal Reads\tClusters\tSingletons\t% Singletons\n")
            f.write(f"{wildcards.sample}\t{reads}\t{clusters}\t{singletons}\t{perc_sing}\n")



rule multiqc_aggregate:
  message:
    """
    > Generate a multiqc report in HTML format for all samples
    > Input: {input.files}
    """ 
  conda:
    MULTIQC
  input:
    files = get_multiqc_inputs,
    config = "config/multiqc_config.yaml"
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
      --force  \
      --outdir $(dirname {output.report}) \
      --filename $(basename {output.report}) \
      --config {input.config} \
      {params.extra} \
      {input.files} > {log} 2>&1 
    """

