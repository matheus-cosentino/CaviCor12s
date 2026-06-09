######################################################################
#                            Rule DB Phylo                           #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################


rule get_fasta_from_db:
  message:
    """
    > BLAST >> Extract Fasta Sequences from Mitochondrial BLAST Database <<
    > Output >> {output.fasta} <<
    """
  input:
    db = "resources/blast_db/mito"
  output:
    fasta = "resources/blast_db/FastaTaxid/{phylo_target}.fasta"
  threads:
    1
  conda:
    BLAST
  log:
    "logs/Databases_log/get_{phylo_target}_fasta.log"
  shell:
    """
    export BLASTDB="resources/blast_db/mito:${{BLASTDB:-}}"
    blastdbcmd -db resources/blast_db/mito/mito -taxids {wildcards.phylo_target} -out {output.fasta}
    """

rule blast_gene_hits:
  message:
    """
    > BLAST >> Find gene hits for {wildcards.gene} in mitochondrial database (taxid: {wildcards.taxid}) <<
    > Output >> {output.blast_results} <<
    """
  input:
    gene_query = "resources/genes/{gene}.fasta",
    db = "resources/blast_db/mito"
  output:
    blast_results = "resources/blast_db/{gene}/{gene}_{taxid}_blast_results.tsv"
  threads:
    4
  conda:
    BLAST
  log:
    "logs/Databases_log/blast_{gene}_{taxid}.log"
  shell:
    """
    export BLASTDB="resources/blast_db/mito:${{BLASTDB:-}}"
    blastn -query {input.gene_query} \
           -db resources/blast_db/mito/mito \
           -taxids {wildcards.taxid} \
           -outfmt "6 sseqid sstart send sstrand" \
           -evalue 1e-10 \
           -max_target_seqs 200 \
           -num_threads {threads} \
           -out {output.blast_results} 2>> {log}
    """

rule extract_gene_raw_fasta:
  message:
    """
    > BLASTDBCMD >> Extract raw gene FASTA for {wildcards.gene} (taxid: {wildcards.taxid}) <<
    > Output >> {output.raw_fasta} <<
    """
  input:
    blast_results = "resources/blast_db/{gene}/{gene}_{taxid}_blast_results.tsv",
    db = "resources/blast_db/mito"
  output:
    raw_fasta = "resources/blast_db/{gene}/{taxid}.raw"
  threads:
    1
  conda:
    BLAST
  log:
    "logs/Databases_log/extract_{gene}_{taxid}.log"
  shell:
    """
    export BLASTDB="resources/blast_db/mito:${{BLASTDB:-}}"
    mkdir -p $(dirname {output.raw_fasta})
    > {output.raw_fasta}
    if [ -s {input.blast_results} ]; then
      while IFS=$'\t' read -r seqid start end strand; do
        if [[ -z "$seqid" ]]; then continue; fi
        strand_flag="plus"
        [[ "$strand" == "minus" ]] && strand_flag="minus"
        blastdbcmd -db resources/blast_db/mito/mito \
                   -entry "$seqid" \
                   -range "$start-$end" \
                   -strand "$strand_flag" >> {output.raw_fasta} 2>> {log}
      done < {input.blast_results}
    fi
    """

rule normalize_gene_fasta:
  message:
    """
    > Format FASTA headers for {wildcards.gene} gene sequences (taxid: {wildcards.taxid}) <<
    > Output >> {output.gene_fasta} <<
    """
  input:
    raw_fasta = "resources/blast_db/{gene}/{taxid}.raw"
  output:
    gene_fasta = "resources/blast_db/{gene}/{taxid}.fasta"
  threads:
    1
  shell:
    """
    mkdir -p $(dirname {output.gene_fasta})
    awk 'BEGIN {{ OFS = "" }}
         /^>/ {{
           header = substr($0, 2)
           split(header, parts, " ")
           id = parts[1]
           gsub(/[:;=,\\[\\]()|-]/, "_", id)
           species = ""
           if (length(parts) >= 3) {{
             species = parts[2] "_" parts[3]
           }} else if (length(parts) == 2) {{
             species = parts[2]
           }}
           if (species == "") {{
             print ">" id
           }} else {{
             print ">" id "_" species
           }}
           next
         }}
         {{ print }}
    ' {input.raw_fasta} > {output.gene_fasta}
    """