######################################################################
#                            Rules Phylogeny                        #
#                          MSc. Matheus Cosentino                    #
######################################################################
#   \  |  _)   |              ___|                       |           #
#  |\/ |   |   __|    _ \    |        _ \    __ \     _` |    _` |   #
#  |   |   |   |     (   |   |       (   |   |   |   (   |   (   |   #
# _|  _|  _|  \__|  \___/   \____|  \___/   _|  _|  \__,_|  \__,_|   #
#                                                                    #
######################################################################

# Wildcard constraints for phylogenetic analysis
wildcard_constraints:
    phylo_target = r"\d+"  # phylo_target should be numeric (TaxID)


rule extract_taxid_by_lca:
    message:
        """
        > Extract sequences from {wildcards.sample} matching phylo_target {wildcards.phylo_target}
        > Input >> {input.lca}, {input.query}
        > Output >> {output.taxid_fasta}
        """
    input:
        lca = "{out_dir}/{sample}/LCA/{sample}_{pident}_LCA_Lineage.txt",
        query = "{out_dir}/{sample}/Vsearch/{sample}_consenso.fasta",
        ref_db = "resources/blast_db/FastaTaxid/{phylo_target}.fasta"
    output:
        taxid_fasta = "{out_dir}/{sample}/Phylo/{sample}_{pident}_{phylo_target}_raw.fasta"
    log:
        "{out_dir}/logs/{sample}_{pident}_{phylo_target}_extract.log"
    run:
        import os
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.taxid_fasta), exist_ok=True)
        
        # Get the taxon name from config using the taxid (phylo_target)
        try:
            taxid = int(wildcards.phylo_target)
        except ValueError:
            taxid = wildcards.phylo_target
            
        taxon_name = config.get("phylo_names", {}).get(taxid)
        if not taxon_name:
            taxon_name = config.get("phylo_names", {}).get(str(taxid))
            
        if not taxon_name:
            with open(str(log), 'w') as f:
                f.write(f"Error: TaxID {taxid} not found as a key in config['phylo_names'].\n")
            raise ValueError(f"TaxID {taxid} must be mapped to a Taxon Name in config.yaml under 'phylo_names'")

        # Map sequence IDs to their lineages from LCA file
        seq_to_lineage = {}
        with open(input.lca, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                query_id = parts[0]
                lineage = parts[-1]
                seq_to_lineage[query_id] = lineage

        # Extract sequences that contain the target taxon in their lineage
        extracted_count = 0
        with open(output.taxid_fasta, 'w') as out_fh:
            for record in SeqIO.parse(input.query, "fasta"):
                if record.id in seq_to_lineage:
                    lineage = seq_to_lineage[record.id]
                    # Check if the desired taxon name is anywhere in the lineage string
                    if taxon_name in lineage.split(';'):
                        SeqIO.write(record, out_fh, "fasta")
                        extracted_count += 1
                        
        with open(str(log), 'w') as f:
            f.write(f"Extracted {extracted_count} sequences for phylo_target {taxon_name} (TaxID {taxid})\n")


rule rename_fasta_headers:
    message:
        """
        > Rename FASTA headers with sample name
        > Input >> {input.fasta}
        > Output >> {output.renamed_fasta}
        """
    input:
        fasta = "{out_dir}/{sample}/Phylo/{sample}_{pident}_{phylo_target}_raw.fasta"
    output:
        renamed_fasta = "{out_dir}/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta"
    run:
        import os
        import re
        from Bio import SeqIO

        os.makedirs(os.path.dirname(output.renamed_fasta), exist_ok=True)
        with open(output.renamed_fasta, 'w') as out_fh:
            for record in SeqIO.parse(input.fasta, "fasta"):
                base_id = record.id.split()[0]
                
                # Replace special characters that break Newick parsing (like semicolons from VSEARCH)
                clean_id = re.sub(r'[;=:,()[\]]', '_', base_id)
                
                original_desc = record.description if record.description else base_id
                record.id = f"{wildcards.sample}_{clean_id}"
                record.description = f"{wildcards.sample}|{original_desc}"
                SeqIO.write(record, out_fh, "fasta")


rule aln_ref_taxid:
    message:
        """
        > Align reference database for phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> resources/blast_db/{wildcards.gene}/{wildcards.phylo_target}.fasta
        > Output >> {output.ref_aln}
        """
    input:
        ref_fasta = "resources/blast_db/{gene}/{phylo_target}.fasta"
    output:
        ref_aln = "{out_dir}/Phylo/{phylo_target}_{gene}_ref_aln.fasta"
    conda:
        MAFFT
    threads:
        4
    log:
        "{out_dir}/logs/mafft_ref_{phylo_target}_{gene}.log"
    shell:
        """
        mafft --thread {threads} --auto {input.ref_fasta} > {output.ref_aln} 2> {log}
        """





rule addfragments_align:
    message:
        """
        > Add sample fragments to reference alignment for sample {wildcards.sample}, phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> {input.fragments} & {input.ref_aln}
        > Output >> {output.aligned}
        """
    input:
        fragments = "{out_dir}/{sample}/Phylo/{sample}_{pident}_{phylo_target}_renamed.fasta",
        ref_aln = "{out_dir}/Phylo/{phylo_target}_{gene}_ref_aln.fasta"
    output:
        aligned = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned.fasta"
    conda:
        MAFFT
    threads:
        4
    log:
        "{out_dir}/logs/mafft_addfragments_{sample}_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        if [ -s {input.fragments} ]; then
            mafft --addfragments {input.fragments} --reorder --thread {threads} {input.ref_aln} > {output.aligned} 2> {log}
        else
            cp {input.ref_aln} {output.aligned}
            echo "Input fragments file is empty. Copied reference alignment." > {log}
        fi
        """

rule clean_alignment:
    message:
        """
        > Clean alignment for sample {wildcards.sample}, phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> {input.aligned}
        > Output >> {output.cleaned_alignment}
        """
    input:
        aligned = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned.fasta"
    output:
        cleaned_alignment = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    conda:
        TRIMAL
    threads:
        4
    log:
        "{out_dir}/logs/trimal_{sample}_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        trimal -in {input.aligned} -out {output.cleaned_alignment} -automated1 > {log} 2>&1
        """


rule build_tree_with_fasttree:
    message:
        """
        > Build phylogenetic tree from alignment for sample {wildcards.sample}, phylo_target {wildcards.phylo_target} and gene {wildcards.gene}
        > Input >> {input.aligned}
        > Output >> {output.tree}
        """
    input:
        aligned = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_Aligned_cleaned.fasta"
    output:
        tree = "{out_dir}/{sample}/Phylo/{sample}_{phylo_target}_{pident}_{gene}_aligned_tree.nwk"
    conda:
        FASTTREE
    threads:
        4
    log:
        "{out_dir}/logs/fasttree_{sample}_{phylo_target}_{pident}_{gene}.log"
    shell:
        """
        fasttree -nt -gtr -gamma < {input.aligned} > {output.tree} 2> {log}
        """

# Note: Generated tree files (.nwk) can be visualized using external tools like:
# - iTOL (Interactive Tree of Life): https://itol.embl.de/
# - FigTree: http://tree.bio.ed.ac.uk/software/figtree/
# - R packages: ape, ggtree, phytools
# - Python packages: dendropy, toytree
