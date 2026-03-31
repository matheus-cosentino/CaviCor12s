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
rule multiqc_aggregate:
  message:
    """
    > Generate a multiqc report in HTML format all samples
    > Input: {input.files}
    > Output Directory: {output.report}
    """ 
  conda:
    MULTIQC
  input:
    # Use a dedicated function to get only files that MultiQC can parse.
    files = get_multiqc_inputs()
  output:
    report = os.path.join(OUT_DIR, "multiqc_all", "{pident}_multiqc_report.html"),
    data_dir = directory(os.path.join(OUT_DIR, "{pident}_multiqc_report.html"))
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
      --outdir {output.data_dir} \
      --filename multiqc_report.html \
      {params.extra} \
      {input.files} > {log} 2>&1 
    """
