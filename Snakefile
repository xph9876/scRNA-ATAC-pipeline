# Load config file
configfile: "config.yaml"


rule all:
  input:
    config['output']


# Generate count matrix from input
rule cellranger_count:
  input:
    csv=config['input_csv'],
    ref=config['ref_path']
  output:
    directory(config['cellranger_output_path'])
  shell:
    "cellranger-arc count --id {output} --reference {input.ref} --libraries {input.csv}"

  
# Analyzing using R
rule analysis:
  input:
    config['cellranger_output_path']
  output:
    directory(config["output"])
  shell:
    "Rscript analyze.R"


  
