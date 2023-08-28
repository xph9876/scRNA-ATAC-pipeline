#!/usr/bin/env Rscript

# Set a CRAN mirror
options(repos = "https://cloud.r-project.org")

# Import functions
source("utils.R")

# Argument parser
parser <- OptionParser()
parser <- add_option(parser, c("-c", "--config"), default="config.yaml",
  help="Config file in YAML format (config.yaml)"
)
args <- parse_args(parser)

# Fetch configs from yaml file
configs <- yaml.load_file(args$config)

# Output configs
cat("Current configs:\n")
str(configs)

# Create output folders if needed
if (!dir.exists(configs$output)) {
  dir.create(configs$output)
}

# Preprocess data
data <- preprocess(
      name=configs$sample_name,
      config=configs
    )
cat("[INFO] Data preprocessed!\n")

# Cluster and annotation
data <- cluster(
  data,
  config=configs,
  output=configs$output
)
cat("[INFO] Clustering figures generated!\n")

