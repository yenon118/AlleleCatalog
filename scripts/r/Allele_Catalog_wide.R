rm(list=ls())

library(dplyr)


input_file <- file.path("/storage/htc/joshilab/yenc/projects/Allele_Catalog/output/Nebraska_Chr04.txt")

output_file <- file.path("/storage/htc/joshilab/yenc/projects/Allele_Catalog/completed/Nebraska_Chr04.txt")


if(!file.exists(input_file)){
  print("Input file does not exists!!!")
  quit(1)
}

if(!dir.exists(dirname(output_file))){
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  if(!dir.exists(dirname(output_file))){
    print("Output directory cannot be created!!!")
    quit(1)
  }
}


dat <- read.table(
  file = input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

df <- dat %>%
  arrange(Position) %>%
  group_by(Classification, Improvement_Status, Maturity_Group, Country, State, Accession, Chromosome, Gene) %>%
  mutate(
    Position = paste0(Position, collapse = " "),
    Genotype = paste0(Genotype, collapse = " "),
    Genotype_with_Description = paste0(Genotype_with_Description, collapse = " ")
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)


write.table(
  x = df,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)
