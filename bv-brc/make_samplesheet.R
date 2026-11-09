# make the samplesheet for the genomes in required format for
# nf-core pipeline genomeQC
# Alexis Vargas Richards, Nov. 2026

library(dplyr)

acc <- read.csv(
  here::here(
    "bv-brc",
    "accession_genomes.acc"
  ),
  header = FALSE
)

acc <- acc |>
  dplyr::mutate(assembly = paste("assembly_",
    as.character(row_number()),
    sep = ""
  )) |>
  dplyr::rename(ncbi = V1)

write.table(
  x = acc, file = here::here(
    "bv-brc",
    "samplesheet.csv"
  ), row.names = FALSE,
  sep = ",", quote = F
)
