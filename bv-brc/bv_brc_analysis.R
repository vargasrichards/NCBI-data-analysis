# analysis of Streptococcus pneumoniae BV-BRC genomes
# Alexis Vargas Richards
# Taxon ID 1313

here::here("bv-brc/bv_brc_analysis.R")
library(here)
library(dplyr)
library(data.table)
library(ggplot2)

library(countries)

path_to_strep <- here::here("bv-brc", "BVBRC_genome (1).csv")

bvbrc_spn <- data.table::fread(path_to_strep, sep = ",") |>
  dplyr::select(!c(Kingdom, Phylum, Class,
   Order, Family, Genus, `NCBI Taxon ID`, County,
   `H1 Clade Global`, `H1 Clade US`, `H5 Clade`,
   `pH1N1-like`, H_type, N_type, Segment, Passage, Authors)) |>
     dplyr::filter(`Genome Quality` != "Poor") |> # backticks!!
  dplyr::mutate(
    isolation_country =
      forcats::fct_infreq(
        countries::country_name(`Isolation Country`,
        verbose = TRUE
      ))
  )
head(bvbrc_spn)


# genome_qualities <- bvbrc_spn |>
#   dplyr::pull("Genome Quality") |>
#   table()

# # Poor quality genomes to be removed
# bvbrc_spn_filt <- bvbrc_spn |>
#   dplyr::filter(`Genome Quality` != "Poor") # backticks!!

# barplot of countries 

bvbrc_spn |> dplyr::pull(isolation_country) |> table()

ggplot(data = bvbrc_spn, aes(y = isolation_country)) + geom_bar() +
 theme_bw(base_size = 12)


# num high quality genomes
num_presumptive_hq_genomes <- nrow(bvbrc_spn)

# accessions for these genomes
num_assemblies_empty <- bvbrc_spn |>  
    dplyr::filter(`Assembly Accession` == "") |>
    nrow()

num_assemblies_nonempty <- num_presumptive_hq_genomes - num_assemblies_empty
pct_empty <- 100 * num_assemblies_empty / num_presumptive_hq_genomes

# which ones have assembled genomes??

hq_with_assembly <- bvbrc_spn |>
    dplyr::filter(!(`Assembly Accession` == ""))

assembly_needed <- hq_with_assembly |> dplyr::pull(`Assembly Accession`)

# can use ncbi datasets tool here.
write.table(x = assembly_needed,
 file = here::here("bv-brc", "accession_genomes.acc"),
 row.names = FALSE, quote = FALSE, col.names = FALSE)

fetch <- FALSE
if (fetch) {
system(glue("/Users/alexis/miniforge3/envs/ncbi_datasets/bin/datasets download genome accession --inputfile bv-brc/accession_genomes.acc --include gff3,genome,seq-report"))
}

