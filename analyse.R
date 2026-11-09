# exploratory data analysis
#

library(dplyr)
library(data.table)
library(ggplot2)

if (exists("strep_isolates_unfilt") == FALSE) {
strep_isolates_unfilt <- as.data.frame(read.csv("data/strep_isolates.tsv", sep ="\t"))
}
organism_group <-unique( strep_isolates $X.Organism.group)

strep_species <- c("Streptococcus pneumoniae", "Streptococcus iniae"
, "Streptococcus agalactiae", "Streptococcus pyogenes" ,  "Streptococcus suis"
,"Streptococcus equi"   ,    "Streptococcus mutans" )


strep_isolates <- strep_isolates_unfilt |>
  dplyr::filter(grepl("Streptococcus", X.Organism.group)) |> #filter to only correct genus
  dplyr::mutate(date_corr  = date)

locations <- strep_isolates |>
  dplyr::pull(Location) |>
  as.factor() |>
  table()


# first let's restrict ourselves to Streptococcus mutans
strep_mutans <- strep_isolates |>
  dplyr::filter(X.Organism.group == "Streptococcus mutans")

strep_pyogenes <- strep_isolates |> # c 68K records
  dplyr::filter(X.Organism.group == "Streptococcus pyogenes")
par(las=1)

# histogram of species

strep_isolates_conc <- data.frame(lapply(strep_isolates, function(x) {
                      gsub("Streptococcus ", "", x)
                  }))

count_species <- strep_isolates_conc |>
  dplyr::pull(X.Organism.group) |>
  table()

# should build a georeferenced
#
