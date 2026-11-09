here::i_am("mobile-elements/mobile_element_analysis.R")

library(data.table)
library(dplyr)
library(forcats)
library(ggplot2)
library(AMR)
library(cowplot)
library(ggpubr)
library(here)

# read the MicroBIGG-E export from next to this script, wherever R was started
mb <- data.table::fread(input = here::here("mobile-elements", "microbigge.tsv"))

message("read in file microbigge.tsv, data have ", nrow(mb), " rows")

## filter to Streptococcus ------------------------------------------------------

# rename first, then filter on the clean names
mb_streptococcus <- mb |>
  dplyr::rename(
    scientific_name        = "Scientificname",
    element_name           = "Element name",
    element_symbol         = "Element symbol",
    pct_coverage_reference = "% Coverage of reference",
    pct_ident_to_reference = "% Identity to reference"
  ) |>
  dplyr::filter(grepl("Streptococcus", scientific_name, fixed = TRUE)) |>
  dplyr::mutate(
    Type = forcats::as_factor(Type),
    element_length = abs(Stop - Start) + 1
  )

message("Streptococcus rows: ", nrow(mb_streptococcus))

mb_streptococcus_pneumo <- mb_streptococcus |>
  dplyr::filter(grepl("Streptococcus pneumoniae", scientific_name, fixed = TRUE))

message("S. pneumoniae rows: ", nrow(mb_streptococcus_pneumo))

## counts by element type -------------------------------------------------------

types_mobile <- levels(droplevels(mb_streptococcus$Type))

# Type is categorical, so this is a bar chart of counts, not a histogram
by_type <- ggpubr::ggbarplot(
  data = dplyr::count(mb_streptococcus_pneumo, Type),
  x = "Type", y = "n",
  fill = "thistle3"
) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(y = "Count (log scale)") +
  cowplot::theme_cowplot()

print(by_type)

## AMR genes: which drugs they relate to ---------------------------------------

amr_elements <- mb_streptococcus_pneumo |>
  dplyr::filter(Type == "AMR")

#' Shorten mobile element names
#'
#' Replaces antimicrobial names with their AMR-package abbreviations,
#' "Streptococcus pneumoniae" with "S.pn." and "beta" with the Greek letter.
#'
#' @param df data frame with an `element_name` column
#' @return `df` with `element_name` shortened (character)
shorten_element_names <- function(df) {
  stopifnot("element_name" %in% names(df))

  antimicrobial_names <- unique(c(
    as.character(AMR::antimicrobials$name),
    unlist(AMR::antimicrobials$synonyms, use.names = FALSE)
  ))
  antimicrobial_names <- antimicrobial_names[
    !is.na(antimicrobial_names) & nchar(antimicrobial_names) >= 4
  ]
  # longest first, so "amoxicillin/clavulanic acid" wins over "amoxicillin"
  antimicrobial_names <- antimicrobial_names[
    order(nchar(antimicrobial_names), decreasing = TRUE)
  ]
  antimicrobial_abbreviations <- as.character(
    AMR::as.ab(antimicrobial_names, info = FALSE)
  )
  keep <- !is.na(antimicrobial_abbreviations) & nzchar(antimicrobial_abbreviations)
  antimicrobial_names <- antimicrobial_names[keep]
  antimicrobial_abbreviations <- antimicrobial_abbreviations[keep]

  shorten_label <- function(label) {
    # exact, whole-word, case-insensitive matches only: fuzzy matching
    # rewrote unrelated words that happened to be close to a drug name
    for (index in seq_along(antimicrobial_names)) {
      pattern <- paste0("\\b", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1",
                                    antimicrobial_names[[index]]), "\\b")
      label <- gsub(pattern, antimicrobial_abbreviations[[index]], label,
                    ignore.case = TRUE, perl = TRUE)
    }
    label
  }

  original_labels <- as.character(df$element_name)
  # species name first, so a drug pattern can never land inside it
  original_labels <- gsub("Streptococcus pneumoniae", "S.pn.",
                          original_labels, fixed = TRUE)

  unique_labels <- unique(original_labels)
  shortened <- vapply(unique_labels, shorten_label, character(1),
                      USE.NAMES = FALSE)
  labels <- shortened[base::match(original_labels, unique_labels)]
  labels <- gsub("\\bbeta\\b", "β", labels, ignore.case = TRUE, perl = TRUE)

  df$element_name <- labels
  df
}

amr_elements <- shorten_element_names(amr_elements) |>
  # most frequent at the top of the horizontal bar chart
  dplyr::mutate(element_name = forcats::fct_rev(forcats::fct_infreq(element_name)))

by_name <- ggplot2::ggplot(
  data = amr_elements,
  ggplot2::aes(y = element_name)
) +
  ggplot2::geom_bar(fill = "thistle3") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(x = "Count (log scale)", y = NULL) +
  ggplot2::theme_bw()

print(by_name)

ggplot2::ggsave(here::here("mobile-elements", "by_type.pdf"), by_type,
                width = 6, height = 4)
ggplot2::ggsave(here::here("mobile-elements", "amr_by_name.pdf"), by_name,
                width = 8, height = max(4, 0.2 * nlevels(amr_elements$element_name)))