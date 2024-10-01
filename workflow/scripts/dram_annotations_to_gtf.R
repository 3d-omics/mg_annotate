library(tidyverse)

read_tsv(snakemake@input[[1]]) %>%
  select(
    seqname = scaffold,
    gene_id = `...1`,
    start = start_position,
    end = end_position,
    strandedness
  ) %>%
  mutate(
    source = "DRAM",
    feature = "gene",
    score = ".",
    strand = if_else(
      strandedness == 1, "+", "-"
    ),
    frame = ".",
    attribute="."
  ) %>%
  select(
    seqname,
    source,
    feature,
    start,
    end,
    score,
    strand,
    frame,
    attribute
  ) %>%
  write_tsv(snakemake@output[[1]])
