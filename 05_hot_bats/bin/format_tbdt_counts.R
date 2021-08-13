library(here)
library(tidyverse)

## Read timeseries dataset
ts.counts <- read_tsv(here("data", "siderophore_abundance_hot_bats.tsv"))

ts <- read_tsv(here("data", "hot_bats_metadata.tsv"))

ts2export <- left_join(ts.counts, ts) %>%
  select(location=section, date, depth, RA) %>%
  mutate(RA=case_when(RA < 0.05 ~ 0.05,
                      RA > 0.5 ~ 0.5, 
                      TRUE ~ RA)) %>%
  mutate(logRA=log10(RA))

write_tsv(ts2export, here("data", "timeseries2plot.tsv"))
