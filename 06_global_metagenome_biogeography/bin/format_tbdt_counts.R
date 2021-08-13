library(here)
library(tidyverse)

## Formatting global dataset

global <- read_tsv(here::here("data", "siderophore_abundance_global.tsv"))
meta <- read_tsv(here::here("data", "geotraces_tara_combined_imputed.tsv"))

### Simons
simons2export <- semi_join(global, meta) %>% 
  left_join(meta) %>%
  filter(!(str_detect(sampleID, "tara_"))) %>%
  mutate(dcm.dummy=as.numeric(dcm.layer)) %>%
  group_by(lat, lon) %>%
  filter(dcm.dummy==max(dcm.dummy)) %>%
  filter(RA==max(RA)) %>%
  mutate(RA=case_when(RA < 0.05 ~ 0.05,
                      RA > 0.5 ~ 0.5, 
                      TRUE ~ RA)) %>%
  mutate(logRA=log10(RA)) %>%
  select(lat, lon, section, RA, logRA)

write_tsv(simons2export, here("data", "simons2plot.tsv"))


### Tara
tara2export <- semi_join(global, meta) %>% 
  left_join(meta) %>%
  filter(str_detect(sampleID, "tara_")) %>%
  filter(dcm.layer==1) %>%
  group_by(lat, lon) %>%
  filter(RA==max(RA)) %>%
  mutate(RA=case_when(RA < 0.05 ~ 0.05,
                      RA > 0.5 ~ 0.5, 
                      TRUE ~ RA)) %>%
  mutate(logRA=log10(RA)) %>%
  select(lat, lon, section, RA, logRA)

write_tsv(tara2export, here("data", "tara2plot.tsv"))
