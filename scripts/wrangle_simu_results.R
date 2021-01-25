#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

path <- args[1]
files <- list.files(path, pattern = "group1-vs-group2.diffexp.tsv", full.names = TRUE, recursive = TRUE)
de_list <- files %>% 
  map(read_table2, skip = 1, col_names = FALSE)
names(de_list) <- dirname(files) %>% 
  str_split("/") %>% 
  map_chr(4)

de_raw <- de_list %>% 
  bind_rows(.id = "id")

colnames(de_raw) <- c("id", "transcriptid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

de_raw <- de_raw %>% 
  mutate_at("transcriptid", str_remove_all, '\"')

de <- de_raw %>% 
  separate("id", c("reps", "n_eff"), sep = "_") %>% 
  mutate_at(c("reps", "n_eff"), str_extract, "\\d+") %>% 
  mutate_at(c("reps", "n_eff"), as.numeric) %>% 
  mutate_at("n_eff", ~ 2 * n_eff)

de_simulation <- de %>% 
  mutate(
    reps = factor(str_c("N=", reps), levels = c("N=3", "N=6", "N=10")),
    n_eff = factor(str_c("Effects=", n_eff), levels = c("Effects=100", "Effects=200", "Effects=400", "Effects=800"))
  )

de_simulation %>% 
  write_csv("output/de_simulation_results.csv")
