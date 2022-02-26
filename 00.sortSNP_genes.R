rm(list = ls())
library(data.table)
library(tidyverse)
setwd("D:/yefang/meta_msea/")
dir.create("inputs/marker_dis", showWarnings = FALSE)
options(stringsAsFactors = F)

rawsnpfiles = list.files("F://projects/pain/", recursive = T, full.names = T, 
                         pattern = "ld_checked")
for (rsf in rawsnpfiles){
  print(paste0("Working on ", rsf))
  outname = file.path("inputs/marker_dis", 
                      gsub(".top50perc.ld_checked.summary", 
                           "marker_dis.tsv", 
                           strsplit(rsf, "/")[[1]][7]))
  fread(rsf) %>% select(snpName, `-log10p`) %>% 
    rename(MARKER=snpName, VALUE = `-log10p`) %>% 
    write.table(., file = outname, sep = "\t", row.names = F)
}

snpfiles = list.files("inputs/marker_dis/", pattern = "marker_dis")
files = list.files("F://ResourceData/GTEx/GTEx_Analysis_v8_eQTL/", pattern = "var_gene.pairs")
regulome = fread("F://ResourceData/ENCODE_regulome/ENCFF297XMQ.snps.highConfidence.2genes", header = F)
head(regulome)

merged = data.frame()
finalunmapped = c()
for(sf in snpfiles){
  snps = fread(file.path("inputs/marker_dis/", sf))
  unmapped = snps$MARKER
  for(f in files){
    print(f)
    inf = fread(file.path("F://ResourceData/GTEx/GTEx_Analysis_v8_eQTL/", f))
    head(inf)
    tmp = merge(snps, inf, by.x = "MARKER", by.y = "var") %>% 
      dplyr::select(MARKER, gene_name) %>% dplyr::rename(GENE = gene_name)
    unmapped = setdiff(unmapped, tmp$MARKER)
    merged = rbind(merged, tmp) %>% unique()
  }
  finalunmapped = c(finalunmapped, unmapped) %>% unique()
}

finalunmapped
head(finalunmapped)
regulome %>% filter(V1 %in% finalunmapped) %>% 
  rename(MARKER = V1, GENE=V2) %>%
  rbind(merged, .) %>%
  write.table(.,
              file = file.path("inputs/", "snps_genes.tsv"), 
              row.names = F, quote = F, sep = "\t")






