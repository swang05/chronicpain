rm(list = ls())
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
library(data.table)
library(ReactomePA)
setwd("D://yefang/meta_msea/")
hs <- org.Hs.eg.db

results = rbind(fread("Results/msea/neuro_pain.results.txt"),
                fread("Results/msea/non_neuro_pain.results.txt")) %>%
  filter(!grepl("_", MODULE) & FDR<=0.05)
unique.modules = unique(results$MODULE)

modfile = fread("inputs/coexpress.modules.modfile.tsv", 
                colClasses = rep("character",2)) %>% 
  filter(MODULE %in% unique.modules)
length(unique(modfile$MODULE)) == length(unique.modules)

if(0){
  findSymbol <- function(symbol){
    r <- GET(paste0("http://rest.genenames.org/search/prev_symbol/",
                    symbol), content_type("application/json"))
    r = XML:::xmlToList(XML::xmlParse(r))
    return(unname(unlist(r$result$doc)[3]))
  }
  
  symbols = modfile %>% pull(GENE) %>% unique()
  part1 = AnnotationDbi::select(hs, keys = symbols,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
  remained.symbols = part1[is.na(part1$ENTREZID),]$SYMBOL
  remained.symbols2 = unlist(lapply(1:length(remained.symbols), function(x){
    print(paste0(x, "/", length(remained.symbols)))
    g = findSymbol(symbol = remained.symbols[x])
    if (is.null(g)){
      return(NA)
    }else{
      return(g)
    }
  }))
  remained = data.frame(prev = remained.symbols, current = remained.symbols2)
  part2 = AnnotationDbi::select(hs, keys = remained[!is.na(remained$current),]$current,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
  part2 = merge(remained, part2, by.x = "current", by.y = "SYMBOL") %>% 
    dplyr::select(prev, ENTREZID) %>% rename(SYMBOL=prev)
  symbols = rbind(part1, part2)
  symbols %>% write.table("symbols2entrezID.tsv", sep = "\t", row.names = F, quote = F)
}

symbols = fread("symbols2entrezID.tsv")
modfile = modfile %>% merge(., symbols, by.x = "GENE", by.y = "SYMBOL")

annDF = data.frame()
for(sel.um in unique.modules){
  print(sel.um)
  sel.mod = modfile[modfile$MODULE==sel.um & !is.na(modfile$ENTREZID),]
  
  kk = enrichKEGG(gene = sel.mod$ENTREZID, 
                  organism = "hsa", pvalueCutoff = 0.05,
                  pAdjustMethod = "bonferroni")
  
  reactome = enrichPathway(gene = sel.mod$ENTREZID,
                           pvalueCutoff = 0.05, pAdjustMethod = "bonferroni",
                           readable = FALSE)
  annDF = rbind(kk %>% as.data.frame() %>% mutate(set = "kegg"),
                reactome %>% as.data.frame() %>% mutate(set = "reactome")) %>%
    mutate(MODULE = sel.um) %>%
    rbind(annDF, .)
}
annDF %>% 
  write.table("Results/msea/modules.kegg.reactome.tsv", sep = "\t", quote = F, 
              row.names = F)

topAnn = annDF %>% filter(p.adjust<=0.05) %>% 
  group_by(MODULE) %>% 
  filter(p.adjust == min(p.adjust))

neuropain = fread("Results/msea/neuro_pain.results.txt") %>% filter(!grepl("_", MODULE) & FDR<=0.05)
nonneuropain = fread("Results/msea/non_neuro_pain.results.txt") %>% filter(!grepl("_", MODULE) & FDR<=0.05)

ggvenn::ggvenn(list(neuropain = neuropain$MODULE,
                    nonneuropain = nonneuropain$MODULE))

intersected = topAnn %>% 
  dplyr::filter(MODULE %in% intersect(neuropain$MODULE, nonneuropain$MODULE))

neuropainOnly = topAnn %>% 
  dplyr::filter(MODULE %in% setdiff(neuropain$MODULE, nonneuropain$MODULE))

nonneuropainOnly = topAnn %>% 
  dplyr::filter(MODULE %in% setdiff(nonneuropain$MODULE, neuropain$MODULE))

ggvenn::ggvenn(list(neuropain = neuropainOnly$Description,
                    nonneuropain = nonneuropainOnly$Description))

intersect(neuropainOnly$Description,
          nonneuropainOnly$Description)

