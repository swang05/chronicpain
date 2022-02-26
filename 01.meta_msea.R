rm(list = ls())
library(tidyverse)
library(Mergeomics)

setwd("D://yefang/meta_msea/")

if(0){
  test = read.table("inputs/coexpress.modules.tsv", header = T)
  test = test %>% 
    separate(., col ="MODULE", sep = ":", 
             into = c("name1", "name2"), remove = F)
  test.tissues = data.frame(tissue = unique(test$name1)) %>% 
    mutate(ix = formatC(1:nrow(.), width = 2, flag = 0))
  test = test %>% mutate(name2 = formatC(as.numeric(name2), width = 3, flag = "0")) %>% 
    merge(test.tissues, by.x ="name1", by.y = "tissue") %>% 
    mutate(MODULE = paste0(as.character(ix), as.character(name2)),
           DESCR = paste0(MODULE, ":", name1), SOURCE = "human")
  test %>% dplyr::select(MODULE, GENE) %>% 
    write.table("inputs/coexpress.modules.modfile.tsv", sep = "\t", row.names = F, quote = F)
  test %>% dplyr::select(MODULE, SOURCE, DESCR) %>% unique() %>% 
    write.table("inputs/coexpress.modules.info.tsv", sep = "\t", row.names = F, quote = F)
  rm(list = c("test", "test.tissues"))
}

runMSEA <- function(label, folder = "Results", marfile, nperm = 2000){
  job.msea <- list()
  genfile = "inputs/snps_genes.tsv"
  modfile = "inputs/coexpress.modules.modfile.tsv"
  inffile = "inputs/coexpress.modules.info.tsv"
  job.msea$label <- label
  job.msea$folder <- folder
  job.msea$genfile <- genfile
  job.msea$marfile <- marfile
  job.msea$modfile <- modfile
  job.msea$inffile <- inffile
  job.msea$nperm <- nperm
  job.msea <- ssea.start(job.msea)
  job.msea <- ssea.prepare(job.msea)
  job.msea <- ssea.control(job.msea)
  job.msea <- ssea.analyze(job.msea)
  job.msea <- ssea.finish(job.msea)
}

neuro_painFiles = list.files("inputs/marker_dis/neuro_pain/", full.names = T)
neuro_pain1 = runMSEA(label = "neuro_pain1", marfile = neuro_painFiles[1])
neuro_pain2 = runMSEA(label = "neuro_pain2", marfile = neuro_painFiles[2])
neuro_pain3 = runMSEA(label = "neuro_pain3", marfile = neuro_painFiles[3])
neuro_pain4 = runMSEA(label = "neuro_pain4", marfile = neuro_painFiles[4])
neuro_pain5 = runMSEA(label = "neuro_pain5", marfile = neuro_painFiles[5])
neuro_pain6 = runMSEA(label = "neuro_pain6", marfile = neuro_painFiles[6])
neuro_pain7 = runMSEA(label = "neuro_pain7", marfile = neuro_painFiles[7])
job.metamsea = list()
job.metamsea$job1 = neuro_pain1
job.metamsea$job2 = neuro_pain2
job.metamsea$job3 = neuro_pain3
job.metamsea$job4 = neuro_pain4
job.metamsea$job5 = neuro_pain5
job.metamsea$job6 = neuro_pain6
job.metamsea$job7 = neuro_pain7
job.metamsea = ssea.meta(job.metamsea, "neuro_pain", "Results")
save(neuro_pain1, neuro_pain2, neuro_pain3, neuro_pain4, neuro_pain5,
     neuro_pain6, neuro_pain7, job.metamsea,
     file = "neuro_pain.msea.Rdata")

non_neuro_painFiles = list.files("inputs/marker_dis/non_neuro_pain/", full.names = T)
non_neuro_pain1 = runMSEA(label = "non_neuro_pain1", marfile = non_neuro_painFiles[1])
non_neuro_pain2 = runMSEA(label = "non_neuro_pain2", marfile = non_neuro_painFiles[2])
non_neuro_pain3 = runMSEA(label = "non_neuro_pain3", marfile = non_neuro_painFiles[3])
non_neuro_pain4 = runMSEA(label = "non_neuro_pain4", marfile = non_neuro_painFiles[4])
non_neuro_pain5 = runMSEA(label = "non_neuro_pain5", marfile = non_neuro_painFiles[5])
non_neuro_pain6 = runMSEA(label = "non_neuro_pain6", marfile = non_neuro_painFiles[6])
non_neuro_pain7 = runMSEA(label = "non_neuro_pain7", marfile = non_neuro_painFiles[7])
non_neuro_pain8 = runMSEA(label = "non_neuro_pain8", marfile = non_neuro_painFiles[8])
non_neuro_pain9 = runMSEA(label = "non_neuro_pain9", marfile = non_neuro_painFiles[9])
non_neuro_pain10 = runMSEA(label = "non_neuro_pain10", marfile = non_neuro_painFiles[10])
non_neuro_pain11 = runMSEA(label = "non_neuro_pain11", marfile = non_neuro_painFiles[11])
non_neuro_pain12 = runMSEA(label = "non_neuro_pain12", marfile = non_neuro_painFiles[12])
non_neuro_pain13 = runMSEA(label = "non_neuro_pain13", marfile = non_neuro_painFiles[13])
non_neuro_pain14 = runMSEA(label = "non_neuro_pain14", marfile = non_neuro_painFiles[14])
job.metamsea2 = list()
job.metamsea2$job1 = non_neuro_pain1
job.metamsea2$job2 = non_neuro_pain2
job.metamsea2$job3 = non_neuro_pain3
job.metamsea2$job4 = non_neuro_pain4
job.metamsea2$job5 = non_neuro_pain5
job.metamsea2$job6 = non_neuro_pain6
job.metamsea2$job7 = non_neuro_pain7
job.metamsea2$job8 = non_neuro_pain8
job.metamsea2$job9 = non_neuro_pain9
job.metamsea2$job10 = non_neuro_pain10
job.metamsea2$job11 = non_neuro_pain11
job.metamsea2$job12 = non_neuro_pain12
job.metamsea2$job13 = non_neuro_pain13
job.metamsea2$job14 = non_neuro_pain14
job.metamsea2 = ssea.meta(job.metamsea2, "non_neuro_pain", "Results")
save(non_neuro_pain1, non_neuro_pain2, non_neuro_pain3, non_neuro_pain4, 
     non_neuro_pain5, non_neuro_pain6, non_neuro_pain7, 
     non_neuro_pain8, non_neuro_pain9, non_neuro_pain10, non_neuro_pain11, 
     non_neuro_pain12, non_neuro_pain13, non_neuro_pain14, 
     job.metamsea, file = "non_neuro_pain.msea.Rdata")
 
library(ggvenn)
options(stringsAsFactors = F)
neuro_pain_df = read.table("Results/msea/neuro_pain.results.txt", header = T, sep = "\t") %>% filter(!grepl("_", MODULE))
non_neuro_pain_df = read.table("Results/msea/non_neuro_pain.results.txt", header = T, sep = "\t") %>% filter(!grepl("_", MODULE))

ggvenn(list(
  neuro_pain = neuro_pain_df %>% filter(FDR<0.05) %>% pull(MODULE),
  non_neuro_pain = non_neuro_pain_df %>% filter(FDR<0.05) %>% pull(MODULE)
))

nsig_neuro = neuro_pain_df %>% filter(FDR<0.05) %>% nrow()
nsig_non_neuro = non_neuro_pain_df %>% filter(FDR<0.05) %>% nrow()
nperm = 100000
permut = unlist(lapply(1:nperm, function(x){
  neuro = sample(neuro_pain_df$MODULE, nsig_neuro, replace = F)
  nonneuro = sample(non_neuro_pain_df$MODULE, nsig_non_neuro, replace = F)
  length(intersect(neuro, nonneuro))
}))
length(which(permut>=34))/nperm
34 / mean(permut)

