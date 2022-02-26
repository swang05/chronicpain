rm(list = ls())
library(data.table)
library(tidyverse)
library(ggvenn)
library(stringr)

setwd("D://yefang/meta_msea/")

neuro_pain = read.table("Results/msea/neuro_pain.results.txt", header = T, sep = "\t") %>% 
  filter(!grepl("_", MODULE)) %>% filter(FDR<0.05)
nonneuro_pain = read.table("Results/msea/non_neuro_pain.results.txt", header = T, sep = "\t") %>%
  filter(!grepl("_", MODULE)) %>% filter(FDR<0.05)
ann = fread("inputs/pathways/modules.kegg.reactome.ann.tsv", colClasses = c(module ="character")) %>%
  filter(p_adj<0.05) %>% arrange(p_adj)
head(ann)


ggvenn(
  list(NP = neuro_pain$MODULE,
       OA = nonneuro_pain$MODULE)
)
#ggsave("../figures/np_oa_modules_overlap.pdf", width = 12.27, height = 8.69, units = "in")

ann %>%
  filter(module %in% intersect(neuro_pain$MODULE, nonneuro_pain$MODULE)) %>% 
  dplyr::select(module, pathway, fisher_exact_or, p_adj) %>% 
  mutate(resource = ifelse(grepl("KEGG", pathway), "kegg", "reactome")) %>%
  group_by(module, resource) %>% 
  filter(p_adj==min(p_adj)) %>% 
  mutate(pathway = tolower(gsub("KEGG_|REACTOME_", "", pathway))) %>% 
  arrange(p_adj)

ggvenn(list(
  `NP only` = ann %>% filter(module %in% setdiff(neuro_pain$MODULE, nonneuro_pain$MODULE)) %>% pull(pathway) %>% unique(),
  `OA only` = ann %>% filter(module %in% setdiff(nonneuro_pain$MODULE, neuro_pain$MODULE)) %>% pull(pathway) %>% unique())
  )

over.path = tolower(gsub("KEGG_|REACTOME_", "", 
                         intersect(
                           ann %>% filter(module %in% setdiff(neuro_pain$MODULE, nonneuro_pain$MODULE)) %>% pull(pathway) %>% unique(),
                           ann %>% filter(module %in% setdiff(nonneuro_pain$MODULE, neuro_pain$MODULE)) %>% pull(pathway) %>% unique()
                         )))

merge(
  neuro_pain %>% merge(., ann, by.x = "MODULE", by.y = "module"),
  nonneuro_pain %>% merge(., ann, by.x = "MODULE", by.y = "module"),
  by = c("MODULE", "pathway", "DESCR")
) %>% 
  select(MODULE, pathway, DESCR) %>% 
  separate(DESCR, into = c("module", "tissue"), sep = ":") %>% 
  select(-module) %>% group_by(pathway) %>% 
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ","))) %>% 
  rowwise() %>% 
  mutate(Nmod = length(unlist(strsplit(MODULE, ","))),
         Ntissue = length(unlist(strsplit(tissue, ",")))) %>%
  mutate(pathway = tolower(gsub("KEGG_|REACTOME_", "", pathway))) %>% 
  arrange(desc(Nmod)) %>%
  mutate(ifspecific.over = ifelse(pathway %in% over.path, "yes", "no")) %>%
  as.data.frame()

ggvenn(list(
  `NP all` = ann %>% filter(module %in% neuro_pain$MODULE) %>% pull(pathway) %>% unique(),
  `OA all` = ann %>% filter(module %in% nonneuro_pain$MODULE) %>% pull(pathway) %>% unique())
)

modToGenes = fread("inputs/coexpress.modules.modfile.tsv",
                   colClasses = c(MODULE ="character"))
np.path = ann %>% filter(module %in% neuro_pain$MODULE) %>% 
  dplyr::select(pathway, module) %>% 
  merge(., modToGenes, by.x = "module", by.y = "MODULE") %>% 
  group_by(pathway) %>%
  summarise(modules = paste(module, collapse = ","),
            genes = paste(GENE, collapse = ","))
oa.path = ann %>% filter(module %in% nonneuro_pain$MODULE) %>% 
  dplyr::select(pathway, module) %>% 
  merge(., modToGenes, by.x = "module", by.y = "MODULE") %>% 
  group_by(pathway) %>%
  summarise(modules = paste(module, collapse = ","),
            genes = paste(GENE, collapse = ","))

shared = merge(np.path, oa.path, by = "pathway") %>% 
  dplyr::select(-c("modules.x", "modules.y"))

shared$overlaprate = NA
for (i in 1:nrow(shared)){
  print(i)
  np.tmp = unique(unlist(strsplit(shared$genes.x[i], ",")))
  oa.tmp = unique(unlist(strsplit(shared$genes.y[i], ",")))
  n.share = length(intersect(np.tmp, oa.tmp))
  shared$overlaprate[i] = min(n.share / length(np.tmp), n.share / length(oa.tmp))
}

superNets = shared %>% 
  mutate(MODULE = str_pad(1:nrow(.), 3, pad = "0")) %>%
  filter(overlaprate>0.15) %>% 
  unite("NODE", genes.x:genes.y, sep = ",", remove = T) %>%
  separate_rows(NODE, sep = ",") %>% dplyr::select(-overlaprate)
length(table(superNets$MODULE))


library(Mergeomics)
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

superNets %>% 
  dplyr::select(MODULE, NODE) %>%
  write.table("inputs/wkda.module.txt", sep = "\t", row.names = F, quote = F)
superNets %>% 
  mutate(SOURCE = "human",
         DESCR = paste0(MODULE, ":", pathway)) %>%
  dplyr::select(MODULE, SOURCE, DESCR) %>%
  write.table("inputs/wkda.module.info.txt", sep = "\t", row.names = F, quote = F)

runkda <- function(network, modfile, inffile, label, idcon = T){
  if(idcon){
    net = fread(network, colClasses = "character")
    symbols = AnnotationDbi::select(hs, keys = unique(c(net$TAIL, net$HEAD)),
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "ENTREZID") %>% 
      dplyr::filter(!is.na(SYMBOL))
    net %>% merge(., symbols, by.x = "TAIL", by.y = "ENTREZID") %>% 
      mutate(TAIL=SYMBOL) %>% dplyr::select(-SYMBOL) %>%
      merge(., symbols, by.x = "HEAD", by.y = "ENTREZID") %>% 
      mutate(HEAD=SYMBOL) %>% dplyr::select(-SYMBOL) %>% 
      dplyr::select(TAIL, HEAD, WEIGHT) %>% 
      write.table(paste0(network, ".1"), quote = F, sep = "\t", row.names = F)
    network = paste0(network, ".1")
  }
  
  job.kda = list()
  job.kda$label = label
  job.kda$folder = "Results/"
  job.kda$netfile = network
  job.kda$modfile = modfile
  job.kda$inffile = inffile
  job.kda$edgefactor = .5
  job.kda$depth = 1
  job.kda$direction = 0
  
  job.kda = kda.configure(job.kda)
  job.kda = kda.start(job.kda)
  job.kda = kda.prepare(job.kda)
  job.kda = kda.analyze(job.kda)
  job.kda = kda.finish(job.kda)
  
  job.kda = kda2cytoscape(job.kda, node.list = NULL, modules = NULL,
                          ndrivers = 5, depth = 1)
  
  save(job.kda, file = paste0("Results/kda/", label,".Rdata"))
}

modfile = "inputs/wkda.module.txt"
inffile = "inputs/wkda.module.info.txt"

runkda("inputs/nets/giant/amygdala.dat", modfile, inffile, "amygdala")
runkda("inputs/nets/giant/blood.dat", modfile, inffile, "blood")
runkda("inputs/nets/giant/caudate_nucleus.dat", modfile, inffile, "caudate_nucleus")
runkda("inputs/nets/giant/caudate_putamen.dat", modfile, inffile, "caudate_putamen")
runkda("inputs/nets/giant/cerebellar_cortex.dat", modfile, inffile, "cerebellar_cortex")
runkda("inputs/nets/giant/cerebellum.dat", modfile, inffile, "cerebellum")
runkda("inputs/nets/giant/frontal_lobe.dat", modfile, inffile, "frontal_lobe")
runkda("inputs/nets/giant/hippocampus.dat", modfile, inffile, "hippocampus")
runkda("inputs/nets/giant/hypothalamus.dat", modfile, inffile, "hypothalamus")
runkda("inputs/nets/giant/nucleus_accumbens.dat", modfile, inffile, "nucleus_accumbens")
runkda("inputs/nets/giant/spinal_cord.dat", modfile, inffile, "spinal_cord")
runkda("inputs/nets/giant/substantia_nigra.dat", modfile, inffile, "substantia_nigra")
# runkda("inputs/nets/PCNet-V1.3_interactions.txt", modfile, inffile, "PCNet", idcon = F)

###=============================================================###
tissues = c("amygdala", "blood", "caudate_nucleus", "caudate_putamen", "cerebellar_cortex", "cerebellum", "frontal_lobe", "hippocampus", "hypothalamus", "nucleus_accumbens", "spinal_cord", "substantia_nigra")
knownPain = fread("Results/cytoscape/knownPainGenes", header = F)

readF <- function(x){
  tissue = gsub(".tophits.txt", "", basename(x))
  fread(x) %>% mutate(TISSUE = tissue)
}
kds = do.call(rbind,
              lapply(paste0("Results/kda/", tissues, ".tophits.txt"),
                     readF)) %>% filter(FDR<0.05) %>%
  separate(col = "DESCR", into = c("node", "PATHWAY"), sep = ":") %>%
  dplyr::select(NODE, TISSUE, PATHWAY) %>% unique() %>% 
  group_by(NODE) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= ";"))) %>% 
  rowwise() %>%
  mutate(Ntissue = length(unlist(strsplit(TISSUE, ";"))),
         ifknown = ifelse(NODE %in% knownPain$V1, "yes", "no"),
         PATHWAY = tolower(gsub("REACTOME_|KEGG_", "", PATHWAY))) %>% 
  arrange(desc(Ntissue))

kds %>% filter(Ntissue>1) %>% filter(ifknown=="yes")
kds %>% filter(Ntissue>1) %>% pull(NODE)

### build pseduo-modules using KD subnetworks based on PCNet 
nets = fread("inputs/nets/PCNet-V1.3_interactions.txt")
snp2gene = fread("inputs/snps_genes.tsv")

neuroGWASsig = unlist(lapply(list.files("inputs/marker_dis/neuro_pain/", full.names = T),
                             function(x){
                               cat(paste0(x, "\n"))
                               merge(fread(x) %>% filter(VALUE>-log10(5e-8)), snp2gene) %>% pull(GENE)
                             })) %>% unique()

nonneuroGWASsig = unlist(lapply(list.files("inputs/marker_dis/non_neuro_pain/", full.names = T),
                                function(x){
                                  cat(paste0(x, "\n"))
                                  merge(fread(x) %>% filter(VALUE>-log10(5e-8)), snp2gene) %>% pull(GENE)
                                })) %>% unique()


network = nets[(nets$TAIL %in% kds[kds$Ntissue>1,]$NODE &
                  nets$HEAD %in% c(neuroGWASsig, nonneuroGWASsig)) |
                 (nets$HEAD %in% kds[kds$Ntissue>1,]$NODE &
                    nets$TAIL %in% c(neuroGWASsig, nonneuroGWASsig)),]

length(neuroGWASsig)
length(intersect(neuroGWASsig, unique(c(network$TAIL, network$HEAD))))
length(nonneuroGWASsig)
length(intersect(nonneuroGWASsig, c(network$TAIL, network$HEAD)))

network %>% 
  write.table("Results/cytoscape/pcnet.kd.edge.txt", sep = "\t", row.names = F, quote = F)


data.frame(GENE = unique(c(network$TAIL, network$HEAD))) %>%
  mutate(
    whichSig = case_when(
      GENE %in% neuroGWASsig & !GENE %in% nonneuroGWASsig ~ "NP",
      GENE %in% nonneuroGWASsig & !GENE %in% neuroGWASsig ~ "OA",
      GENE %in% nonneuroGWASsig & GENE %in% neuroGWASsig ~ "NP & OA",
      TRUE ~ "na"
    ),
    ifknown = ifelse(GENE %in% knownPain$V1, "yes", "na"),
    label = ifelse(GENE %in% kds[kds$Ntissue>1,]$NODE, "KD", "na")
  ) %>%
  mutate(iffun = ifelse(whichSig=="na" & ifknown=="na" & label=="na",
                        "na", "yes")) %>%
  write.table("Results/cytoscape/pcnet.kd.node.txt", sep = "\t", row.names = F, quote = F)

fread("inputs/pathways/conservative.kds.nodes.kegg.reactome.ann.tsv") %>%
  mutate(pathway = tolower(gsub("KEGG_|REACTOME_", "", pathway))) %>%
  ggplot(., aes(x = reorder(pathway, fisher_exact_or),
                y = log2(fisher_exact_or), fill = -log10(p_adj)))+
  geom_bar(stat = "identity")+coord_flip()+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(y = "Log2(Enrichment fold)", x = "")
#ggsave("../figures/conservative.kds.nodes.kegg.reactome.ann.pdf",
#       width = 12.27, height = 8.69, units = "in")

allgenes = unique(c(nets$HEAD, nets$TAIL))
set.seed(123)
permut <- unlist(lapply(1:10000, function(x){
  tmp = sample(allgenes, size = 24)
  length(intersect(tmp, knownPain$V1))
}))
length(which(permut>=2))/10000
2/mean(permut)

####========================================########
## infer the directions
snp2gene = fread("inputs/snps_genes.tsv") %>% 
  filter(GENE %in% kds[kds$Ntissue>1,]$NODE)
table(snp2gene$GENE)

readSNPs <- function(files){
  out = fread(files[1]) %>% 
    mutate(study = strsplit(basename(files[1]), "\\.")[[1]][1])
  for (i in 2:length(files)){
    out = rbind(
      out,
      fread(files[i]) %>% 
        mutate(study = strsplit(basename(files[i]), "\\.")[[1]][1])
    )
  }
  return(out)
}
neuropainSNPs = readSNPs(
  list.files("F://projects/pain/neuro_pain/", 
             recursive = T, full.names = T, 
             pattern = "ld_checked")
)
nonneuropainSNPs = readSNPs(
  list.files("F://projects/pain/non_neuro_pain/", 
             recursive = T, full.names = T, 
             pattern = "ld_checked")
)
rbind(
  merge(snp2gene, 
        neuropainSNPs %>% filter(`-log10p`>3), 
        by.x = "MARKER", by.y = "snpName") %>% mutate(label = "NP"),
  merge(snp2gene, 
        nonneuropainSNPs %>% filter(`-log10p`>3), 
        by.x = "MARKER", by.y = "snpName") %>% mutate(label = "OA")
) %>% 
  dplyr::select(GENE, MARKER, effect_size, label) %>% 
  mutate(effect_size = sign(effect_size)) %>% 
  group_by(GENE, label) %>% 
  mutate(effect_size = sum(effect_size), n = n()) %>% 
  filter(effect_size+n==0) %>% 
  dplyr::select(GENE, label, effect_size) %>% unique() %>%
  reshape2::dcast(GENE ~ label)


##=============================================================##
columns = fread("f://ResourceData/brainspan/columns_metadata.csv")
rows = fread("f://ResourceData/brainspan/rows_metadata.csv")
expr = fread("f://ResourceData/brainspan/expression_matrix.csv")
head(expr)
rows %>% filter(gene_symbol %in% kds[kds$Ntissue>1,]$NODE)

head (columns)
cols = columns %>% 
  mutate(stage = case_when(
    age %in% c("4 mos") ~ 6,
    age %in% c("10 mos", "1 yrs") ~ 7,
    age %in% c("2 yrs", "3 yrs", "4 yrs") ~ 8,
    age %in% c("8 yrs", "11 yrs") ~ 9,
    age %in% c("13 yrs", "15 yrs", "18 yrs", "19 yrs") ~ 10,
    age %in% c("21 yrs", "23 yrs", "30 yrs",
               "36 yrs", "37 yrs", "40 yrs") ~ 11,
    TRUE ~ 0
  )) %>% 
  filter(stage>0) %>%
  filter(structure_acronym %in% c("AMY", "CBC", "CB", "HIP", "MFC", "DFC", "M1C", "VFC", "STR"))

cols$expr_col = colnames(expr)[cols$column_num+1]
int.expr = expr[rows[rows$gene_symbol %in% kds[kds$Ntissue>1,]$NODE,]$row_num,
                c(1, (cols$column_num+1)), with=FALSE]

merge(rows[,c("gene_symbol", "row_num")],
      int.expr, by.x = "row_num", by.y = "V1") %>% 
  dplyr::select(-row_num) %>%
  t() %>%
  as.data.frame(stringsAsFactors = F) %>% 
  rownames_to_column("gene_symbol") %>%
  `colnames<-`(.[1,]) %>%
  .[-1,] %>%
  `rownames<-`(NULL) %>% 
  merge(cols[,c("expr_col", "stage", "structure_acronym", "gender")], ., 
        by.x = "expr_col", by.y = "gene_symbol") %>% 
  mutate_at(c(5:27), as.numeric) %>% 
  dplyr::select(-expr_col) %>%
  group_by(stage, structure_acronym, gender) %>% 
  summarise_each(funs(mean(., na.rm=TRUE))) %>% 
  reshape2::melt(., id.vars = 1:3) %>% 
  filter(structure_acronym=="HIP") %>%
  ggplot(., aes(x = stage, y = log2(value),
                color = gender))+
  geom_point()+facet_wrap( ~ variable)+geom_line()

##=============================================================##
library(readxl)
expr = rbind(
  read_xlsx("../public_data/whole_blood_cfa_deg.xlsx") %>% mutate(label = "wholeBlood_cfa"),
  read_xlsx("../public_data/whole_blood_sni_deg.xlsx") %>% mutate(label = "wholeBlood_sni"),
  read_xlsx("../public_data/whole_brain_cfa_deg.xlsx") %>% mutate(label = "wholeBrain_cfa"),
  read_xlsx("../public_data/whole_brain_sni_deg.xlsx") %>% mutate(label = "wholeBrain_sni"),
  read_xlsx("../public_data/dorsal_root_ganglia_cfa_deg.xlsx") %>% mutate(label = "dorsalRootGanglia_cfa"),
  read_xlsx("../public_data/dorsal_root_ganglia_sni_deg.xlsx") %>% mutate(label = "dorsalRootGanglia_sni"),
  read_xlsx("../public_data/spinal_cord_cfa_deg.xlsx") %>% mutate(label = "spinalCord_cfa"),
  read_xlsx("../public_data/spinal_cord_sni_deg.xlsx") %>% mutate(label = "spinalCord_sni")
) %>%
  filter(toupper(gene_symbol) %in% kds[kds$Ntissue>1,]$NODE) %>% 
  #filter(toupper(gene_symbol) %in% knownPain$V1) %>%
  dplyr::select(gene_symbol, label, value_CTR, value_CFA)
forplot = expr %>% 
  separate(col = "label", into = c("region", "model"), sep = "_") %>%
  reshape2::melt(id.vars = 1:3) %>%
  mutate(label = paste0(region, "_", variable)) %>% 
  dplyr::select(gene_symbol, model, label, value) %>%
  reshape2::dcast(gene_symbol+model ~ label) %>% 
  mutate(wholeBrain = (wholeBrain_value_CFA/wholeBrain_value_CTR)/(wholeBlood_value_CFA/wholeBlood_value_CTR),
         dorsalRootGanglia = (dorsalRootGanglia_value_CFA/dorsalRootGanglia_value_CTR)/(wholeBlood_value_CFA/wholeBlood_value_CTR),
         spinalCord = (spinalCord_value_CFA/spinalCord_value_CTR)/(wholeBlood_value_CFA/wholeBlood_value_CTR)
  ) %>% 
  dplyr::select(gene_symbol, model, 
                wholeBrain, dorsalRootGanglia, spinalCord) %>%
  reshape2::melt(id.vars = 1:2) %>% 
  dplyr::rename(enrichFactor=value) %>%
  mutate(Ntissue = unlist(lapply(gene_symbol, function(x){
    kds[kds$NODE==toupper(x),]$Ntissue
  }))) %>%
  filter(!is.na(log(enrichFactor))&!is.infinite(log(enrichFactor)))

orders = forplot %>% 
  filter(model=="cfa" & variable=="wholeBrain") %>% 
  arrange(desc(enrichFactor)) %>% pull(gene_symbol)

forplot %>% 
  mutate(gene_symbol = factor(gene_symbol,
                              levels = rev(c(orders, setdiff(.$gene_symbol, orders)))),
         variable = factor(variable, 
                           levels = c("wholeBrain", "spinalCord", "dorsalRootGanglia"))) %>%
  ggplot(., aes(x = variable, y = gene_symbol, 
                fill = log(enrichFactor)))+
  geom_tile(color = "black")+facet_wrap(~model)+
  scale_fill_gradient2(low = "blue", high = "red")+
  theme_bw()+labs(x = "", y = "")
#ggsave("../figures/caf_sni_mouse_models_expr.pdf",
#       width = 12.27, height = 8.69, units = "in")







