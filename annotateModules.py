#!/usr/bin/python
import sys
import pandas as pd
import pandas as np
from scipy.stats import fisher_exact

modfile = open("../coexpress.modules.modfile.tsv", "r")
modfile.readline()
modToGenes = {}
allgenes = set()
for line in modfile:
	module, gene = line.rstrip().split()
	if module in modToGenes:
		modToGenes[module].add(gene)
	else:
		modToGenes[module] = set(gene)
	allgenes.add(gene)
modfile.close()
print(f"Found {len(modToGenes)} modules, including {len(allgenes)} genes.", file = sys.stderr)

kegg = open("c2.cp.kegg.v7.5.1.symbols.gmt", "r")
keggDict = {}
allKEGGgenes = set()
for line in kegg:
	terms = line.rstrip().split()
	name, url = terms[0:2]
	genes = terms[2:]
	keggDict[(name, url)] = set(genes)
	allKEGGgenes |= set(genes)
kegg.close()
print(f"Found {len(keggDict)} canonical pathways from KEGG, including {len(allKEGGgenes)} genes.", file = sys.stderr)
rawKEGGinput = allKEGGgenes.intersection(allgenes)
# print(len(rawKEGGinput))

reactome = open("c2.cp.reactome.v7.5.1.symbols.gmt", "r")
reactomeDict = {}
allReactomegenes = set()
for line in reactome:
	terms = line.rstrip().split()
	name, url = terms[0:2]
	genes = terms[2:]
	reactomeDict[(name, url)] = set(genes)
	allReactomegenes |= set(genes)
reactome.close()
print(f"Found {len(reactomeDict)} canonical pathways from Reactome, including {len(allReactomegenes)} genes.", file = sys.stderr)
rawReactomeinput = allReactomegenes.intersection(allgenes)
# print(len(rawReactomeinput))

module, pathway, url, inHitSymbol, Ninhit, Nin, Ntotalhit, Ntotal = [], [], [], [], [], [], [], []
pathAnns = {
	"module":[], "pathway":[], "url":[], "inHitSymbol":[], 
	"Ninhit":[], "Nin":[], "Ntotalhit":[], "Ntotal":[]
}

def getPW(modKey, pathDict, background):
	for key in pathDict.keys():
		## remove the ones not included in pathDict at all from the module gene list
		ingenes = modToGenes[modKey].intersection(background) 
		## check overlap between the input and the pathway
		overlap = ingenes.intersection(pathDict[key])
		if len(overlap)==0: continue
		## use all genes having kegg annotation as input to define the background
		backgroundhit = background.intersection(pathDict[key])

		pathAnns["module"].append(modKey); pathAnns["pathway"].append(key[0])
		pathAnns["url"].append(key[1]); pathAnns["inHitSymbol"].append(",".join(list(overlap)))
		pathAnns["Ninhit"].append(len(overlap)); pathAnns["Nin"].append(len(ingenes))
		pathAnns["Ntotalhit"].append(len(backgroundhit)); pathAnns["Ntotal"].append(len(background))

flag = 0
for key in modToGenes.keys():
	flag+=1
	print(f"Running annotation for {flag}/{len(modToGenes)}...", file = sys.stderr)
	getPW(modKey=key, pathDict=keggDict, background=rawKEGGinput)
	getPW(modKey=key, pathDict=reactomeDict, background=rawReactomeinput)

print("Annotation finished, now calculating the enrichment and saving to output.", file = sys.stderr)
df = pd.DataFrame(pathAnns)
df2 = df.apply(lambda x: fisher_exact([[x["Ninhit"], x["Nin"]], [x["Ntotalhit"], x["Ntotal"]]], alternative='greater'), axis = 1)
df[['fisher_exact_or','fisher_exact_p']] = pd.DataFrame(df2.tolist(), index= df.index)

df = df.assign(p_adj = lambda x: x["fisher_exact_p"] * (len(modToGenes)*(len(keggDict)+len(reactomeDict)))).query("p_adj<0.05")
df.to_csv("modules.kegg.reactome.ann.tsv", index = None, sep = "\t")