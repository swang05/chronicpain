#!/usr/bin/python3
import sys, gzip, os
import pandas as pd

vcfName = sys.argv[1]
infile = gzip.open(vcfName, "rt")
tmpFileName = "tmp"
tmpFile = open(tmpFileName, "w")
print("chrom\tpos\tname\tref\talt\taf\t-log10p\teffect_size", file = tmpFile)
nowchrom = ""
for line in infile:
	if line.startswith("#"): continue
	terms = line.rstrip().split("\t")
	chrom, pos, name, ref, alt = terms[0:5]
	if not chrom==nowchrom:
		print(f"Now scanning on Chr{chrom}...", file = sys.stderr)
		nowchrom = chrom

	try:
		af = float(terms[7].replace("AF=", ""))
		# Required: minor allele frequency > 0.05
		if af<=0.05 or af>=0.95: continue
	except:
		af = "NA"

	vals = dict(zip(terms[-2].split(":"), terms[-1].split(":")))
	if not "AF" in vals: vals["AF"]=af
	
	# Required: imputation quality > 0.3
	if "SI" in vals and not float(vals["SI"])>0.3: continue

	print(chrom, pos, name, ref, alt, vals["AF"], vals["LP"], vals["ES"], sep = "\t", file = tmpFile)
tmpFile.close()

snps = pd.read_csv(tmpFileName, sep = "\t")
snps = snps.query('not (chrom=="X" or chrom=="Y" or chrom=="MT")')
snps = snps.sort_values(by = ["-log10p"], ascending = False)

snps = snps.drop_duplicates(subset = ["chrom", "pos", "ref", "alt"], keep = False)

snps.to_csv(vcfName.replace(".vcf.gz", "")+".summary", index = False, sep = "\t")

#ranked in the bottom 50% (weaker association) based on disease association p-values 
snps = snps.head(n=int(snps.shape[0]/2))

snps.to_csv(vcfName.replace(".vcf.gz", "")+".top50perc.summary", index = False, sep = "\t")
os.remove(tmpFileName)