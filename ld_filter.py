import sys, subprocess, logging, os, time
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.DEBUG, 
					format = '%(asctime)s - %(levelname)s - %(message)s')
logging.getLogger().setLevel(logging.INFO)
logging.debug('Close this by uncommenting above line')

plink = "/mnt/f/ResourceData/1000g_plink/plink_1.9/plink"
bedfiles = "/mnt/f/ResourceData/1000g_plink/1000g_eur/"
ld_r_sq_cutoff = 0.5


## file header: chrom\tpos\tname\tref\talt\taf\t-log10p\teffect_size
filename = sys.argv[1]
outname = filename.replace(".summary", "")+".ld_checked.summary"

infile = pd.read_csv(filename, sep = "\t")
infile["snpName"] = infile["chrom"].astype(str)+":"+infile["pos"].astype(str)+":"+infile["ref"].astype(str)+":"+infile["alt"].astype(str)

chroms = set(infile["chrom"])

for chrom in chroms:
	tag = "RUN"
	logging.info(f"[{tag}] RUN FOR CHR{chrom}...")

	tmpDF = infile.query(f"chrom=={chrom}").sort_values(by = "pos")
	tmpDF.snpName.to_csv(f"/tmp/ld_filter_chr{chrom}", index = False, header = False)
	snpnames = tmpDF.snpName.to_list()


	while True:
		logging.info(f"[{tag}] Extract SNPs from chromosome {chrom}...")
		subprocess.run(f"{plink} --bfile {bedfiles}/chr{chrom} --extract /tmp/ld_filter_chr{chrom} --make-bed --out /tmp/ld_filter_chr{chrom}", shell = True, stdout = subprocess.DEVNULL)
		subprocess.run(f"cut -f2 /tmp/ld_filter_chr{chrom}.bim|sort|uniq -d > /tmp/dups", shell = True, stdout = subprocess.DEVNULL)
		if os.path.getsize("/tmp/dups")>0:
			subprocess.run(f"{plink} --bfile {bedfiles}/chr{chrom} --extract /tmp/ld_filter_chr{chrom} --make-bed --exclude /tmp/dups --out /tmp/ld_filter_chr{chrom}", shell = True, stdout = subprocess.DEVNULL)

		logging.info(f"[{tag}] Calculate pair-wise LD...")
		subprocess.run(f"{plink} --bfile /tmp/ld_filter_chr{chrom} --r2 --ld-snp-list /tmp/ld_filter_chr{chrom} --ld-window-r2 0.5 --out /tmp/ld_filter_chr{chrom}", shell = True, stdout = subprocess.DEVNULL)
		check = pd.read_csv(f"/tmp/ld_filter_chr{chrom}.ld", delim_whitespace = True)
		check = check.query('SNP_A!=SNP_B')
		logging.debug(check.shape)

		if check.shape[0]==0:
			break
		else:
			check = check.merge(tmpDF[["snpName", "-log10p"]], left_on = "SNP_A", right_on = "snpName").merge(tmpDF[["snpName", "-log10p"]], left_on = "SNP_B", right_on = "snpName", suffixes = ("_A", "_B"))
			toRemove = set(np.where(check["-log10p_A"]>=check["-log10p_B"], check["snpName_B"], check["snpName_A"]))
			logging.info(f"[{tag}] Totally {len(toRemove)} SNPs to be removed")
			tmpDF = tmpDF[~tmpDF.snpName.isin(toRemove)]
			tmpDF.snpName.to_csv(f"/tmp/ld_filter_chr{chrom}", index = False, header = False)
			snpnames = tmpDF.snpName.to_list()
			tag = "Re-Run"
	logging.info(f"[DONE] Done for chr{chrom}!")
	subprocess.run(f"rm -f /tmp/ld_filter_chr{chrom}*", shell = True)

	tmpDF.to_csv(outname, index = False, sep = "\t", header = (not os.path.exists(outname) or os.path.getsize(outname) == 0), mode = "a")