#!/lustre/home/xdzou/src/anaconda3/bin/python
import argparse
import numpy as np


# -- functions

def pick_outliers(filename):
	fh = open(filename,'r')
	gene2inds = {}
	for line in fh.readlines()[1:]:
		line = line.strip()
		w = line.split("\t")
		gene = w[0].split("|")[0].split(".")[0]
		if gene not in gene2inds:
			gene2inds[gene] = [w[1]]
		else:
			gene2inds[gene].append(w[1])
	fh.close()
	return gene2inds

def pick_controls(filename):
	fh = open(filename,'r')
	gene2controls = {}
	inds = np.array(fh.readline().strip().split("\t")[1:])
	for line in fh.readlines():
		line = line.strip()
		w = line.split("\t")
		gene = w[0].split("|")[0].split(".")[0]
		zscores = w[1:]
		ctrl_idx = []
		for i in range(len(zscores)):
			if zscores[i] != "NA":
				if abs(float(zscores[i]))<1:
					ctrl_idx.append(i)
				else:
					continue
			else:
				continue
		if gene not in gene2controls:
			gene2controls[gene] = list(inds[ctrl_idx])
	fh.close()
	return gene2controls

def extract_gene_body(filename):
	fh = open(filename,'r')
	gene2locs = {}
	for line in fh.readlines():
		line = line.strip()
		w = line.split("\t")
		locs = ":".join(w[0:3])
		gene = w[3].split("|")[0].split(".")[0]
		gene2locs[gene] = locs
	fh.close()
	return gene2locs

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='pick outliers and non-outlier controls')
	parser.add_argument('--outliers',help="a file contains the list of outliers")
	parser.add_argument('--all_medz',help="a file contains medZ matrix")
	parser.add_argument('--gene_loc',help="a file contains locs of gene body")
	parser.add_argument('--out_prefix',help="prefix of output file")

	args = parser.parse_args()

	outliers = pick_outliers(args.outliers)
	print("There are %d aOutlier genes." % (len(outliers)))
	controls = pick_controls(args.all_medz)
	print("Load %d genes with controls." % (len(controls)))
	gb_locs = extract_gene_body(args.gene_loc)
	print("Load gene locs of %d genes." % (len(gb_locs)))
	# iterate on each aOutlier gene and pick controls
	
	for gene in outliers:
		genename = gene.split("|")[0].split(".")[0]
		if genename in gb_locs:
			chrom,start,end = gb_locs[gene].split(":")
		
			outfile = args.out_prefix + "_" + genename + ".txt"
			fho = open(outfile,'w')
			print("%d outliers and %d controls in %s." % (len(outliers[gene]),len(controls[gene]),gene))
			for ind_out in outliers[gene]:
				print("%s\t%s\t%s\t%s\t%s\t%s" % (chrom,start,end,gene,ind_out,"outlier"),file=fho)
			for ind_out in controls[gene]:
				print("%s\t%s\t%s\t%s\t%s\t%s" % (chrom,start,end,gene,ind_out,"control"),file=fho)

			fho.close()
	print("Done!")

