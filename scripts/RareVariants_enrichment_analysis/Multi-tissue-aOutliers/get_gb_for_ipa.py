import sys

#load ensembl to gene body

fh = open(sys.argv[1]) # gencode.v26.genes.v8.patched_contigs.bed
en2gb = {}
for line in fh.readlines():
	line = line.strip()
	w = line.split("\t")
	gene = w[3].split(".")[0]
	gb = ":".join(w[0:3])
	en2gb[gene] = gb
fh.close()

# load ensembl to genename

fh = open(sys.argv[2]) #~/data/ref/gencode.v26.geneid2name.txt
genename2ens = {}
for line in fh.readlines():
	line = line.strip()
	w = line.split("\t")
	ensem = w[0].split(".")[0]
	if w[1] not in genename2ens:
		genename2ens[w[1]] = ensem
fh.close()

# main

fh = open(sys.argv[3]) # input/All_IPA_events.ipa_location.bed

for line in fh.readlines():
	line = line.strip()
	w = line.split("\t")
	gene = w[3].split("|")[0]
	if gene in genename2ens:
		if genename2ens[gene] in en2gb:
			gb_chr,gb_s,gb_e = en2gb[genename2ens[gene]].split(":")
			print("%s\t%s\t%s\t%s" % (gb_chr,gb_s,gb_e,w[3]))
		else:
			print("%s\t%s\t%s\t%s" % ("-","-","-",w[3]))
	else:
		print("%s\t%s\t%s\t%s" % ("-","-","-",w[3]))

fh.close()
