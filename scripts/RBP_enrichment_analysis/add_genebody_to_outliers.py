import sys

ref2gb = {}


fh = open(sys.argv[1]) #hg38_Refseq_Genebody.medz_2.apa.bed
for line in fh.readlines():
	line = line.strip()
	w = line.split(" ")
	ref2gb[w[3]] = "\t".join(w[0:3])
fh.close()


fh = open(sys.argv[2])# apa_medz.picked.Z_2.rm_global_chrX_Y.txt
for line in fh.readlines()[1:]:
	line = line.strip()
	w = line.split("\t")
	refid = w[0].split("|")[0].split(".")[0]

	if refid in ref2gb:
		print("%s\t%s" % (ref2gb[refid],line))
fh.close()
