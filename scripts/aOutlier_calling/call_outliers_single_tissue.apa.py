#!/usr/bin/env python

## pick single tissue outliers by parsing flat file
## for each tissue-gene pair, we pick individuals with |Z| >= 2
## for each tissue also make file where the most extreme individual is picked with no threshold requirement
## also for each tissue make file with matrix of which individuals were tested for which genes
## these last two sets of files will be used by the control picking script to make features

## outputs outliers as tab-delimited file with the following columns:
## GENE INDS TISSUE Z

# modified by @zouxd,2021-08-31
import os

dir = "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA"
infile = dir + '/output/gtex_aOutlier_v8_normalized_pdui.peer.txt'
outfile = dir + '/output/apa/aOutliers_singlez_picked.Z3.txt'
#indivfile = dir + '/output/apa/apa_medz_1.picked_count_per_inds.txt' # @xdzou 2023-1-24
noThreshPrefix = dir + '/output/apa/single_tissue_z1/aOutliers_singlez_nothreshold_'


zscores = open(infile, 'r')
outliers = open(outfile, 'w')
#indivCounts = open(indivfile,'r')
nothresh = dict() # for file handle; make these on the fly

# get list of individuals to keep (remove global outliers)
#print "Getting list of individuals to keep..."
#indivsToKeep = []
#for indivline in indivCounts.readlines()[1:]:
#	indivline = indivline.strip()
#	w = indivline.split("\t")
#	if w[-1] == "NO":
#		indivsToKeep.append(w[0])
#indivCounts.close()

header = zscores.readline().strip().split()
individuals = header[2:]

# restrict to individuals in keep list
#keepInd = [ind in indivsToKeep for ind in individuals]
#individuals = [ind for keep,ind in zip(keepInd,individuals) if keep]

outliers.write('GENE\tINDS\tTISSUE\tZ\n')

## function to create file handles for a new tissue
## adds them to the nothresh dict
## note this function relies on pre-existing variables
## modifies the nothresh dict directly
def addHandles(tissuename):
    countFile = noThreshPrefix + tissuename + "_counts.txt"
    count = open(countFile, 'w')
    count.write('GENE\t' + '\t'.join(individuals) + "\n")
    nothresh[tissuename] = count

# each line is a tissue-gene pair. iterate through them
print "Identifying single-tissue outliers..."
for line in zscores:
    line = line.strip().split()
    tissue = line[0]
    gene = line[1]
    line = line[2:]
    #line = [val for keep, val in zip(keepInd,line[2:]) if keep]
    # get list of individuals in that tissue and then replace NA by 0
    countstring = gene + '\t' + '\t'.join(['0' if i == "NA" else '1' for i in line]) + '\n'
    values = [0 if i == "NA" else float(i) for i in line]

    # get all individuals with |z|>2/3/4/5
    out_idx = [(ind,z) for ind,z in zip(individuals,values) if abs(z)>3]
    if len(out_idx)>0:
        for inds,zvalue in out_idx:
            outstring = '\t'.join([gene, inds, tissue, str(zvalue)]) + '\n'
# write to combined outlier file if the individual passes the Z threshold
            outliers.write(outstring)

    # write to files specific for this tissue
    if tissue not in nothresh:
        addHandles(tissue)
    count = nothresh[tissue]
    count.write(countstring)

zscores.close()
outliers.close()
for handle in nothresh.values():
    handle.close()

print "Done!"
