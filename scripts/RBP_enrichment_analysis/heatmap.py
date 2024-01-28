import sys
import pandas

def make_heatmap_table(enrichment,output):
	data=pandas.read_table(enrichment,header=None)
	pvalue=[]
	genes=[]
	tissues=list(data[1][0:27])
	for i in range(162):
		a=list(data[7][27*i:27*(i+1)])
		gene=data[0][27*i]
		pvalue.append(a)
		genes.append(gene)
	df=pandas.DataFrame(pvalue,index=genes,columns=tissues)
	df.to_csv(output,sep='\t')

if __name__ =='__main__':
	enrichment=sys.argv[1]
	output=sys.argv[2]
	make_heatmap_table(enrichment,output)
