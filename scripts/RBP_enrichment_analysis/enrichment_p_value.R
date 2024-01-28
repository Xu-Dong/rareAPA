args = commandArgs(trailingOnly=TRUE)
file=args[1]
output=args[2]

con=file(file,'r')
line=readLines(con)
for (i in 1:length(line)){
	peak=as.numeric(strsplit(line[i],"\t")[[1]][3])
	shuffle_peak=as.numeric(strsplit(line[i],"\t")[[1]][4])
	qtl=as.numeric(strsplit(line[i],"\t")[[1]][5])
	matx=rbind(c(peak,qtl-peak),c(shuffle_peak,qtl-shuffle_peak))
	pvalue=fisher.test(matx,alternative="two.sided")$p.value
	odds_ratio=fisher.test(matx,alternative="two.sided")$estimate
	print(pvalue)
	print(odds_ratio)
	new_line=paste(line[i],"\t",as.character(pvalue),"\t",as.character(odds_ratio),sep="")
	write(new_line,output,append=T)
}
close(con)

