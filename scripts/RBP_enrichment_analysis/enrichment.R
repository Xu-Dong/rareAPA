args = commandArgs(trailingOnly=TRUE)

library(dplyr)

dat <- read.table(args[1],header=F,sep="\t")
names(dat) <- c("RBP","peak","shuffle","total")
cat("total number of RBPs processed: ",dim(dat)[1],"\n")

dat.f <- dat %>% filter(peak+shuffle>5)
cat("Number of RBPs passed filtering: ",dim(dat.f)[1],"\n")

OR <- c()
Low <- c()
High <- c()
Pvalue <- c()
for(i in 1:dim(dat.f)[1]){
	a <- dat.f[i,2]
	b <- dat.f[i,3]
	total <- dat.f[i,4]
	c <- total - a
	d <- total - b
	ft <- fisher.test(matrix(c(a,b,c,d),ncol=2))
	OR[i] <- ft$estimate
	Low[i] <- ft$conf.int[1]
	High[i] <- ft$conf.int[2]
	Pvalue[i] <- ft$p.value
}

outdf <- data.frame(RBP = dat.f$RBP,odds_ratio=OR,Low_CI=Low,Up_CI=High,P=Pvalue)

write.table(outdf,file=args[2],quote=F,sep="\t",row.names=F)

