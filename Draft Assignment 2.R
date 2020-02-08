
# Q1) ---------------------------------------------------------------------

ceu = t(read.table(file.choose()))
jpt = t(read.table(file.choose()))
yri = t(read.table(file.choose()))
tsi = t(read.table(file.choose()))

# a)

rowsum(df1[,1], as.integer(gl(nrow(df1), 2, nrow(df1))))


#  b)

median(...,na.rm=TRUE)

# between every pairing; 5 populations;
# calculate the median Fst across all SNPs

# each seperate 3 times loop; then calculation of median

# method?
# calculate the Fst for each population?
# find the median Fst for across all populations

# DO - Calculate the Fst, between every pairing, for each column (SNP)
# Then find the median using all SNPs of that population
# Should output 5 medians

# something similar
expected, why, why not?

num.snps=dim(ceu)[2]
fst.snps=rep(NA,num.snps)
for (i in 1:num.snps) 
  fst.snps[i]=fst.calc(ceu[,i],jpt[,i],yri[,i])

summary(fst.snps)



# c)
# where do you think PopX comes from?

# could Popx be a combination of the other populations?



# d)
# first five columns (SNPs); first 4 rows (first two individuals)
# Buld an Ancestral Recombination Graph (ARG) 5 SNP region

DO - research Ancestral Recombination Graph (ARG)

ceu[1:4,1:5]
jpt[1:4,1:5]
yri[1:4,1:5]
tsi[1:4,1:5]


# e)

# Arg as expected, why or why not?
# ten sequences, 5 diploid, time until each of the coalescent events?
# all individuals coalesce?


DO - research time to coalesce (each / all) ?









