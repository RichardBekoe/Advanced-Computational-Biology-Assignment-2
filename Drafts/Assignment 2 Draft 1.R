
# Q1) Section 1 ---------------------------------------------------------------------

ceu = t(read.table(file.choose()))
jpt = t(read.table(file.choose()))
yri = t(read.table(file.choose()))
tsi = t(read.table(file.choose()))

PopX = t(read.table(file.choose()))


#  a) ---------------------------------------------------------------------

num.snps=dim(ceu)[2]
table_heterozygous=rep(NA,num.snps)

for (i in 1:num.snps) {
  count <- rowsum(ceu[,i], as.integer(gl(nrow(ceu), 2, nrow(ceu))))
  table_heterozygous[i] <- (sum(count == 1) / length(count))
}
ceu_median_heterozygosity <- median(table_heterozygous)



num.snps=dim(jpt)[2]
table_heterozygous=rep(NA,num.snps)

for (i in 1:num.snps) {
  count <- rowsum(jpt[,i], as.integer(gl(nrow(jpt), 2, nrow(jpt))))
  table_heterozygous[i] <- (sum(count == 1) / length(count))
}
jpt_median_heterozygosity <- median(table_heterozygous)


num.snps=dim(tsi)[2]
table_heterozygous=rep(NA,num.snps)

for (i in 1:num.snps) {
  count <- rowsum(tsi[,i], as.integer(gl(nrow(tsi), 2, nrow(tsi))))
  table_heterozygous[i] <- (sum(count == 1) / length(count))
}
tsi_median_heterozygosity <- median(table_heterozygous)


num.snps=dim(yri)[2]
table_heterozygous=rep(NA,num.snps)

for (i in 1:num.snps) {
  count <- rowsum(yri[,i], as.integer(gl(nrow(yri), 2, nrow(yri))))
  table_heterozygous[i] <- (sum(count == 1) / length(count))
}
yri_median_heterozygosity <- median(table_heterozygous)



num.snps=dim(PopX)[2]
table_heterozygous=rep(NA,num.snps)

for (i in 1:num.snps) {
  count <- rowsum(PopX[,i], as.integer(gl(nrow(PopX), 2, nrow(PopX))))
  table_heterozygous[i] <- (sum(count == 1) / length(count))
}
PopX_median_heterozygosity <- median(table_heterozygous)


# # adds every other 2nd row and puts the result into a matrix i.e. heterozygosity = 1
# count <- rowsum(ceu[,1], as.integer(gl(nrow(ceu), 2, nrow(ceu))))
# # the number of heterozygous / total number of participants (2 rows per participant)
# proportion_heterozygous <- (sum(count == 1) / length(count))



# # creation of a loop
# heter.calc <- function(pop=ceu)
# {
# num.snps=dim(ceu)[2]
# table_heterozygous=rep(NA,num.snps)
# for (i in 1:num.snps) {
#    count <- rowsum(ceu[,i], as.integer(gl(nrow(ceu), 2, nrow(ceu))))
#    table_heterozygous[i] <- (sum(count == 1) / length(count))
#  }
#    median_heterozygosity <- median(table_heterozygous)
#  }
# # call function
#  ceu_median_heter=heter.calc(pop=ceu)
#  tsi_median_heter=heter.calc(pop=tsi)
#  jpt_median_heter=heter.calc(pop=jpt)
#  yri_median_heter=heter.calc(pop=yri)
#  PopX_median_heter=heter.calc(pop=PopX)


# b) ----------------------------------------------------------------------


fst.calc=function(x,y,z,A,B)
{
  all=c(x,y,z,A,B)
  p=length(all[all==0])/length(all)
  pX=length(x[x==0])/length(x)
  pY=length(y[y==0])/length(y)
  pZ=length(z[z==0])/length(z)
  pA=length(A[A==0])/length(A)
  pB=length(B[B==0])/length(B)

  
  cX=length(x)/length(all)
  cY=length(y)/length(all)
  cZ=length(z)/length(all)
  cA=length(A)/length(all)
  cB=length(B)/length(all)
  Fst=(sum(c(cX,cY,cZ,cA,cB)*(c(pX,pY,pZ,cA,cB)-p)^2))/(p*(1-p))
  return(Fst)
}

num.snps=dim(ceu)[2]
fst.snps=rep(NA,num.snps)
for (i in 1:num.snps) 
  fst.snps[i]=fst.calc(ceu[,i],jpt[,i],yri[,i],tsi[,i],PopX[,i])

median_fst <- median(fst.snps,na.rm=TRUE)

# summary(fst.snps)
 

# c)


# d)

ceu_five_snp <- ceu[1:4,1:5]

jpt_five_snp <- jpt[1:4,1:5]

yri_five_snp <- yri[1:4,1:5]

tsi_five_snp <- tsi[1:4,1:5]

PopX_five_snp <- PopX[1:4, 1:5]


# e)



# Section 2 ---------------------------------------------------------------

ldsel = function(ngen=500,nall=2000,init=c(0.5,0,0,0.5),rho=0,s=0,mu=0){
  hp = matrix(0,ngen,4)
  p = rep(0,4)
  hp[1,] = init/sum(init)
  ldstat = matrix(0,ngen,2)
  p1 = hp[1,1]+hp[1,2]		# allele prop at loc 1
  p2 = hp[1,1]+hp[1,3]		# allele prop at loc 2
  for(i in 2:ngen){
    # haplotype proportions in next generation after recombination:
    pp = hp[i-1,]*(1-rho)+rho*c(p1*p2,p1*(1-p2),(1-p1)*p2,(1-p1)*(1-p2))
    # now let's have some mutation and selection:
    p[1] = sum(pp*c((1-mu)^2,mu*(1-mu),mu*(1-mu),mu^2))*(1+s)
    p[2] = sum(pp*c(mu*(1-mu),(1-mu)^2,mu^2,mu*(1-mu)))*(1+s)
    p[3] = sum(pp*c(mu*(1-mu),mu^2,(1-mu)^2,mu*(1-mu)))
    p[4] = sum(pp*c(mu^2,mu*(1-mu),mu*(1-mu),(1-mu)^2))
    # sample haplotypes in next generation and record counts
    tmp = sample(1:4,nall,repl=T,prob=p)
    for (j in 1:4) hp[i,j]=sum(tmp==j)/nall
    p1 = hp[i,1]+hp[i,2]		# allele prop at loc 1
    p2 = hp[i,1]+hp[i,3]		# allele prop at loc 2
    # compute D'_00
    D00 = hp[i,1]-p1*p2
    if(D00>0) ldstat[i,1] = D00/min(p1*(1-p2),p2*(1-p1))
    if(D00<=0) ldstat[i,1] = -D00/min(p1*p2,(1-p2)*(1-p1))
    # compute r^2_00
    ldstat[i,2] = D00^2/(p1*p2*(1-p1)*(1-p2))
  }
  pA=apply(hp[,1:2],1,sum)
  pB=apply(hp[,c(1,3)],1,sum)
  to.plot='yes'
  if (to.plot=='yes')
  {
    lwd.val=2
    par(mfrow=c(2,1),mar=c(3,2,2,1))
    matplot(hp[,1:4],type="l",xlim=c(0,ngen*1.2),ylim=c(0,1),lty=1,lwd=lwd.val,xlab="",ylab="haplotype proportion",main=paste("LD: rho=",rho,", s=",s,", mu=",mu,sep=''))
    lines(1:dim(hp)[1],pA,lty=2,lwd=lwd.val,col=5)
    lines(1:dim(hp)[1],pB,lty=3,lwd=lwd.val,col=6)
    legend(ngen,1,leg=c("AB","Ab","aB","ab","A","B"),lty=c(1,1,1,1,2,3),col=1:6,lwd=2,cex=1)
    ylim.plot=c(0,max(c(ldstat)))
    if (sum(is.na(ylim.plot))!=0) ylim.plot=c(0,1)
    matplot(ldstat[-1,],type="l",xlim=c(0,ngen*1.2),ylim=ylim.plot,lty=1,lwd=lwd.val,col=3:4,xlab="",ylab="")
    legend(ngen,max(ylim.plot),leg=c("|D'|","r^2"),lty=1,col=3:4,lwd=2,cex=1.1)
  }
  return(c(hp[ngen,],ldstat[ngen,]))
}

ldsel(ngen=5000,nall=2000,init=c(0.5,0.0,0.0,0.5),rho=0,s=0,mu=0)

# (i.e. npop x init / nall)





