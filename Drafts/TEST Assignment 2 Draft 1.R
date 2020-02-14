
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
  pX=length(ceu[,1][ceu[,1]==0])/length(ceu[,1])
  pY=length(y[y==0])/length(y)
  pZ=length(z[z==0])/length(z)
  pA=length(A[A==0])/length(A)
  pB=length(B[B==0])/length(B)

  
  cX=length(x)/length(all)
  cY=length(y)/length(all)
  cZ=length(z)/length(all)
  cA=length(A)/length(all)
  cB=length(B)/length(all)
  Fst=(sum(c(cX,cY,cZ,cA,cB)*(c(pX,pY,pZ,pA,pB)-p)^2))/(p*(1-p))
  return(Fst)
}

num.snps=dim(ceu)[2]
fst.snps=rep(NA,num.snps)
for (i in 1:num.snps) {
  fst.snps[i]=fst.calc(ceu[,i],jpt[,i],yri[,i],tsi[,i],PopX[,i])
}
median_fst <- median(fst.snps,na.rm=TRUE)

median_fst
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
# FUNCTION DEFINITION ONLY
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
  to.plot='no'
  if (to.plot=='yes')
    
  {
    lwd.val=2
    par(mfrow=c(2,1),mar=c(3,2,2,1))
    matplot(hp[,1:4],type="l",xlim=c(0,ngen*1.2),ylim=c(0,1),lty=1,lwd=lwd.val,xlab="",
            ylab="haplotype proportion",main=paste("LD: rho=",rho,", s=",s,", mu=",mu,sep=''))
    lines(1:dim(hp)[1],pA,lty=2,lwd=lwd.val,col=5)
    lines(1:dim(hp)[1],pB,lty=3,lwd=lwd.val,col=6)
    legend(ngen,1,leg=c("AB","Ab","aB","ab","A","B"),lty=c(1,1,1,1,2,3),col=1:6,lwd=2,cex=1)
    
    ylim.plot=c(0,max(c(ldstat)))
    if (sum(is.na(ylim.plot))!=0) ylim.plot=c(0,1)
    matplot(ldstat[-1,],type="l",xlim=c(0,ngen*1.2),ylim=ylim.plot,lty=1,lwd=lwd.val,col=3:4,xlab="",ylab="")
    legend(ngen,max(ylim.plot),leg=c("|D'|","r^2"),lty=1,col=3:4,lwd=2,cex=1.1)
  }
 return(c(hp[ngen,],ldstat[ngen,], which.max(pA)))
 
}

# USING MAX AS pA

maxA_values.df <- data.frame(maxA_values = NA)

# CALLING THE FUNCTION
# Saving proportions and maxA in environment as 'Output' as a numberic list

npop=1000

for (i in 1:npop) { 
  
  Output <- ldsel(ngen=5000,nall=2000,init=c(0.5,0.0,0.0,0.5),rho=0,s=0,mu=0)
  # 1    0    0    0  NaN  NaN 1544 (example output)
  # AB   Ab   aB   ab  A    B   maxA
  if (Output[1] == 1) { maxA_values.df <- rbind(maxA_values.df, Output[7])
  }
  # If Output[1] (AB) equals to 1 (fixation) then add the Output[7] (the max(pA) value to the data.frame maxA_values.df)
}

# Creation of a loop to simulate every iteration of a population (npop); 
# as you can only simulate one population at a time otherwise

median(maxA_values.df$maxA_values,na.rm=TRUE)
# find the median of the values in maxA_values.df hence median fixation time
#  Or median(maxA_values.df[,1],na.rm=TRUE)
# npop = 10,000 is takes too long, tried a smaller number for npop 1000 

# Labels <- c("AB","Ab","aB","ab","A","B","maxA")
# Output_dataframe <- data.frame(Output, Labels) [was going to create a dataframe with labels]



# b) ----------------------------------------------------------------------

maxA_values.df <- data.frame(maxA_values = NA)

# adding some recombination
npop=1000

for (i in 1:npop) { 
  Output <- ldsel(ngen=5000,nall=2000,init=c(0.5,0.0,0.0,0.5),rho=0.01,s=0,mu=0)
    if (Output[1] == 1) { maxA_values.df <- rbind(maxA_values.df, Output[7])
  }
}

median(maxA_values.df$maxA_values,na.rm=TRUE)


# RERUN WITH NEW OUTPUT DERIVED FROM AB INSTEAD OF FROM ALLELE A FREQUENCY
# c) ----------------------------------------------------------------------


# ("AB","Ab","aB","ab")
# setting values so that the selected allele A has initial frequency 0.01

# Varying s, what is the pattern over time? For one value of s, find the median
# time to fixation of the A allele.

maxA_values.df <- data.frame(maxA_values = NA)

npop=1000

for (i in 1:npop) { 
  Output <- ldsel(ngen=5000,nall=2000,init=c(0.005,0.005,0.495,0.495),rho=0.0,s=0.2,mu=0)
  if (Output[1] == 1) { maxA_values.df <- rbind(maxA_values.df, Output[7])
  }
}

median(maxA_values.df$maxA_values,na.rm=TRUE)


# d) ----------------------------------------------------------------------


# Now increase rho and repeat (c). How do patterns change?
# setting values so that the selected allele A has initial frequency 0.01
# ("AB","Ab","aB","ab")

# Varying s, what is the pattern over time? For one value of s, find the median
# time to fixation of the A allele.

maxA_values.df <- data.frame(maxA_values = NA)

npop=1000

for (i in 1:npop) { 
  Output <- ldsel(ngen=5000,nall=2000,init=c(0.005,0.005,0.495,0.495),rho=0.3,s=0.2,mu=0)
  if (Output[1] == 1) { maxA_values.df <- rbind(maxA_values.df, Output[7])
  }
}

median(maxA_values.df$maxA_values,na.rm=TRUE)






# MAX AB OUTPUT; instead of A allele frequency--------------------------------------------------------------------

# FUNCTION DEFINITION ONLY
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
    matplot(hp[,1:4],type="l",xlim=c(0,ngen*1.2),ylim=c(0,1),lty=1,lwd=lwd.val,xlab="",
            ylab="haplotype proportion",main=paste("LD: rho=",rho,", s=",s,", mu=",mu,sep=''))
    lines(1:dim(hp)[1],pA,lty=2,lwd=lwd.val,col=5)
    lines(1:dim(hp)[1],pB,lty=3,lwd=lwd.val,col=6)
    legend(ngen,1,leg=c("AB","Ab","aB","ab","A","B"),lty=c(1,1,1,1,2,3),col=1:6,lwd=2,cex=1)
    
    ylim.plot=c(0,max(c(ldstat)))
    if (sum(is.na(ylim.plot))!=0) ylim.plot=c(0,1)
    matplot(ldstat[-1,],type="l",xlim=c(0,ngen*1.2),ylim=ylim.plot,lty=1,lwd=lwd.val,col=3:4,xlab="",ylab="")
    legend(ngen,max(ylim.plot),leg=c("|D'|","r^2"),lty=1,col=3:4,lwd=2,cex=1.1)
  }
  return(c(hp[ngen,],ldstat[ngen,], which.max(hp[,1:1])))
}


maxA_values.df <- data.frame(maxA_values = NA)

# c (MAX AB OUTPUT; SAME CODE DIFFERENT FUNCTION OUTPUT) -----------------------------------------------------------------------

# setting values so that the selected allele A has initial frequency 0.01

# Varying s, what is the pattern over time? For one value of s,
# find the median time to fixation of the A allele.
# ("AB","Ab","aB","ab")

maxA_values.df <- data.frame(maxA_values = NA)

npop=1000

for (i in 1:npop) { 
  Output <- ldsel(ngen=5000,nall=2000,init=c(0.005,0.005,0.495,0.495),rho=0.0,s=0.02,mu=0)
  if (Output[1] == 1) { maxA_values.df <- rbind(maxA_values.df, Output[7])
  }
}

median(maxA_values.df$maxA_values,na.rm=TRUE)


# d) ----------------------------------------------------------------------

# Now increase rho and repeat (c). How do patterns change?

# setting values so that the selected allele A has initial frequency 0.01

# Varying s, what is the pattern over time? For one value of s, find the median
# time to fixation of the A allele.
# ("AB","Ab","aB","ab")

maxA_values.df <- data.frame(maxA_values = NA)

npop=1000

for (i in 1:npop) { 
  Output <- ldsel(ngen=5000,nall=2000,init=c(0.005,0.005,0.495,0.495),rho=0.3,s=0.001,mu=0)
  if (Output[1] == 1) { maxA_values.df <- rbind(maxA_values.df, Output[7])
  }
}

median(maxA_values.df$maxA_values,na.rm=TRUE)

