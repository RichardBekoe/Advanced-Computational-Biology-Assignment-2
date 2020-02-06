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
  	legend(ngen,1,leg=c("AB","Ab","aB","ab","A","B"),lty=c(1,1,1,1,2,3),col=1:6,lwd=2,cex=1.4)
  	ylim.plot=c(0,max(c(ldstat)))
  	if (sum(is.na(ylim.plot))!=0) ylim.plot=c(0,1)
  	matplot(ldstat[-1,],type="l",xlim=c(0,ngen*1.2),ylim=ylim.plot,lty=1,lwd=lwd.val,col=3:4,xlab="",ylab="")
  	legend(ngen,max(ylim.plot),leg=c("|D'|","r^2"),lty=1,col=3:4,lwd=2,cex=1.4)
  }
  return(c(hp[ngen,],ldstat[ngen,]))
}
