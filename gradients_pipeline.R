######################################################################
## R script to demonstrate partial derivatives of l and n*          ##
######################################################################
## James Liley
## 1 Oct

######################################################################
## R script to demonstrate partial derivatives of l and n*          ##
######################################################################
##
## James Liley
## October 2021
##
## This script should be ran in the directory in which it is saved, or in
##  some directory with subdirectories 'data', 'figures'.
## Not all figures are necessarily used in the manuscript.
##
## For simplicity of demonstration, this script does not use the R
##  package OptHoldoutSize and is stand-alone.


######################################################################
## Switches                                                         ##
######################################################################

# Set to TRUE to save plot to PDF, FALSE to draw to R graphics device
save_plot=FALSE



######################################################################
## Parameters                                                       ##
######################################################################

Ns=10000 # value of N we will usually use
nr=1:floor(Ns/20) # values of n we will consider
ns=100 # test value of n

# Parameters of learning curve; k2(n,theta=c(a,b,c)) = a n^(-b) + c
as=1.5; bs=1.5; cs=0.25;
# can assume cs<k1s since cs is the minimum cost for a sample attainable with a model in place

# value of k1 we will usually use
k1s=0.8



######################################################################
## Functions; worked in full for demonstration's sake               ##
######################################################################

### Cost function

# cost function; usually called l
l=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (a*n^(-b) + c)*(N-n) + k1*n


### Inverse of cost function (n*; optimal holdout size)

# nstar=n*={n:l'(n;t)=0} where t={a,b,c,k1,N}
nstar=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  return(uniroot(function(n) dldn(n,a,b,c,k1,N),c(1,N-1))$root)
}


### Minimum cost function (l(n*))
mincost=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  l(nstar(a,b,c,k1,N),a,b,c,k1,N)
}


### Partials of cost function w.r.t n

# partial derivative of l with respect to n
dldn=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (-b)*a*N*n^(-b-1) - (-b+1)*a*n^(-b) + k1 - c

# second partial derivative of l with respect to n
dldn2=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (-b)*(-b-1)*a*N*n^(-b-2) - (-b)*(-b+1)*a*n^(-b-1)

# del_ l /  del_ a
dlda=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (N-n)*n^(-b)

# del_ l /  del_ b
dldb=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (-log(n))*(N-n)*a*n^(-b)

# del_ l /  del_ c
dldc=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (N-n)

# del_ l /  del_ k1
dldk1=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  n

# del_ l /  del_ N
dldN=function(n,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (a*n^(-b) + c)







### Partials of cost function w.r.t n, theta

# del^2 l / del_n del_a
dda=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (-b)*N*n^(-b-1) - (-b+1)*n^(-b)

# del^2 l / del_n del_b
ddb=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  a*N*(  (-1)*n^(-b-1) + (-b)*(-log(n))*n^(-b-1))-a*((-1)*n^(-b) + (-b+1)*(-log(n))*n^(-b) )

# del^2 l / del_n del_c
ddc=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  return(-1)

# del^2 l / del_n del_k1
ddk1=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  return(1)

# del^2 l / del_n del_N
ddN=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (-b)*a*n^(-b-1)







### Partials of n*

# partial derivative of n* wrt a
dnstarda=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -dda(n,a,b,c,k1,N)/dldn2(n,a,b,c,k1,N)

# partial derivative of n* wrt a, expanded
dnstarda_spec=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  (b*N*n - (b-1)*(n^2))/(b*(b+1)*a*N - b*(b-1)*a*n)


# partial derivative of n* wrt b
dnstardb=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -ddb(n,a,b,c,k1,N)/dldn2(n,a,b,c,k1,N)

# partial derivative of n* wrt b, expanded
dnstardb_spec=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -(n*N*( -1 + b*log(n))-(n^2)*(- 1 + (b-1)*log(n)))/
  (b*(b+1)*N - b*(b-1)*n)



# partial derivative of n* wrt c
dnstardc=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -ddc(n,a,b,c,k1,N)/dldn2(n,a,b,c,k1,N)

# partial derivative of n* wrt c, expanded
dnstardc_spec=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  1/(b*(b+1)*a*N*n^(-b-2) - b*(b-1)*a*n^(-b-1))



# partial derivative of n* wrt k1
dnstardk1=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -ddk1(n,a,b,c,k1,N)/dldn2(n,a,b,c,k1,N)

# partial derivative of n* wrt k1, expanded
dnstardk1_spec=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -1/(b*(b+1)*a*N*n^(-b-2) - b*(b-1)*a*n^(-b-1))



# partial derivative of n* wrt N
dnstardN=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  -ddN(n,a,b,c,k1,N)/dldn2(n,a,b,c,k1,N)

# partial derivative of n* wrt N, expanded
dnstardN_spec=function(n=ns,a=as,b=bs,c=cs,k1=k1s,N=Ns)
  b*n/(b*(b+1)*N - b*(b-1)*n)





## Partials of minimum cost

# partial derivative of l(n*) wrt a
dmcda=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  nstar=nstar(a,b,c,k1,N)
  dlda(nstar,a,b,c,k1,N)
}

# partial derivative of l(n*) wrt b
dmcdb=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  nstar=nstar(a,b,c,k1,N)
  dldb(nstar,a,b,c,k1,N)
}

# partial derivative of l(n*) wrt c
dmcdc=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  nstar=nstar(a,b,c,k1,N)
  dldc(nstar,a,b,c,k1,N)
}

# partial derivative of l(n*) wrt k1
dmcdk1=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  nstar=nstar(a,b,c,k1,N)
  dldk1(nstar,a,b,c,k1,N)
}

# partial derivative of l(n*) wrt N
dmcdN=function(a=as,b=bs,c=cs,k1=k1s,N=Ns) {
  nstar=nstar(a,b,c,k1,N)
  dldN(nstar,a,b,c,k1,N)
}


######################################################################
## Plot partials of n*                                              ##
######################################################################

# Setup
if (save_plot) pdf(paste0("figures/partials_demo.pdf"),width=7.5,height=5)
par(mfrow=c(2,3),mar=c(4,4,1,1))

## Plot cost function and minimum as reference
plot(nr,l(nr),type="l",ylim=c(0,2*l(max(nr))),xlab="n (h.o. set size)", ylab="l (cost)")
nm=nstar() # Equivalent to optimal_holdout_size(Ns,k1s,c(as,bs,cs))$size
mc=mincost() # Equivalent to optimal_holdout_size(Ns,k1s,c(as,bs,cs))$cost
abline(v=nm,col="red")
abline(h=mc,col="blue")
legend("topright",lty=1,
       c(expression(paste("l(n"["*"],")")),
         expression(paste("l(n"["*"],")")),
         "Cost"),
       col=c("red","blue","black"))


aseq=seq(0.5,5,length=200)
xseq=rep(0,length(aseq)); for (i in 1:length(xseq)) xseq[i]=nstar(aseq[i],bs,cs,k1s,Ns)
plot(aseq,xseq,type="l",ylim=c(0,200),xlab="a",ylab="n*");
abline(v=as)
dma=dnstarda(nm,as,bs,cs,k1s,Ns) ## Equivalent to dma=grad_nstar_powerlaw(Ns,k1s,c(as,bs,cs))[3]
abline(nm - dma*as,dma,col="red")
legend("bottomright",c(expression("n"["*"]),expression(paste(delta,"n"["*"],"/",delta,"a"))),lty=1,col=c("black","red"))


bseq=seq(1.01,5,length=200)
xseq=rep(0,length(bseq)); for (i in 1:length(xseq)) xseq[i]=nstar(as,bseq[i],cs,k1s,Ns)
plot(bseq,xseq,type="l",ylim=c(0,200),xlab="b",ylab="n*");
abline(v=bs)
dmb=dnstardb(nm,as,bs,cs,k1s,Ns)
abline(nm - dmb*bs,dmb,col="red")
legend("topright",c(expression("n"["*"]),expression(paste(delta,"n"["*"],"/",delta,"b"))),lty=1,col=c("black","red"))


cseq=seq(0,0.99*k1s,length=200)
xseq=rep(0,length(cseq)); for (i in 1:length(xseq)) xseq[i]=nstar(as,bs,cseq[i],k1s,Ns)
plot(cseq,xseq,type="l",ylim=c(0,200),xlab="c",ylab="n*");
abline(v=cs)
dmc=dnstardc(nm,as,bs,cs,k1s,Ns)
abline(nm - dmc*cs,dmc,col="red")
legend("bottomright",c(expression("n"["*"]),expression(paste(delta,"n"["*"],"/",delta,"c"))),lty=1,col=c("black","red"))



k1seq=seq(1.01*cs,2,length=200)
xseq=rep(0,length(k1seq)); for (i in 1:length(xseq)) xseq[i]=nstar(as,bs,cs,k1seq[i],Ns)
plot(k1seq,xseq,type="l",ylim=c(0,200),xlab=expression("k"[1]),ylab="n*");
abline(v=k1s)
dmk1=dnstardk1(nm,as,bs,cs,k1s,Ns)
abline(nm - dmk1*k1s,dmk1,col="red")
legend("topright",c(expression("n"["*"]),expression(paste(delta,"n"["*"],"/",delta,"k"[1]))),lty=1,col=c("black","red"))


Nseq=seq(200,15000,length=200)
xseq=rep(0,length(Nseq)); for (i in 1:length(xseq)) xseq[i]=nstar(as,bs,cs,k1s,Nseq[i])
plot(Nseq,xseq,type="l",ylim=c(0,200),xlab="N",ylab="n*");
abline(v=Ns)
dmN=dnstardN(nm,as,bs,cs,k1s,Ns)
abline(nm - dmN*Ns,dmN,col="red")
legend("topleft",c(expression("n"["*"]),expression(paste(delta,"n"["*"],"/",delta,"N"))),lty=1,col=c("black","red"))

if (save_plot) dev.off()




######################################################################
## Plot partials of l(n*)                                           ##
######################################################################

if (save_plot) pdf(paste0("figures/partials_demo_ln.pdf"),width=7.5,height=5)


par(mfrow=c(2,3),mar=c(4,4,1,1))
yr=c(0,2800)

plot(nr,l(nr),type="l",ylim=c(0,2*l(max(nr))),xlab="n (h.o. set size)", ylab="l (cost)")
nm=nstar()
mc=mincost()
abline(v=nm,col="red")
abline(h=mc,col="blue")
legend("topright",lty=1,c(expression(paste("l(n"["*"],")")),expression(paste("l(n"["*"],")")),"Cost"),
       col=c("red","blue","black"))

aseq=seq(0.5,5,length=200)
xseq=rep(0,length(aseq)); for (i in 1:length(xseq)) xseq[i]=mincost(aseq[i],bs,cs,k1s,Ns)
plot(aseq,xseq,type="l",ylim=yr,xlab="a",ylab="n*");
abline(v=as)
dma=dmcda(as,bs,cs,k1s,Ns)
abline(mc - dma*as,dma,col="red")
legend("bottomright",c(expression(paste("l(n"["*"],")")),expression(paste(delta,"l(n"["*"],")/",delta,"a"))),lty=1,col=c("black","red"))


bseq=seq(1.01,5,length=200)
xseq=rep(0,length(bseq)); for (i in 1:length(xseq)) xseq[i]=mincost(as,bseq[i],cs,k1s,Ns)
plot(bseq,xseq,type="l",ylim=yr,xlab="b",ylab="n*");
abline(v=bs)
dmb=dmcdb(as,bs,cs,k1s,Ns)
abline(mc - dmb*bs,dmb,col="red")
legend("bottomright",c(expression(paste("l(n"["*"],")")),expression(paste(delta,"l(n"["*"],")/",delta,"b"))),lty=1,col=c("black","red"))


cseq=seq(0,0.99*k1s,length=200)
xseq=rep(0,length(cseq)); for (i in 1:length(xseq)) xseq[i]=mincost(as,bs,cseq[i],k1s,Ns)
plot(cseq,xseq,type="l",ylim=yr,xlab="c",ylab="n*");
abline(v=cs)
dmc=dmcdc(as,bs,cs,k1s,Ns)
abline(mc - dmc*cs,dmc,col="red")
legend("bottomright",c(expression(paste("l(n"["*"],")")),expression(paste(delta,"l(n"["*"],")/",delta,"c"))),lty=1,col=c("black","red"))



k1seq=seq(1.01*cs,2,length=200)
xseq=rep(0,length(k1seq)); for (i in 1:length(xseq)) xseq[i]=mincost(as,bs,cs,k1seq[i],Ns)
plot(k1seq,xseq,type="l",ylim=yr,xlab=expression("k"[1]),ylab="n*");
abline(v=k1s)
dmk1=dmcdk1(as,bs,cs,k1s,Ns)
abline(mc - dmk1*k1s,dmk1,col="red")
legend("bottomright",c(expression(paste("l(n"["*"],")")),expression(paste(delta,"l(n"["*"],")/",delta,"k"[1]))),lty=1,col=c("black","red"))


Nseq=seq(200,15000,length=200)
xseq=rep(0,length(Nseq)); for (i in 1:length(xseq)) xseq[i]=mincost(as,bs,cs,k1s,Nseq[i])
plot(Nseq,xseq,type="l",ylim=yr,xlab="N",ylab="n*");
abline(v=Ns)
dmN=dmcdN(as,bs,cs,k1s,Ns)
abline(mc - dmN*Ns,dmN,col="red")
legend("topleft",c(expression(paste("l(n"["*"],")")),expression(paste(delta,"l(n"["*"],")/",delta,"N"))),lty=1,col=c("black","red"))

if (save_plot) dev.off()


