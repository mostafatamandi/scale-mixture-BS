#W has beta(nu,1)
library(ghyp)

pgi=function(x,lambda,chi0,psi0){
  #ab=sqrt(del1*del2);
  #k=sqrt(del1/del2)*besselK(ab,lambda+1)/besselK(ab,lambda);
  return(pgig(x,lambda=lambda,chi=chi0,psi=psi0))
  #return(pghyp(x,object=ghyp(mu=0,sigma=k,lambda=lambda,alpha.bar=ab,gamma=k*b)))
}
##Vectorizing pgig based on observes and chi and psi;==========
PGIG=Vectorize(pgi,vectorize.args=c("x","chi0","psi0"))

dbetbs=function(Y,alp,bet,nu){ 
  psi0=Y/(bet*alp^2)
  chi0=bet/(Y*alp^2)
  C=nu/(2*besselK(alp^(-2),1/2)*sqrt(bet*Y^3))
  f=C*(Y*(besselK(sqrt(chi0*psi0),nu+1/2)/(sqrt(psi0/chi0))^(nu+1/2))*PGIG(1,lambda=nu+1/2,chi0,psi0)+
         bet*(besselK(sqrt(chi0*psi0),nu-1/2)/(sqrt(psi0/chi0))^(nu-1/2))*PGIG(1,lambda=nu-1/2,chi0,psi0))
  return(f)
}

pbetbs=function(x, alp, bet, nu)
{
  
  f = function(t, alp, bet, nu){dbetbs(t,alp, bet, nu)}
  cdf=integrate(f, lower=0, upper=x, alp=alp, bet=bet, nu=nu)$value
  return(cdf)
}
#=============================================
BETBS.EM = function(Y, alp, bet, nu, tol = 1e-4, max.iter = 500, per = 200,gof=F, n.ks=1e3)
{
  begin=proc.time()[1]
  iter = 1
  n=length(Y)
  b.fn=function(bet, alp, nu, Y) {-sum(log(dbetbs(Y,alp,bet,nu)))}
  #nu.fn=function(nu, hs1, hs3) {log(nu)+1-digamma(nu)+mean(hs3-hs1) }
  nu2.fn=function(nu,bet, alp, Y) {-sum(log(dbetbs(Y,alp,bet,nu)))}
  #h1=function(a,nu,chi0,psi0) {besselK(sqrt(chi0*psi0),a+nu-1/2)}
  #h2=function(a,nu,chi0,psi0) {besselK(sqrt(chi0*psi0),a+nu+1/2)}
  ind.den = dbetbs(Y,alp,bet,nu)
  lk=LL=logli.old = sum(log(ind.den))
  cat(paste(rep("=",60),sep="",collapse=""),"\n")
  cat('iter =', iter, '\t logli =', logli.old, '\n')
  
  #logt = expression(lgamma((nu+1)/2)-lgamma(nu/2)-log(nu*pi)/2-(nu+1)*log(1+u^2/nu)/2)
  #nu.fn2=function(nu,u) -sum(eval(logt))
  epsilon= Inf
  while((iter <=max.iter) && (epsilon > tol)) {
    psi0=Y/(bet*alp^2)
    chi0=bet/(Y*alp^2)
    #a1=a2=rep(0,n)
    #for (i in 1:n){
      #a1[i]=grad(h1,0,nu=nu,chi0=chi0[i],psi0=psi0[i])
      #a2[i]=grad(h2,0,nu=nu,chi0=chi0[i],psi0=psi0[i])
    #}
    
    C.inv=(bet*besselK(sqrt(chi0*psi0),nu-1/2)/((sqrt(psi0/chi0))^(nu-1/2))*PGIG(1,lambda=nu-1/2,chi0,psi0)
       +Y*besselK(sqrt(chi0*psi0),nu+1/2)/((sqrt(psi0/chi0))^(nu+1/2))*PGIG(1,lambda=nu+1/2,chi0,psi0))
    
    C=1/C.inv
    
    hs1 = C*(bet*besselK(sqrt(chi0*psi0),nu+1/2)/((sqrt(psi0/chi0))^(nu+1/2))*PGIG(1,lambda=nu+1/2,chi0,psi0)
             +Y*besselK(sqrt(chi0*psi0),nu+3/2)/((sqrt(psi0/chi0))^(nu+3/2))*PGIG(1,lambda=nu+3/2,chi0,psi0))
    
    hs2 = C*(bet*besselK(sqrt(chi0*psi0),nu-3/2)/((sqrt(psi0/chi0))^(nu-3/2))*PGIG(1,lambda=nu-3/2,chi0,psi0)
             +Y*besselK(sqrt(chi0*psi0),nu-1/2)/((sqrt(psi0/chi0))^(nu-1/2))*PGIG(1,lambda=nu-1/2,chi0,psi0))
    
    #hs3 = C*(bet*(sqrt(chi0/psi0))^(nu-1/2)*a1+bet*(nu-1/2)*(sqrt(chi0/psi0))^(nu-1/2)*log(sqrt(chi0/psi0))*besselK(sqrt(chi0*psi0),nu-1/2)
            # +(Y*sqrt(chi0/psi0))^(nu+1/2)*a2+Y*(nu+1/2)*(sqrt(chi0/psi0))^(nu+1/2)*log(sqrt(chi0/psi0))*besselK(sqrt(chi0*psi0),nu+1/2))
    # M-step
    
    alp=sqrt(mean((hs1*Y/bet)+(hs2*bet/Y)-2))
    bet=nlminb(start=bet, b.fn, lower=1e-5, upper=200, alp=alp,nu=nu,Y=Y)$par
    #nu=uniroot(nu.fn,hs1=hs1, hs3=hs3, lower=1, upper=100)$root
    nu=nlminb(start=nu, nu2.fn, lower=1e-4, upper=200, alp=alp,bet=bet,Y=Y)$par
    #nu=(-mean(hs3))^(-1)
    ind.den = dbetbs(Y,alp,bet,nu)
    newLL=logli.new = sum(log(ind.den))
    lk = c(lk, newLL)
    # Aikten's method
    if (iter < 2) epsilon = abs(LL-newLL)/abs(newLL)
    else 
    {
      tmp  = (newLL - LL)/(LL-lk[iter-1])
      tmp2 =  LL + (newLL-LL)/(1-tmp)
      epsilon = abs(tmp2-newLL)
    }
    diff = newLL - LL
    if(diff < 0)   cat('iter=', iter, "logL is decreasing!\n")
    LL= newLL
    if(iter %% per == 0 | is.na(diff))
    {
      cat('iter =', iter, "\n")
      cat("newLL =", newLL, '\t diff =', diff, "\t Aikten's diff =", epsilon, "\n")
      cat(paste(rep("-", 60), sep = "", collapse = ""), "\n")
    }
    if(is.na(diff)) break
    iter=iter+1
  }
  cat(paste(rep("=",60),sep="",collapse=""),"\n")
  cat('\niter =', iter, '\t logli =', logli.new, '\t diff =', diff, '\n\n')
  para = c(alp=alp,bet=bet, nu=nu)
  m = 3
  aic = -2 * logli.new + 2 * m
  bic = -2 * logli.new + log(n) * m
  mds = c(logli = logli.new, m = m, aic = aic, bic = bic)
  end = proc.time()[1]
  if (gof ==T){
    f = function(t, alp, bet, nu){dbetbs(t,alp, bet, nu)}
    cdf=c()
    for (i in 1:n){cdf[i]=integrate(f, lower=0, upper=Y[i], alp=alp, bet=bet, nu=nu)$value}
    order.cdf=sort(cdf)
    j=1:n
    d=max(j/n-order.cdf,order.cdf-(j-1)/n)
    g=c()
    for (i in 1:n){ g[i]=((2*i)-1)*log(order.cdf[i]*(1-order.cdf[n+1-i]))}
    A2=-n-(1/n)*sum(g)
    I=c()
    for(i in 1: n.ks){
      u=sort(runif(n))
      I[i]=max(j/n-u, u-(j-1)/n)
    }
    mds=round(c(mds,ks.statistic=d,ks.p.value=mean(I>d),ad.statistic=A2,ad.p.value=mean(I>A2)),3)
  }
  
  
  cat(paste(rep("=", 80), sep = "", collapse = ""), "\n")
  return(list(CT=end - begin, mds = mds, para=para))
  #cat("SMBS.EM takes", end - begin, "seconds\n\n")
}


#logl=function(t,y) -sum(log(dbetbs(y,t[1],t[2],t[3])))
#fopt=nlminb(c(1,1,2),logl,lower=c(1e-4,1e-4,1e-4),upper=c(Inf,Inf,Inf),y=T)



