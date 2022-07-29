dgambs=function(Y,alp,bet,nu){
  psi0=Y/(bet*alp^2)+2*nu
  chi0=bet/(Y*alp^2)
  C=(nu^nu)/(2*gamma(nu)*besselK(alp^(-2),1/2)*sqrt(bet*Y^3))
  f=C*(Y*besselK(sqrt(chi0*psi0),nu+1/2)/(sqrt(psi0/chi0))^(nu+1/2)+
         bet*besselK(sqrt(chi0*psi0),nu-1/2)/(sqrt(psi0/chi0))^(nu-1/2))
  return(f)
}

pgambs=function(x, alp, bet, nu)
{
  
  f = function(t, alp, bet, nu){dgambs(t,alp, bet, nu)}
  cdf=integrate(f, lower=0, upper=x, alp=alp, bet=bet, nu=nu)$value
  return(cdf)
}

GAMBS.EM = function(Y, alp, bet, nu, tol = 1e-3, max.iter = 500, per = 200,gof=F, n.ks=1e3)
{
  begin=proc.time()[1]
  iter = 1
  n=length(Y)
  b.fn=function(bet, alp, nu, Y) {-sum(log(dgambs(Y,alp,bet,nu)))}
  nu.fn=function(nu, hs1, hs3) {log(nu)+1-digamma(nu)+mean(hs3-hs1) }
  nu2.fn=function(nu,bet, alp, Y) {-sum(log(dgambs(Y,alp,bet,nu)))}
  h1=function(a,nu,chi0,psi0) {besselK(sqrt(chi0*psi0),a+nu-1/2)}
  h2=function(a,nu,chi0,psi0) {besselK(sqrt(chi0*psi0),a+nu+1/2)}
  ind.den = dgambs(Y,alp,bet,nu)
  lk=LL=logli.old = sum(log(ind.den))
  cat(paste(rep("=",60),sep="",collapse=""),"\n")
  cat('iter =', iter, '\t logli =', logli.old, '\n')
  
  #logt = expression(lgamma((nu+1)/2)-lgamma(nu/2)-log(nu*pi)/2-(nu+1)*log(1+u^2/nu)/2)
  #nu.fn2=function(nu,u) -sum(eval(logt))
  epsilon= Inf
  while((iter <=max.iter) && (epsilon > tol)) {
    psi0=Y/(bet*alp^2)+2*nu
    chi0=bet/(Y*alp^2)
    a1=a2=rep(0,n)
    for (i in 1:n){
      a1[i]=grad(h1,0,nu=nu,chi0=chi0[i],psi0=psi0[i])
      a2[i]=grad(h2,0,nu=nu,chi0=chi0[i],psi0=psi0[i])
    }
    
    C=(bet*besselK(sqrt(chi0*psi0),nu-1/2)/((sqrt(psi0/chi0))^(nu-1/2))
       +Y*besselK(sqrt(chi0*psi0),nu+1/2)/((sqrt(psi0/chi0))^(nu+1/2)))^-1
    
    
    
    hs1 = C*(bet*(sqrt(chi0/psi0))^(1+nu-1/2)*besselK(sqrt(chi0*psi0),1+nu-1/2)
             +Y*(sqrt(chi0/psi0))^(1+nu+1/2)*besselK(sqrt(chi0*psi0),1+nu+1/2))
    
    hs2 = C*(bet*(sqrt(chi0/psi0))^(-1+nu-1/2)*besselK(sqrt(chi0*psi0),-1+nu-1/2)
             +Y*(sqrt(chi0/psi0))^(-1+nu+1/2)*besselK(sqrt(chi0*psi0),-1+nu+1/2))
    
    hs3 = C*(bet*(sqrt(chi0/psi0))^(nu-1/2)*a1+bet*(nu-1/2)*(sqrt(chi0/psi0))^(nu-1/2)*log(sqrt(chi0/psi0))*besselK(sqrt(chi0*psi0),nu-1/2)
             +(Y*sqrt(chi0/psi0))^(nu+1/2)*a2+Y*(nu+1/2)*(sqrt(chi0/psi0))^(nu+1/2)*log(sqrt(chi0/psi0))*besselK(sqrt(chi0*psi0),nu+1/2))
    # M-step
    
    alp=sqrt(mean((hs1*Y/bet)+(hs2*bet/Y)-2))
    bet=nlminb(start=bet, b.fn, lower=1e-5, upper=200, alp=alp,nu=nu,Y=Y)$par
    #nu=uniroot(nu.fn,hs1=hs1, hs3=hs3, lower=1, upper=100)$root
    nu=nlminb(start=nu, nu2.fn, lower=1, upper=100, alp=alp,bet=bet,Y=Y)$par
    
    ind.den = dgambs(Y,alp,bet,nu)
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
    f = function(t, alp, bet, nu){dgambs(t,alp, bet, nu)}
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


#=========================================================================
BS.est=function(Y,alp,bet,gof=F,n.ks=1000){
  begin = proc.time()[1]
  n=length(Y)
  fitbs=nlminb(c(alp,bet),logbs,lower=c(0,0),upper=c(100,100),y=Y)
  al0=fitbs$par[1];be0=fitbs$par[2]
  logl.bs=-fitbs$objective
  aic.bs = 2 * fitbs$objective + 2 * 2
  bic.bs = 2 * fitbs$objective + log(n) * 2
  
  mds = c(logli = logl.bs, m = 2, aic = aic.bs, bic = bic.bs)
  para=c(al0,be0)
  end = proc.time()[1]
  if (gof ==T){
    f = function(t, alp, bet){dbs(t,alp, bet)}
    cdf=c()
    for (i in 1:n){cdf[i]=integrate(f, lower=0, upper=Y[i], alp=al0, bet=be0)$value}
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
  #cat("BS.EM takes", end - begin, "seconds\n\n")
  return(list(CT=end - begin, mds = mds, para=para))
  
}
