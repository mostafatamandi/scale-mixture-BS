#W has contaminated(nu1,nu2)


dbs=function(Y,al,be){
  a=(1/al)*(sqrt(Y/be)-sqrt(be/Y))
  dnorm(a)*(Y+be)/(2*al*sqrt(be)*Y^(3/2))
}

dconbs=function(Y,alp,bet,nu1,nu2){ 
  f=(1-nu1)*dbs(Y,alp,bet)+nu1*dbs(Y,alp,bet/nu2)
  return(f)
}

pconbs=function(x, alp, bet, nu1,nu2)
{
  
  f = function(t, alp, bet, nu1,nu2){dconbs(t,alp, bet, nu1,nu2)}
  cdf=integrate(f, lower=0, upper=x, alp=alp, bet=bet, nu1=nu1,nu2=nu2)$value
  return(cdf)
}


#=============================================
CONBS.EM = function(Y, alp, bet, nu1,nu2, tol = 1e-3, max.iter = 500, per = 200,gof=F, n.ks=1e3)
{
  begin=proc.time()[1]
  iter = 1
  n=length(Y)
  b.fn=function(bet, alp, nu, Y) {-sum(log(dconbs(Y,alp,bet,nu1,nu2)))}
  #nu.fn=function(nu, hs1, hs3) {log(nu)+1-digamma(nu)+mean(hs3-hs1) }
  nu1.fn=function(nu1,nu2,bet, alp, Y) {-sum(log(dconbs(Y,alp,bet,nu1,nu2)))}
  nu2.fn=function(nu2,nu1,bet, alp, Y) {-sum(log(dconbs(Y,alp,bet,nu1,nu2)))}
  #h1=function(a,nu,chi0,psi0) {besselK(sqrt(chi0*psi0),a+nu-1/2)}
  #h2=function(a,nu,chi0,psi0) {besselK(sqrt(chi0*psi0),a+nu+1/2)}
  ind.den = dconbs(Y,alp,bet,nu1,nu2)
  lk=LL=logli.old = sum(log(ind.den))
  cat(paste(rep("=",60),sep="",collapse=""),"\n")
  cat('iter =', iter, '\t logli =', logli.old, '\n')
  
  #logt = expression(lgamma((nu+1)/2)-lgamma(nu/2)-log(nu*pi)/2-(nu+1)*log(1+u^2/nu)/2)
  #nu.fn2=function(nu,u) -sum(eval(logt))
  epsilon= Inf
  while((iter <=max.iter) && (epsilon > tol)) {
    
    
    hs1 = (nu2*nu1*dbs(Y,alp,bet/nu2)+(1-nu1)*dbs(Y,alp,bet))/dconbs(Y,alp,bet,nu1,nu2)
    hs2 = (nu2^(-1)*nu1*dbs(Y,alp,bet/nu2)+(1-nu1)*dbs(Y,alp,bet))/dconbs(Y,alp,bet,nu1,nu2)
    #hs3 = C*(bet*(sqrt(chi0/psi0))^(nu-1/2)*a1+bet*(nu-1/2)*(sqrt(chi0/psi0))^(nu-1/2)*log(sqrt(chi0/psi0))*besselK(sqrt(chi0*psi0),nu-1/2)
    # +(T*sqrt(chi0/psi0))^(nu+1/2)*a2+T*(nu+1/2)*(sqrt(chi0/psi0))^(nu+1/2)*log(sqrt(chi0/psi0))*besselK(sqrt(chi0*psi0),nu+1/2))
    # M-step
    
    alp=sqrt(mean((hs1*Y/bet)+(hs2*bet/Y)-2))
    bet=nlminb(start=bet, b.fn, lower=1e-5, upper=200, alp=alp,nu=nu,Y=Y)$par
    #nu=uniroot(nu.fn,hs1=hs1, hs3=hs3, lower=1, upper=100)$root
    nu1=nlminb(start=nu1, nu1.fn, lower=1e-5, upper=1, alp=alp,bet=bet,nu2=nu2,Y=Y)$par
    nu2=nlminb(start=nu2, nu2.fn, lower=1e-5, upper=1, alp=alp,bet=bet,nu1=nu1,Y=Y)$par
    #nu=(-mean(hs3))^(-1)
    ind.den = dconbs(Y,alp,bet,nu1,nu2)
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
  para = c(alp=alp,bet=bet, nu1=nu1,nu2=nu2)
  m = 4
  aic = -2 * logli.new + 2 * m
  bic = -2 * logli.new + log(n) * m
  mds = c(logli = logli.new, m = m, aic = aic, bic = bic)
  end = proc.time()[1]
  if (gof ==T){
    f = function(t, alp, bet, nu1,nu2){dconbs(t,alp, bet, nu1,nu2)}
    cdf=c()
    for (i in 1:n){cdf[i]=integrate(f, lower=0, upper=Y[i], alp=alp, bet=bet, nu1=nu1,nu2=nu2)$value}
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


