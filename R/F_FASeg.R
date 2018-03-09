F_FASeg <- function(Y,uniKmax,multiKmax,qmax=M-1,selection=TRUE,WithoutCorr=FALSE){
  
n = max(Y$position) 
M = max(Y$series) 
N=n*M
 

#The conditions to be fulfilled
  if (qmax>=M){cat("qmax", qmax," needs to be lower or equal to M-1", M-1)}
  if (uniKmax >n) {cat("The maximal number of segments per series uniKmax", uniKmax," needs to be  lower than the length of the series n" ,n,"\n")} 
  if (multiKmax<M) {cat("The total number of segments multiKmax", multiKmax, "needs to be greater than the number of series M", M, "\n")}
  if (N/M !=n){print("The series needs to have the same length n")}
  if (M==1){print("A simple segmentation can be used")}
  if (length(which(is.na(Y$signal)))>=1) {print("Missing data is forbidden")}

#The used functions
######################
FASeg_kq<- function(Y, K, q, uniKmax,n,M){
  N=n*M
  tol = 1e-6  
  maxIter = 1000
  lv    = c()  
  #Segmentation initialization
  para.seg      = c()
  para.seg      = MultiSeg(Y,n,M,uniKmax, K)
  T.mu          = para.seg$mu
  mu.fit = data.frame(  series      = Y$series,
                        position    = Y$position,
                        mean        = rep(T.mu$mean, T.mu$end-T.mu$begin+1) )
  res    = data.frame(  series      = Y$series,
                        position    = Y$position,
                        resid       = (Y$signal-mu.fit$mean))
  Ytilde = matrix(res$resid, n, M)
  E      = Ytilde  
  # Covariance initialization
  Iq    = diag(rep(1, q))
  Einit =t(t(Ytilde) - rowMeans((t(Ytilde))))
  Ytilde = matrix(res$resid, n, M)
  S     = t(Einit)%*%Einit/n
  EigS  = eigen(S)
  s2    = 0
  s2    = mean(EigS$values[(q+1):M])
  EigS$values = EigS$values - s2
  B     = EigS$vectors[, 1:q] %*% diag(sqrt(EigS$values[1:q]), nrow=q, ncol=q)
  Psi   = s2*diag(rep(1, M))
  Diff  = 2*tol
  Iter  = 0  
  #EM
  while ((Diff  > tol) & (Iter < maxIter))
  {
    Iter = Iter +1
    
    # E step
    invPsi = diag(1/diag(Psi))
    W = (Iq + t(B)%*%invPsi%*%B)
    W = solve(W)
    Z = Ytilde%*%invPsi%*%B%*%W
    ZpZ = n*W + t(Z)%*%Z
    
    # M step: Factor Analysis
    B.tmp = (t(Ytilde) %*% Z) %*% solve(ZpZ)
    E = Ytilde-Z%*%t(B.tmp)
    sigma2 = (sum(diag(E%*%t(E)))/n+ sum(diag((B.tmp%*%W%*%t(B.tmp))))) / M
    Psi.tmp = diag(rep(sigma2, M))
    
    # M step: Segmentation
    Y_ZB = data.frame(  series      = Y$series,
                        position    = Y$position,
                        signal      = Y$signal - matrix(Z%*%t(B.tmp), n*M, 1) )
    para.seg     = MultiSeg(Y_ZB,n,M,uniKmax,K)
    T.mu   = para.seg$mu
    mu.fit.tmp = data.frame(series      = Y$series,
                            position    = Y$position,
                            mean        = rep(T.mu$mean, T.mu$end-T.mu$begin+1) )
    res = data.frame( series      = Y$series,
                      position    = Y$position,
                      resid       = (Y$signal-mu.fit.tmp$mean))
    Ytilde = matrix(res$resid, n, M)
    E = Ytilde-Z%*%t(B.tmp)
    # Test
    if (Iter == 2)
    {
      t2 = c(as.vector(B), diag(Psi), mu.fit$mean)
    }
    if (Iter == 3)
    {
      t1 = c(as.vector(B), diag(Psi), mu.fit$mean)
      t0 = c(as.vector(B.tmp), diag(Psi.tmp), mu.fit.tmp$mean)
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
    }
    if (Iter > 3)
    {
      t2 = t1
      t1 = t0
      t0 = c(as.vector(B.tmp), diag(Psi.tmp), mu.fit.tmp$mean)
      tp1 = tp0
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
      Diff = sum((tp0-tp1)^2)
    }    
    # Update
    Psi = Psi.tmp
    B = B.tmp
    mu.fit = mu.fit.tmp
  }  
  S=Psi + B%*%t(B)
  lv=-0.5*N*log(2*pi)-0.5*n*log(det(S))-0.5*sum(diag(Ytilde%*%solve(S)%*%t(Ytilde)))
  invisible(list(B=B, Psi=Psi, Z=Z, T.mu = T.mu, lv=lv))
}
######################
MultiSeg <- function(Y,n,M,uniKmax, multiKmax){
  N=n*M  
  breakpoints  = list()
  contrast     = data.frame()
  m = 0  
  for (ell in (1:M)){
    m                = m+1
    x                = Y[as.factor(Y$series)==ell,]$signal
    matD             = Gsegmentation(x)
    out              = DynProg(matD,uniKmax)
    th               = out[[2]]
    contrast         = rbind(contrast,cbind(rep(ell,uniKmax),c(1:uniKmax), out[[1]]))
    breakpoints[[m]] = th
    rm(matD,out,th,x)
  }   
  colnames(contrast) = c("series", "K","J.est")  
  out.ibp  = DynProg.ibp(contrast,M,multiKmax)
  J.ibp    = out.ibp[[1]]
  seg.rep  = out.ibp[[2]]
  rm(out.ibp)
  Kh.ibp=multiKmax  
  mu         = data.frame()
  m          = 0
  bp.global  = c()
  
  for (ell in (1:M)){
    m         = m+1
    k         = seg.rep[Kh.ibp,m]
    bp.global = c(bp.global,breakpoints[[m]][k,1:k]+(m-1)*n) 
    rupt      = matrix(Inf,ncol = 2 , nrow= k)
    bp        = breakpoints[[m]][k,1:k]    
    rupt[,2]  = bp
    bp        = bp +1
    rupt[,1]  = c(1, bp[1:k-1])
    mu  = rbind(mu,cbind(rep(ell,k), rupt))
    rm(k,rupt,bp)
  } 
  colnames(mu) = c("series", "begin","end")
  rownames(mu) = c(1:dim(mu)[1])
  
  rupt.global     = matrix(Inf,ncol = 2 , nrow= Kh.ibp)
  rupt.global[,2] = bp.global
  bp.global       = bp.global +1
  rupt.global[,1] = c(1, bp.global[1:(Kh.ibp-1)])  
  mu$mean   = apply(rupt.global,1,FUN=function(y) mean(Y$signal[y[1]:y[2]], na.rm=T))
  rm(breakpoints,seg.rep,ell)
  
  list(J.ibp=J.ibp,mu=mu)
  
}
######################
Gsegmentation<-function(x){
  n = length(x)
  matD=matrix(Inf,n,n)
  x2=x^2
  x2i=cumsum(x2)
  xi=cumsum(x)
  x2i=x2i[1:n]
  xi=xi[1:n]
  
  matD[1,1:n]=x2i-((xi^2)/(1:n))
  nl=n
  for (i in 2:nl){
    ni=n-i+2
    x2i=x2i[2:ni]-x2[i-1] 
    xi=xi[2:ni]-x[i-1]
    deno<-(i:n)-i+1
    matD[i,(i:n)]=x2i-((xi^2)/deno)
  }
  invisible(matD)
}
######################
DynProg<-function(matD,Kmax)
{
  n<-dim(matD)[1]
  I<-matrix(Inf,Kmax,n)
  t<-matrix(0,Kmax,n)   
  I[1,]=matD[1,]
  matD=t(matD) 
  if (Kmax>2){
    for (k in 2:(Kmax-1)){
      for (L in k:n){
        I[k,L]<-min(I[(k-1),1:(L-1)]+matD[L,2:L])
        if(I[k,L]!=Inf)
          t[k-1,L]<-which.min(I[(k-1),1:L-1]+matD[L,2:L])
      }
    }
  }
  I[Kmax,n]<-min(I[Kmax-1,1:(n-1)]+matD[n,2:n])
  if(I[Kmax,n]!=Inf)
    t[Kmax-1,n]<-which.min(I[(Kmax-1),1:n-1]+matD[n,2:n]) 
  t.est<-matrix(0,Kmax,Kmax)
  diag(t.est)<-n
  for (K in 2:Kmax){
    for (k in seq(K-1,1,by=-1)){
      if(t.est[K,k+1]!=0)
        t.est[K,k]<-t[k,t.est[K,k+1]]
    }
  }
  list(J.est = I[,n],t.est = t.est)
}
######################
DynProg.ibp<-function(contrast,M,multiKmax){  
  Km        = max(contrast$K)
  contrast  = matrix(contrast$J.est, ncol = M, nrow = Km)
  contrast  = rbind(contrast, matrix(Inf,ncol=M,nrow=(M-1)*Km))
  I         = matrix(Inf,M, multiKmax)
  I[1,1:multiKmax] = contrast[1:multiKmax,1]
  t         = matrix(0,M,multiKmax)  
  for (m in (2:M)){
    for (k in (m:multiKmax)){
      I[m, k] = min( I[(m-1), (1:(k-1))] + contrast[(k-1):1,m])
      if(I[m,k]!=Inf){
        t[m,k]<-which.min( I[(m-1), (1:(k-1))] + contrast[(k-1):1,m])
      }
    }
  }  
  t.est  = matrix(NA,M,multiKmax)
  for (k in (multiKmax:M)){
    tmp1 = t[M,k]
    t.est[M,k] = k - tmp1
    for (m in ((M-1):1)){
      tmp2 = t[m,tmp1]
      t.est[m,k] = tmp1-tmp2
      tmp1=tmp2
    }
  }  
  list(J.est = I[M,],seg.rep = t(t.est))  
} 
######################
######################


if (selection==FALSE){
  SelectedK = multiKmax
  Selectedq = qmax
  if (Selectedq!=0){
  FA =c()
  FA = FASeg_kq(Y, SelectedK, Selectedq, uniKmax,n,M)
  SelectedSigma=FA$B%*%t(FA$B)  + FA$Psi
  SelectedPsi=FA$Psi
  SelectedB=FA$B
  SelectedZ=FA$Z
  SelectedSeg=FA$T.mu
  } else {
    FA =c()
    FA = MultiSeg(Y,n,M,uniKmax, SelectedK)
    SelectedSigma=NULL
    SelectedPsi=NULL
    SelectedB=NULL
    SelectedZ=NULL
    SelectedSeg=FA$mu
  }
} else{

  Kseq=seq(M,multiKmax,1)
  qseq=seq(1,qmax,1)
  BICq = matrix(-Inf, max(Kseq), max(qseq));  BICK = BICq;  SSall = BICq; SSwg = BICq; SSbg = BICq; logLg = BICq;Lv=BICq;      
  Ytmp = matrix(Y$signal, n, M)
  Ybar = matrix(mean(Ytmp), n, M)
  D =  qseq*(2*M-qseq+1)/2 +1
  
  
  for (q in qseq){
    for (Km in Kseq){
        FA=c()
        FA = FASeg_kq(Y, Km, q, uniKmax,n,M)
        Sigma = FA$B%*%t(FA$B)  + FA$Psi
        invSigma = solve(Sigma) 
        Lg = FA$T.mu$end-FA$T.mu$begin + 1
        logLg[Km, q] = sum(log(Lg))
        Lv[Km, q]= FA$lv
        Tmutmp = matrix(0, n, M)
        for (i in (1:length(FA$T.mu$mean))){
          Tmutmp[FA$T.mu$begin[i]:FA$T.mu$end[i],FA$T.mu$series[i]] = FA$T.mu$mean[i]
        }
        SSbg[Km, q]  = sum(diag((Tmutmp - Ybar)%*%invSigma%*%t(Tmutmp - Ybar)))
        SSwg[Km, q]  = sum(diag((Ytmp - Tmutmp)%*%invSigma%*%t(Ytmp - Tmutmp)))
        SSall[Km, q] = sum(diag((Ytmp - Ybar)%*%invSigma%*%t(Ytmp - Ybar)))
        BICK[Km, q] = ((Km-M)/2)*log(SSall[Km, q]/2) + ((N-Km)/2+1)*log(1+(SSbg[Km, q]/SSwg[Km,q]))+ lgamma((N-Km)/2+1) - 0.5*logLg[Km,q] -(Km-M)*log(N)
      }
    BICq[Kseq,q] = 2*Lv[Kseq,q] - D[q]*log(n)
   

  }

qest = apply(BICq, 1, which.max)
BICKest = matrix(-Inf,max(Kseq),1)
for (Km in Kseq){BICKest[Km] = BICK[Km,qest[Km]]}

if (WithoutCorr==FALSE){  
  SelectedK = which.max(BICKest)
  Selectedq = qest[SelectedK]
  FA =c()
  FA = FASeg_kq(Y, SelectedK, Selectedq, uniKmax,n,M)
  SelectedSigma=FA$B%*%t(FA$B)  + FA$Psi
  SelectedPsi=FA$Psi
  SelectedB=FA$B
  SelectedZ=FA$Z
  SelectedSeg=FA$T.mu
} else{
  BICK_only = matrix(-Inf, max(Kseq),1); logLg_only=BICK_only; SSall_only = BICK_only; SSwg_only = BICK_only; SSbg_only = BICK_only;  logLg_only = BICK_only;    
  for (Km in Kseq){
    FA = c()
    FA = MultiSeg(Y,n,M,uniKmax, Km)
    FA$T.mu=FA$mu
    Lg = FA$T.mu$end-FA$T.mu$begin + 1
    logLg_only[Km] = sum(log(Lg))
    Tmutmp = matrix(0, n, M)
    for (i in (1:length(FA$T.mu$mean)))
    {
      Tmutmp[FA$T.mu$begin[i]:FA$T.mu$end[i],FA$T.mu$series[i]] = FA$T.mu$mean[i]
    }

    SSbg_only[Km]  = sum(diag((Tmutmp - Ybar)%*%t(Tmutmp - Ybar)))
    SSwg_only[Km]  = sum(diag((Ytmp - Tmutmp)%*%t(Ytmp - Tmutmp)))
    si2=c()
    si2   = SSwg_only[Km]/N
    SSall_only[Km] = (sum(diag((Ytmp - Ybar)%*%t(Ytmp - Ybar))))/si2
    BICK_only[Km]=((Km-M)/2)*log(SSall_only[Km]/2) + ((N-Km)/2+1)*log(1+(SSbg_only[Km]/SSwg_only[Km])) + lgamma((N-Km)/2+1) - 0.5*logLg_only[Km] -(Km-M)*log(N)
  }

  Kest_q0=which.max(BICK_only)
  Kest= which.max(BICKest)
  rg=which.max(c(BICK_only[Kest_q0],BICKest[Kest]))

  if (rg==1) {
    SelectedK=Kest_q0
    Selectedq=0
    FA = c()
    FA = MultiSeg(Y,n,M,uniKmax, SelectedK)
    SelectedSigma=NULL
    SelectedPsi=NULL
    SelectedB=NULL
    SelectedZ=NULL
    SelectedSeg=FA$mu
    } else{
      SelectedK=Kest
      Selectedq=qest[SelectedK]
      FA =c()
      FA = FASeg_kq(Y, SelectedK, Selectedq, uniKmax,n,M)
      SelectedSigma=FA$B%*%t(FA$B)+ FA$Psi
      SelectedPsi=FA$Psi
      SelectedB=FA$B
      SelectedZ=FA$Z
      SelectedSeg=FA$T.mu
    }
}
}


  
list(SelectedK=SelectedK,Selectedq=Selectedq,SelectedSeg=SelectedSeg,SelectedSigma=SelectedSigma,SelectedPsi=SelectedPsi,SelectedB=SelectedB,SelectedZ=SelectedZ)
 }