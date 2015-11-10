rw2.prec <- function(x)
{
  n = length(x)
  d <- c(Inf,diff(x))
  dm1 = c(d[2:n],Inf)
  dm2 = c(d[3:n],Inf,Inf)
  d1  = c(Inf,d[1:(n-1)])
  d2  = c(Inf,Inf,d[1:(n-2)])
  Q = cBind(2/(dm1^2*(dm2+dm1)) + (2/(dm1*d))*(1/dm1+1/d) + 2/(d^2*(d+d1)),
             -(2/dm1^2)*(1/dm2 + 1/d), 2/(dm2*dm1*(dm2+dm1)))
  return(bandSparse(n=n,m=n,k=c(0,1,2),diagonals=Q,symmetric=1))
}

rw2.A <- function(x,loc)
{
  if(min(loc)< min(x) || max(loc) > max(x))
    stop("locations outside support of basis")

  n.x  <- length(x)
  n.loc <- length(loc)
  i <- as.vector(cBind(1:n.loc,1:n.loc))
  j <- matrix(0,n.loc,2)
  vals <- matrix(1,n.loc,2)
  for(ii in seq_len(n.loc)){
    j[ii,1] <- sum(sum((loc[ii] - x)>=0))
    vals[ii,1] <- loc[ii] - x[j[ii,1]]
    j[ii,2] <- j[ii,1] + 1
    if(j[ii,2]<=n.x){
      vals[ii,2] <- x[j[ii,2]] - loc[ii]
    } else {
      j[ii,2] = j[ii,2] -2
    }
  }
  j <- as.vector(j)
  vals <- as.vector(matrix(1-vals/rowSums(vals)))

  A <- sparseMatrix(i=i,j=j,x=vals, dims=c(n.loc,n.x))
}

rw2.As <- function(x,loc,ind,n)
{
  A <- rw2.A(x,loc)
  v = sparseMatrix(i=1,j=ind,dims=c(1,n))
  return(kronecker(A,matrix(v,1,n)))
}

icar.prec <- function(m,n)
{
  G = - kronecker(toeplitz(sparseVector(1, 2, n)),Diagonal(m,1)) -
        kronecker(Diagonal(n,1),toeplitz(sparseVector(1, 2, m)))
  G = G - Diagonal(m*n,rowSums(G))
return(G)
}
igmrf.prec <- function(m,n)
{
  W1 = sparseMatrix(i=c(2:m),j=c(1:(m-1)),x=rep(1,m-1),
                    dims=c(m, m), symmetric=TRUE)
  W2 = sparseMatrix(i=c(2:n),j=c(1:(n-1)),x=rep(1,n-1),
                    dims=c(n, n), symmetric=TRUE)
  W1[1,2] = 0; W1[m,m-1]=0; W2[1,2] = 0; W2[n,n-1]=0
  W = kronecker(Diagonal(n,1),W1) + kronecker(W2,Diagonal(m,1))

Wc = sparseMatrix(i=c(rep(1,4),rep(m,4),rep(m*(n-1)+1,4),rep(m*n,4)),
                  j=c(1,2,m+1,m+2,
                      m-1,m,2*m-1,2*m,
                      m*(n-2)+1,m*(n-2)+2,m*(n-1)+1,m*(n-1)+2,
                      m*(n-1)-1,m*(n-1),m*n-1,m*n),
                  x=c(1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1),dims=c(n*m,m*n))
  diag(W) <- -rowSums(W)
  W=W + Wc
  return(list(Q=t(W)%*%W,W=W))
}

build.A <- function(Xc,obs.ind,m,n)
{
  T <- dim(Xc)[2]
  D.obs = Diagonal(m*n,1)
  D.obs = D.obs[obs.ind,]
  A <- cBind(Xc[,1]*D.obs,D.obs)
  for(t in 2:T){
    At <- cBind(Xc[,t]*D.obs,D.obs)
    A <- rBind(A,At)
  }
  return(A)
}

llike.mat <- function(p,obj)
{
  tau = exp(p[1])
  kappa = exp(p[2])
  sigma2 = exp(p[3])
  if(tau>1e-16 && kappa>1e-16 && sigma2>1e-16)
  {
    Q = tau*(obj$M0q + kappa*obj$M1q + kappa^2*obj$M2q)
    Rp = update(obj$Rp,Q + obj$AtA/sigma2)
    v = forwardsolve(Rp,obj$AtY/sigma2)
    R = update(obj$R,obj$M0 + kappa*obj$M1 + kappa^2*obj$M2)
    l = log(tau)*(obj$d-2)*obj$n + 2*obj$d*sum(log(diag(R))) -2*sum(log(diag(Rp))) -obj$N*log(sigma2) + t(v)%*%v -obj$yty/sigma2
    return(-l/2)
  } else {
    return(-Inf)
  }
}


llike <- function(p,obj)
{
  if(exp(p[1])>1e-16 && exp(p[2])>1e-16)
  {
    Rp <- chol(exp(p[1])*obj$Q + obj$AtA/exp(p[2]))
    #Rp = update(obj$Rp,exp(p[1])*obj$Q + obj$AtA/exp(p[2]))
    #v = forwardsolve(Rp,obj$AtY/exp(p[2]))
    v <- solve(t(Rp),obj$AtY/exp(p[2]))
    l = p[1]*obj$n -2*sum(log(diag(Rp))) -obj$N*p[2] + t(v)%*%v -obj$yty/exp(p[2])
    return(-as.double(l)/2)
  } else {
    return(-Inf)
  }
}


llike_reverse <- function(p,obj)
{
  if(exp(p[1])>1e-16 && exp(p[2])>1e-16)
  {
    AtQ <- exp(p[1])*t(obj$A)%*%obj$Q
    Qhat <-  obj$I/exp(p[2])  + AtQ%*%obj$A
    Rp <- chol(Qhat)
    v = forwardsolve(Rp,AtQ%*%obj$Y)
    l = p[1]*obj$N -2*sum(log(diag(Rp))) -obj$n*p[2] + t(v)%*%v - exp(p[1])*t(obj$Y)%*%obj$Q%*%obj$Y
    return(-l/2)
  } else {
    return(-Inf)
  }
}

#llike2 defines tau with covariates!
# add covariates to tau
llike2 <- function(p,obj)
{
  #cat(p,'\n')
  d_tau  = dim(obj$Btau)[2]
  Diag_tau = obj$Btau%*%p[1:d_tau]
  if(max(exp(Diag_tau))>1e-14 && exp(p[2])>1e-14)
  {
    Q_tau = .sparseDiagonal(length(Diag_tau),x = as.vector(exp(Diag_tau/2)))
    Q_post = Q_tau%*%obj$Q%*%Q_tau
    Rp = update(obj$Rp,as.spam.dgCMatrix(Q_post) + obj$AtA/exp(p[d_tau+1]), cholupdatesingular="null")
    v = forwardsolve(Rp,obj$AtY/exp(p[d_tau+1]))
    l = sum(Diag_tau) -2*sum(log(diag(Rp))) -obj$N*p[d_tau+1] + t(v)%*%v -obj$yty/exp(p[d_tau+1])
    return(-l/2)
  } else {
    return(-Inf)
  }
}
# model with covariets on:
# tau   - Btau
# sigma - Bsigma
llike3 <- function(p,obj)
{
  d_tau  = dim(obj$Btau)[2]
  Diag_tau = obj$Btau%*%p[1:d_tau]
  d_sigma = dim(obj$Bsigma)[2]
  Diag_sigma = obj$Bsigma%*%exp(p[(d_tau+1):(d_tau + d_sigma)])

  if(min(exp(Diag_tau))>1e-14 && min(exp(Diag_sigma))>1e-14)
  {
    Q_tau = .sparseDiagonal(length(Diag_tau),x = as.vector(exp(Diag_tau/2)))
    Q_post = Q_tau%*%obj$Q%*%Q_tau
    Q_eps = .sparseDiagonal(length(Diag_sigma),x = 1/as.vector((Diag_sigma)))
    AtQ = t(obj$A)%*%Q_eps
    AtQY = AtQ%*%obj$Y
    ytQy = t(obj$Y)%*%Q_eps%*%obj$Y
    AtQA = as.spam.dgCMatrix(as(AtQ%*%obj$A,"dgCMatrix"))
    Rp = update(obj$Rp,as.spam.dgCMatrix(Q_post) + AtQA, cholupdatesingular="null")
    v = forwardsolve(Rp,AtQY)
    l = as.double(sum(Diag_tau) -2*sum(log(diag(Rp))) -sum(log(Diag_sigma)) + t(v)%*%v -ytQy)
    return(-l/2)
  } else {
    return(-Inf)
  }
}

f_quad <- function(X)
{
  F.obs = cBind(X^2,X,rep(1,length(X)))
  return(F.obs)
}
f_lin <- function(X)
{
  F.obs = cBind(X,rep(1,length(X)))
  return(F.obs)
}

#Creating A matrix of regerssion data
# @param FX        - nxm matrix of covariates
# @param obs.index - int location of the indices
# @param N         - int total number of observations
ftoA <- function(FX, obs.index, N)
{
  nm  <- dim(FX)
  i <- rep(1:nm[1], nm[2])
  j <- rep(obs.index, nm[1]*nm[2]) + rep(N * (0:(nm[2]-1)),each=nm[1])
  A <- sparseMatrix(i=i,
                    j=j,
                    x=c(FX),dims=c(nm[1], nm[2]*N))

  return(A)
}

#Creating measurement error from the Delta function in
#Doksum
# assume \hat{F} is emperical distirbution function of X (true is F),
# assume \hat{G} is emperical distribution function of Y (true is G),
# then V[\Delta(x(i)) |x(i)] = V[ \hat{G}^-1(i/m) - x(i)| x(i)]= V[ \hat{G}^-1(i/m)] =
#      V[ Y_1(i/m)] = (i/m)*(1-i/m)/ g^2(G^{-1}(i/m))
# we assume that each square has gamma distribution disregarding the zero observations
# so we then fit for each square Gamma distribution, then get a measurement error for
# each pair (Y_i,X_i) from the equation above (not the above model assumes indepndence
# in variance which obviously is not true but an approximation)
# @param Y         - nxm observations in each square from G
# @param X         - nxm observations in each square from F
QuantileMeasurmentError<- function(Y, X)
{
  if(is.list(X)== 0){
    n <- dim(Y)[1]
    m <- dim(X)[2] + 1
    shape_rate = matrix(0,n,2)
    VX <- matrix(0,dim(X)[1], dim(X)[2])
    for(i in 1:n){
      shape_rate[i,] <- MASS::fitdistr(Y[i,Y[i,]>0],"gamma")$estimate
      index <- sort(X[i,],index.return=T)$ix
      Y_p <- qgamma( index/m,shape = shape_rate[i,1], rate= shape_rate[i,2])
      VX[i, ] <- index/m * (1 - index/m )/dgamma(Y_p, shape = shape_rate[i,1], rate= shape_rate[i,2])^2
    }

    return(list(VX = VX, shape_rate=shape_rate))
  }

  n <- length(Res$x)
  shape_rate = matrix(0,n,2)
  VX <- list()
  for(i in 1:n){
    shape_rate[i,] <- MASS::fitdistr(na.omit(Y[i,Y[i,]>0]),"gamma")$estimate
    index <- sort(X[[i]], index.return=T)$ix
    m <- length(X[[i]]) + 1
    Y_p <- qgamma( index/m,shape = shape_rate[i,1], rate= shape_rate[i,2])
    VX[[i]] <- index/m * (1 - index/m )/dgamma(Y_p, shape = shape_rate[i,1], rate= shape_rate[i,2])^2
  }

  return(list(VX = VX, shape_rate=shape_rate))
}
#Creating measurement error from the Delta function in
#Doksum
# assume \hat{F} is emperical distirbution function of X (true is F),
# assume \hat{G} is emperical distribution function of Y (true is G),
# then V[\Delta(x(i)) |x(i)] = V[ \hat{G}^-1(i/m) - x(i)| x(i)]= V[ \hat{G}^-1(i/m)] =
#      V[ Y_1(i/m)] = (i/m)*(1-i/m)/ g^2(G^{-1}(i/m))
# we assume that each square has gamma distribution disregarding the zero observations
# so we then fit for each square Gamma distribution, then get a measurement error for
# each pair (Y_i,X_i) from the equation above (not the above model assumes indepndence
# in variance which obviously is not true but an approximation)
# @param Y         - nxm observations in each square from G
# @param X         - nxm observations in each square from F
# @param eps       - (double) since the higsht quantile typically has inf we add small bias by reducing by eps
QuantileMeasurmentError2<- function(Y, X,eps = 10^-8)
{

  n <- dim(Y)[1]
  m <- dim(X)[2]
  VX <- matrix(0,dim(X)[1], dim(X)[2])
  Y <- na.omit(c(Y))
  shape_rate <- MASS::fitdistr(Y[Y>0],"gamma")$estimate
  for(i in 1:n){
    index <- sort(X[i,],index.return=T)$ix
    Y_p <- qgamma((index- eps)/m ,shape = shape_rate[1], rate= shape_rate[2])
    VX[i, ] <- (index -  eps)/m * (1 - (index- eps)/m)/dgamma(Y_p, shape = shape_rate[1], rate= shape_rate[2])^2
  }

  return(list(VX = VX, shape_rate=shape_rate))
}

#Creating measurement error from the Delta function in
#Doksum
# assume \hat{F} is emperical distirbution function of X (true is F),
# assume \hat{G} is emperical distribution function of Y (true is G),
# then V[\Delta(x(i)) |x(i)] = V[ \hat{G}^-1(i/m) - x(i)| x(i)]= V[ \hat{G}^-1(i/m)] =
#      V[ Y_1(i/m)] = (i/m)*(1-i/m)/ g^2(G^{-1}(i/m))
# we assume that each square has gamma distribution disregarding the zero observations
# so we then fit for each square Gamma distribution, then fit the full measurment error model
# from the assymtotically distribution of the quantile
# @param Y         - nxm observations in each square from G
# @param X         - nxm observations in each square from F
# @param eps       - (double) since the higsht quantile typically has inf we add small bias by reducing by eps
QuantileMeasurmentError3 <- function(Y, X,eps = 10^-8)
{
  print("WARNING QuantileMeasurmentError3 is very adhoc, and only works for small dimensions and sorted data")
  n <- dim(Y)[1]
  m <- dim(X)[2] + 1 #adding small bias term
  VX <- matrix(0,dim(X)[1], dim(X)[2])
  shape_rate = matrix(0,n,2)
  p <- (1:dim(X)[2])/m
  P <- matrix(0,nrow=length(p),ncol=length(p))
  for(i in 1:length(p))
  {
    for(ii in i:length(p))
    {
      P[i,ii] <- p[i]*(1-p[ii])
      P[ii,i] <- p[i]*(1-p[ii])
    }
  }
  Qe <- solve(P)
  Qe[abs(Qe) <0.01] =0
  Qe <- as(Qe,"dgCMatrix")
  Qe_vec <- list()
  for(i in 1:n){
    shape_rate[i,] <- MASS::fitdistr(Y[i,Y[i,]>0],"gamma")$estimate
    index <- sort(X[i,],index.return=T)$ix
    Y_p <- qgamma(index/m ,shape = shape_rate[i,1], rate= shape_rate[i,2])
    dGamma_vec <- dgamma(Y_p, shape = shape_rate[i,1], rate= shape_rate[i,2])
    VX[i, ] <- index/m * (1 - index/m)/dGamma_vec^2
    Qe_vec[[i]] <-  .sparseDiagonal(dim(X)[2], dGamma_vec)%*%Qe%*%.sparseDiagonal(dim(X)[2], dGamma_vec)
  }

  return(list(VX = VX, shape_rate=shape_rate, Qe_vec=Qe_vec))
}
# creating the data in from of \Delta(x) as observations
# from the duksum article
# @param X             - the data
# @param Y             - the predictive data
# @param remove.unique - makes linear interpolation between the unique points
duksomData <- function(X,Y, minVal = 0, maxVal = Inf, remove.unique=TRUE)
{
  y <- x   <-list()
  for(i in 1:dim(X)[1]){
    x_i = sort(na.omit(X[i,]))
    index = x_i > minVal & x_i < maxVal
    x_i <- x_i[index]
    if(remove.unique==TRUE)
      x_i <- interpolate_duplicate(x_i)
    m = length(x_i)
    index_Y = Y[i,]>minVal & Y[i,] < maxVal & is.na(Y[i,])==F
    Y_i <- sort(Y[i,index_Y])
    if(remove.unique==TRUE)
      Y_i <- interpolate_duplicate(Y_i)
    qy = quantilefun(Y_i)
    #y[[i]] = qy( 0:(m-1) / m ) - x_i alternative form makes no difference!
    y[[i]] = qy((1:m) /(m + 1)) - x_i
    x[[i]] = x_i
  }
  return(list(x = x, y= y) )
}

interpolate_duplicate <- function(x)
{
  index_nondup = which(duplicated(x)==FALSE)
  for(j in 1:(length(index_nondup)-1))
  {

    k <- index_nondup[j+1] - index_nondup[j] + 1
    if(k > 2){
      x_1 <- x[index_nondup[j]]
      x_2 <- x[index_nondup[j+1]]
      x[(index_nondup[j]):(index_nondup[j+1])] <- seq(x_1,x_2,length.out=k)
    }
  }
  return(x)
}

evalDumsom <- function(x, Y, x.grid, f)
{
  Ai <- rw2.A(x.grid, as.matrix(x))
  Delta <- Ai%*%f
  x_ <- x + Delta
  x_sort <- sort(as.vector(x_),index.return=T)
  G <- rep(0,length(x_))
  Y_temp <- na.omit(Y)
  G_p <- 0
  n<- length(Y_temp)
  for( i in 1:length(x_)){
    index = Y_temp < x_sort$x[i]
    G_p <- G_p + sum(index)
    G[x_sort$ix[i]] <- G_p
    Y_temp <- Y_temp[index==FALSE]
  }
  G <- G/n
  return(G)
}
evalDumsom_Empirical<- function(x, Y, Delta.grid,Delta)
{
  Delta_fun <- approxfun(Delta.grid, Delta, method = "linear", rule = 2, f = 0, ties = mean)
  x_ <- x + Delta_fun(x)
  x_sort <- sort(as.vector(x_),index.return=T)
  G <- rep(0,length(x_))
  Y_temp <- na.omit(Y)
  G_p <- 0
  n<- length(Y_temp)
  for( i in 1:length(x_)){
    index = Y_temp < x_sort$x[i]
    G_p <- G_p + sum(index)
    G[x_sort$ix[i]] <- G_p
    Y_temp <- Y_temp[index==FALSE]
  }
  G <- G/n
  return(G)
}
