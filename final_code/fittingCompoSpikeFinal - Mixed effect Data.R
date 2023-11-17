fitting <- function(y, X=Delta, BS1, cliques=clique.list, id=id, Als=Als, Total_itr=1000){
  K <- length(cliques)
  #n <- nrow(X)
  p <- ncol(X)
  
  #+diag(p)/100
  out <- nnls::nnls(crossprod(X), crossprod(X, y))
  b <- out$x#rep(0, p)
  
  nu <- 1
  lambda <- 1
  
  sigmals <- rep(10, p) #runif(p, 0, 10)
  invsigmals <- rep(10, p) #runif(p, 0, 10)
  
  #sig <- 1
  sig <- sd(y-array(X%*%b))
  #gam    <- sapply(1:length(cliques), FUN=function(i) {prod(dnorm(b[cliques[[i]]], sd = sigmals[cliques[[i]]]*sig))})
  #gamvec <- rep(1/length(cliques), length(cliques))#unlist(gam) / sum(unlist(gam))
  
  Amat <- matrix(0, length(cliques), p)
  
  invsigmals <- list()
  
  for(i in 1:nrow(Amat)){
    invsigmals[[i]] <- rep(1, length(cliques[[i]]))
    Amat[i, cliques[[i]]] <- rep(1, length(cliques[[i]]))
  }
  
  sigcls <- 1/sqrt(clique.size)#rep(10, length(cliques))#runif(length(cliques))
  #sigmals <- #unlist(sapply(1:length(cliques), FUN=function(i){rep(sigcls[i], length(cliques[[i]]))}))
  
  meancls <- rep(0, length(cliques))
  
  Agamvec <- colSums(Amat)
  dvec    <- Agamvec #/ sigmals^2
  
  theta <- rep(0.5, length(cliques))
  a = c = 1
  v0 = 0.1
  v1 = 1
  
  # Ei <- function(i){
  #   ni <- which(id==i)
  #   Gamma <- solve(crossprod(Als[[i]])/sig^2 + Omega^{-1})
  #   mu    <- Gamma%*%crossprod(Als[[i]], y[ni] - X[ni, ]%*%b)/sigma^2
  #   
  #   return(mu)
  # }
  # 
  # E2i <- function(i){
  #   ni <- which(id==i)
  #   Gamma <- solve(crossprod(Als[[i]])/sig^2 + Omega^{-1})
  #   mu    <- Gamma%*%crossprod(Als[[i]], y[ni] - X[ni, ]%*%b)/sigma^2
  #   
  #   return(tcrossprod(mu) + Gamma)
  # }
  
  theta <- rep(0.5, length(cliques))
  p1 <- length(unlist(cliques))
  ind <- 1:p
  Omega <- diag(J)
  thetaF <- a
  gamma <- rep(1/length(clique.list), length(clique.list))
  #BS1 <- bsplineS((1:7)/7, breaks = seq(0, 1, 1/knot))
  z <- as.vector(y - array(X%*%b))
  for(itr in 1:Total_itr){
    ##Sample random variables
    #eta <- unlist(lapply(1:n, Ei))
    
    #b <- as.vector(solve(t(X)%*%X)%*%t(X)%*%(y-u))
    z <- as.vector(y - array(X%*%b))
    
    #Mixed effect other parts
    iO <- solve(Omega)
    Ts <- R <- C <- 0
    mu <- u <- NULL
    for (i in 1:n) {
      row.i <- which(id==i)
      Xi <- X[row.i,]
      Ai <- Als[[i]]
      AAi <- t(Ai)%*%Ai
      zi <- z[row.i]
      Gammai <- solve(AAi/sig^2 + iO)
      mui <- (Gammai%*%t(Ai)%*%zi)/sig^2
      mu <- c(mu, mui)
      u <- c(u, Ai%*%mui)
      Si <- Gammai + mui%*%t(mui)
      R <- R + Si
      Ts <- Ts + sum(diag(Si%*%AAi))
      C <- C + t(mui)%*%t(Ai)%*%zi
    }
    #sigma2 <- (sum(z^2) -2*C[1] + Ts)/n
    Omega <- as.matrix(R/n)
    
    yred <- y - u
    
    num <- sum((yred - array(X%*%b))^2) + nu*lambda #+ sum(b^2*dvec)/K#sum(unlist(ret)) #sum(b^2*dvec)
    dem <- length(yred)+nu +1 #+sum(colSums(Amat))
    
    sig = sqrt(num/dem)
    
    #solve(crossprod(X) + diag(dvec))#
    #if(itr>1){
    #ind <- which(dvec>0.01) 
    #}
    #M <- solve(crossprod(X[,ind]) + diag(dvec[ind]))#diag(1/dvec) - diag(1/dvec) %*% t(X) %*% solve(diag(n)+ X%*%diag(1/dvec)%*%t(X)) %*% X %*% diag(1/dvec)
    out <- nnls::nnls(crossprod(X[,ind]) + sig^2*diag(dvec[ind]), crossprod(X[,ind], yred))#+colSums(Amat * gamvec * meancls/sigcls^2))#M %*% crossprod(X[,ind], y)
    b[ind] <- out$x
    #b[-ind] <- 0
    
    #bconti   <- unlist(sapply(1:length(cliques), FUN=function(i) {sum((b[cliques[[i]]])^2)/sigcls[i]^2}))
    
    #ret <- sapply(1:length(cliques), FUN=function(i) {gamvec[i]*sum((b[cliques[[i]]]-meancls[i])^2*invsigmals[cliques[[i]]])})
    
    
    Amat <- matrix(0, length(cliques), p)
    
    for(i in 1:nrow(Amat)){
      Amat[i, cliques[[i]]] <- invsigmals[[i]]
    }
    
    
    Agamvec <- colSums(Amat)
    dvec    <- Agamvec
    #thetaF <- a
    s <- rep(0, length(theta))
    for(i in 1:length(cliques)){
      bi <- b[cliques[[i]]]#-meancls[i]
      p0 <- 1-gamma[i] + (gamma[i])*(1-theta[i])
      p1 <- 1-p0
      temp <- dnorm(bi, sd = sqrt(v1), log=T)-dnorm(bi, sd = sqrt(v0), log=T)
      pstar   <- p1/(p1+p0*exp(-temp))
      if(any(temp>100)) {pstar[which(temp>100)] <- 1}
      if(any(temp< -100)) {pstar[which(temp< -100)] <- 0}
      invsigmals[[i]] <- ((1-pstar)/v0 + pstar/v1)
      s[i] <- sum(pstar)
      #theta[i] <- (sum(pstar) + thetaF-1)/(thetaF+c+length(bi)-2)
    }
    gamma <- s/sum(s)
    theta <- s/gamma
    #thetaF <- 1/mean(log(1/theta), na.rm=T)
    #theta <- rdirichlet(1, s+1)
    #plot(BS1 %*% Delta %*% b)
    #print(which.max(theta))
    #bvec   <- unlist(sapply(1:length(cliques), FUN=function(i) {mean((b[cliques[[i]]])^2)}))#, sd = sigmals[cliques[[i]]]*sig))})
    #plot(bvec)
    #plot(bvec)
    #plot(theta)
    #print(mean((ftf - BS1 %*% Delta %*% b)^2))
  }
  #ind <- which(gamma>1e-2)
  #mean((ftf - BS1 %*% Delta[, unlist(clique.list[ind])] %*% b[unlist(clique.list[ind])])^2)
  return(list(gamma=gamma,sigma=sig,esti=BS1 %*% Delta %*% b))
  #return(list(betaest = b, thetals = theta))
}