fittingUN <- function(y, xt, Ti, splin=T, Total_itr=1000){
  BS1 <- diag(Ti)
  BS <- matrix(0, length(y), Ti)
  BS[cbind(1:length(y), xt)] <- 1#bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
  X1 <- BS
  if(splin){
    X1 <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
    BS1 <- bsplineS((1:Ti)/Ti, breaks = seq(0, 1, 1/knot)) 
  }
  #K <- length(cliques)
  #n <- nrow(X)
  p <- ncol(X1)
  
  b <- solve(crossprod(X1)+diag(p)/100)%*%crossprod(X1, y) #rep(0, p)
  #model <- glmnet::cv.glmnet(X1, y) #rep(0, p)
  #b = as.array(coef(model))
  #b = b[-1]
  nu <- 1
  lambda <- 1
  
  sigmals <- rep(10, p) #runif(p, 0, 10)
  invsigmals <- rep(10, p) #runif(p, 0, 10)
  
  sig <- sd(y-array(X1%*%b))
  
  #gam    <- sapply(1:length(cliques), FUN=function(i) {prod(dnorm(b[cliques[[i]]], sd = sigmals[cliques[[i]]]*sig))})
  #gamvec <- rep(1/length(cliques), length(cliques))#unlist(gam) / sum(unlist(gam))
  
  Amat <- matrix(0, length(cliques), p)
  
  invsigmals <- list()
  
  # for(i in 1:nrow(Amat)){
  #   invsigmals[[i]] <- rep(1, length(cliques[[i]]))
  #   Amat[i, cliques[[i]]] <- rep(1, length(cliques[[i]]))
  # }
  
  # sigcls <- 1/sqrt(clique.size)#rep(10, length(cliques))#runif(length(cliques))
  # #sigmals <- #unlist(sapply(1:length(cliques), FUN=function(i){rep(sigcls[i], length(cliques[[i]]))}))
  # 
  # meancls <- rep(0, length(cliques))
  # 
  # Agamvec <- colSums(Amat)
  # dvec    <- Agamvec #/ sigmals^2
  # 
  # theta <- rep(0.5, length(cliques))
  # a = c = 1
  # v0 = 0.1
  # v1 = 100
  
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
  #BS1 <- bsplineS((1:7)/7, breaks = seq(0, 1, 1/knot))
  z <- as.vector(y - array(X1%*%b))
  for(itr in 1:Total_itr){
    ##Sample random variables
    #eta <- unlist(lapply(1:n, Ei))
    
    #ret <- sapply(1:length(cliques), FUN=function(i) {gamvec[i]*sum((b[cliques[[i]]]-meancls[i])^2*invsigmals[cliques[[i]]])})
    
    #b <- as.vector(solve(t(X)%*%X)%*%t(X)%*%(y-u))
    z <- as.vector(y - array(X1%*%b))
    
    #Mixed effect other parts
    iO <- solve(Omega)
    Ts <- R <- C <- 0
    mu <- u <- NULL
    for (i in 1:n) {
      row.i <- which(id==i)
      Xi <- X1[row.i,]
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
    
    num <- sum((yred - array(X1%*%b))^2) + nu*lambda #+ sum(b^2)/100#sum(unlist(ret)) #sum(b^2*dvec)
    dem <- length(y)+nu + 1 #+sum(colSums(Amat))
    
    sig = sqrt(num/dem)
    
    #solve(crossprod(X) + diag(dvec))#
    #if(itr>1){
    #ind <- which(dvec>0.01) 
    #}
    M <- solve(crossprod(X1[,ind]) + sig^2*diag(p)/100)%*%crossprod(X1[,ind], yred)#diag(1/dvec) - diag(1/dvec) %*% t(X) %*% solve(diag(n)+ X1%*%diag(1/dvec)%*%t(X)) %*% X %*% diag(1/dvec)
    ##out <- nnls::nnls(crossprod(X1[,ind])+1/100, crossprod(X1[,ind], yred))#+colSums(Amat * gamvec * meancls/sigcls^2))#M %*% crossprod(X[,ind], y)
    b[ind] <- M#out$x
    #b[-ind] <- 0
    
    #bconti   <- unlist(sapply(1:length(cliques), FUN=function(i) {sum((b[cliques[[i]]])^2)/sigcls[i]^2}))
    
    
    # Amat <- matrix(0, length(cliques), p)
    # 
    # for(i in 1:nrow(Amat)){
    #   Amat[i, cliques[[i]]] <- invsigmals[[i]]
    # }
    # 
    # 
    # Agamvec <- colSums(Amat)
    # dvec    <- Agamvec
    # 
    # for(i in 1:length(cliques)){
    #   bi <- b[cliques[[i]]]#-meancls[i]
    #   pstar   <- theta[i]*dnorm(bi, sd = sqrt(v1))/(theta[i]*dnorm(bi, sd = sqrt(v1))+(1-theta[i])*dnorm(bi, sd = sqrt(v0)))
    #   invsigmals[[i]] <- ((1-pstar)/v0 + pstar/v1)
    #   theta[i] <- (sum(pstar) + a-1)/(a+c+p-2)
    # }
    
    #plot(BS1 %*% b)
    #print(which.max(theta))
    #bvec   <- unlist(sapply(1:length(cliques), FUN=function(i) {mean((b[cliques[[i]]])^2)}))#, sd = sigmals[cliques[[i]]]*sig))})
    #plot(bvec)
    #plot(bvec)
    #plot(theta)
  }
  #return(mean((ftf - BS1 %*% b)^2))
  return(list(sigma = sig, esti=BS1 %*% b))
  #return(list(betaest = b, thetals = theta))
}