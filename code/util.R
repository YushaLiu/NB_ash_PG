nb_ash_pg <- function(x, sigma2, r=1000, init=NULL, maxiter=100, tol=1e-6, verbose=FALSE){
  n <- length(x)
  K <- length(sigma2)
  
  mu <- init$mu
  if(is.null(mu)){
    mu <- log(sum(x)) - log(sum(x+r))
  }
  
  s2 <- init$s2
  if(is.null(s2)){
    s2 <- 1
  }
  
  pi <- init$pi
  if(is.null(pi)){
    pi <- rep(1/K, K)
  }
  
  x_mat <- x %*% t(rep(1,K))
  sigma2_mat <- rep(1, n) %*% t(sigma2)
  m <- matrix(0, nrow=n, ncol=K)
  v2 <- matrix(1, nrow=n, ncol=K)
  theta <- matrix(0, nrow=n, ncol=K)
  w2 <- matrix(1, nrow=n, ncol=K)

  ELBOs <- c()
  mu.seq <- c()
  s2.seq <- c()
  pi.seq <- matrix(NA, nrow=K, ncol=0)
   
  for(iter in 1:maxiter){
    # update posterior mean of polya gamma latent variables
    psi <- sqrt(mu^2 + m^2 + v2 + theta^2 + w2 + 2*mu*m + 2*mu*theta + 2*m*theta)
    xi <- (x_mat + r)*tanh(psi/2)/(2*psi)
    xi <- pmin(pmax(xi, 1e-6), 1e6)
    xi_inv <- 1/xi
    
    # update posterior mean m_ik and variance v2_ik for beta
    v2 <- 1/(xi + 1/sigma2_mat)
    v2 <- pmax(v2, 1e-6)
    m <- sigma2_mat*(-mu - theta + (x_mat-r)*xi_inv/2)/(sigma2_mat + xi_inv)
    
    # update posterior mean theta_ik and variance w2_ik for u
    w2 <- 1/(xi + 1/s2)
    w2 <- pmax(w2, 1e-6)
    theta <- s2*(-mu - m + (x_mat-r)*xi_inv/2)/(s2 + xi_inv)
    
    # update posterior mean of z_ik
    ELBO.local <- x_mat*(mu+m+theta)-(x_mat+r)*exp(mu+m+theta) -
      0.5*(log(sigma2_mat) - log(v2) + (m^2+v2)/sigma2_mat - 1) - 0.5*(log(s2) - log(w2) + (theta^2+w2)/s2 - 1)
    ELBO.cen <- ELBO.local - apply(ELBO.local, 1, max)
    zeta <- t(t(exp(ELBO.cen)) * pi)
    zeta <- zeta*(1/rowSums(zeta)) 
    zeta <- pmax(zeta, 1e-15)
    
    # update mu
    mu.new <- log(sum(x)) - log(sum(zeta*(x_mat+r)*exp(m+theta)))
    diff.mu <- mu.new - mu
    mu <- mu.new
    mu.seq <- c(mu.seq, mu)
    
    # update s2
    s2.new <- sum(zeta*(theta^2 + w2))/n
    s2.new <- pmax(s2.new, 1e-6)
    diff.s2 <- s2.new - s2
    s2 <- s2.new
    s2.seq <- c(s2.seq, s2)
    
    # update pi
    pi.new <- colMeans(zeta)
    pi.new <- pmax(pi.new, 1e-6)
    diff.pi <- pi.new - pi
    pi <- pi.new
    pi.seq <- cbind(pi.seq, pi)
    
    # compute overall ELBO
    ELBO.local <- x_mat*(mu+m+theta)-(x_mat+r)*exp(mu+m+theta) -
      0.5*(log(sigma2_mat) - log(v2) + (m^2+v2)/sigma2_mat - 1) - 0.5*(log(s2) - log(w2) + (theta^2+w2)/s2 - 1)
    pi_mat <- rep(1, n) %*% t(pi)
    ELBO.overall <- sum(zeta*(log(pi_mat) + ELBO.local - log(zeta)))
    ELBOs <- c(ELBOs, ELBO.overall) 
    
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO.overall))
    }
    
    # if(abs(diff.mu) < tol & max(abs(diff.pi)) < tol) break
  }
  
  return(list(mu=mu, mu.seq=mu.seq, diff.mu=diff.mu, s2=s2, s2.seq=s2.seq, diff.s2=diff.s2, pi=pi, diff.pi=diff.pi, pi.seq=pi.seq,
              m=m, v2=v2, theta=theta, w2=w2, zeta=zeta, ELBO=ELBOs))
}


