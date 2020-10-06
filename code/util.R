### fit the NB ash model based on polya gamma augmentation, with r fixed and large
### no random effects is included in the model
### x: the n x 1 vector of observed counts
### sigma2: the K x 1 vector of prior variances used in the ash prior
### r: the dispersion parameter of the negative binomial distribution, and treated as fixed in this function
### init: a list containing initial values for mu and pi, and can be empty
### maxiter: the maximum number of iterations
### verbose: logical indicator whether to print ELBO at each iteration
nb_ash_pg <- function(x, sigma2, r=1000, init=NULL, maxiter=100, verbose=FALSE){
  n <- length(x)
  K <- length(sigma2)
  
  mu <- init$mu
  if(is.null(mu)){
    mu <- log(sum(x)) - log(sum(x+r))
  }
  
  pi <- init$pi
  if(is.null(pi)){
    pi <- rep(1/K, K)
  }
  
  x_mat <- x %*% t(rep(1,K))
  sigma2_mat <- rep(1, n) %*% t(sigma2)
  m <- matrix(0, nrow=n, ncol=K)
  v2 <- matrix(1, nrow=n, ncol=K)
  
  ELBOs <- c()
  mu.seq <- c()
  pi.seq <- matrix(NA, nrow=K, ncol=0)
  
  for(iter in 1:maxiter){
    # update posterior mean of polya gamma latent variables
    psi <- sqrt(mu^2 + m^2 + v2 + 2*mu*m)
    xi <- (x_mat + r)*tanh(psi/2)/(2*psi)
    xi <- pmin(pmax(xi, 1e-6), 1e6)
    xi_inv <- 1/xi
    
    # update posterior mean m_ik and variance v2_ik of beta
    v2 <- 1/(xi + 1/sigma2_mat)
    v2 <- pmax(v2, 1e-6)
    m <- sigma2_mat*(-mu + (x_mat-r)*xi_inv/2)/(sigma2_mat + xi_inv)
    
    # update posterior mean of z_ik
    ELBO.local <- x_mat*(mu+m)-(x_mat+r)*exp(mu+m) - 0.5*(log(sigma2_mat) - log(v2) + (m^2+v2)/sigma2_mat - 1) 
    ELBO.cen <- ELBO.local - apply(ELBO.local, 1, max)
    zeta <- t(t(exp(ELBO.cen)) * pi)
    zeta <- zeta*(1/rowSums(zeta)) 
    zeta <- pmax(zeta, 1e-15)
    
    # update mu
    mu.new <- log(sum(x)) - log(sum(zeta*(x_mat+r)*exp(m)))
    diff.mu <- mu.new - mu
    mu <- mu.new
    mu.seq <- c(mu.seq, mu)
    
    # update pi
    pi.new <- colMeans(zeta)
    pi.new <- pmax(pi.new, 1e-6)
    diff.pi <- pi.new - pi
    pi <- pi.new
    pi.seq <- cbind(pi.seq, pi)
    
    # compute overall ELBO
    ELBO.local <- x_mat*(mu+m)-(x_mat+r)*exp(mu+m) - 0.5*(log(sigma2_mat) - log(v2) + (m^2+v2)/sigma2_mat - 1)
    pi_mat <- rep(1, n) %*% t(pi)
    ELBO.overall <- sum(zeta*(log(pi_mat) + ELBO.local - log(zeta)))
    ELBOs <- c(ELBOs, ELBO.overall) 
    
    if(verbose){
      print("iter         ELBO")
      print(sprintf("%d:    %f", iter, ELBO.overall))
    }
    
    # if(abs(diff.mu) < tol & max(abs(diff.pi)) < tol) break
  }
  
  return(list(mu=mu, mu.seq=mu.seq, diff.mu=diff.mu, pi=pi, diff.pi=diff.pi, pi.seq=pi.seq, m=m, v2=v2, zeta=zeta, ELBO=ELBOs))
}