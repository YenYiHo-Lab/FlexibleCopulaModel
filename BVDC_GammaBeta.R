############################################################################
##### Bivariate dynamic correlation

library(msm) # for univariate truncated normal

############################################################################
##### Log-likelihood or conditionals needed in the MCMC updates

##### CDF of Y1 ~ Zero-inflated Gamma 
##### [ZIG(mean = mu, dispersion = lambda, zero-inflation = p)]

pgam = function(y,mu,lambda,p){(1-p)*pgamma(y,mu*lambda,lambda) + p}

##### Log-likelihood of Y1

llike1 = function(y,z,mu,lambda,p,rho){
	id0 = which(y==0)
	id1 = which(y!=0)

	u = pgam(y[id1],mu[id1],lambda,p[id1])
	u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
	llik = sum(-.5*(rho[id1]^2*qnorm(u)^2 - 2*qnorm(u)*z[id1,2]*rho[id1])/
		(1-rho[id1]^2),na.rm=T)
	llik = llik + sum(log(1-p[id1]) + 
			dgamma(y[id1],mu[id1]*lambda,lambda,log=T),na.rm=T)

	if(length(id0)>0){
		llik = llik + sum(pnorm(qnorm(p[id0]),mean=rho[id0]*z[id0,2],
			sd=sqrt(1-rho[id0]^2),log=T),na.rm=T)
	}
return(llik)}

##### Log-likelihood of Y2 ~ Beta(mu*phi, (1-mu)*phi)

llike2 = function(y,z,mu,phi,rho){
    u = pbeta(y,mu*phi,(1-mu)*phi)
    u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
    llik = sum(-.5*(rho^2*qnorm(u)^2 - 2*qnorm(u)*z[,1]*rho)/
		(1-rho^2),na.rm=T)
    llik = llik + sum(dbeta(y,mu*phi,(1-mu)*phi,log=T),na.rm=T)
return(llik)}

##### Log-likelihood of (Z1, Z2) given rho

mvdens.z = function(z,rho){
	llik = sum(-.5*(z[,1]^2 - 2*rho*z[,1]*z[,2] + z[,2]^2)/(1-rho^2) - 
		.5*log(1-rho^2),na.rm=T)
return(llik)}

############################################################################
##### Functions that update parameters in MCMC

##### Update parameters associated with 
##### Y1 ~ ZIG(mu,lambda,p)

# 1. Update beta such that mu[i] = exp( x[i]'beta )

beta.update = function(beta,y,x,z,mu.old,lambda,p,rho,B=200,h=200){
    n = nrow(x); d = ncol(x)
    it = if(is.vector(beta)){1} else {nrow(beta)}
    beta.old = if(is.vector(beta)){beta} else{beta[it,]}
    
    if(it <= B){
        cov.beta = diag(d)*.01
    } else {cov.beta = cov(beta[(it-h+1):it,]) + .0001*diag(d)}
    beta.star = beta.old + c(rnorm(d)%*%chol(cov.beta))
    
    mu.star = c(exp(x%*%beta.star))
    
    log.r = llike1(y,z,mu.star,lambda,p,rho) + 
	  sum(dnorm(beta.star,0,100,log=T)) -  
        llike1(y,z,mu.old,lambda,p,rho) - 
	  sum(dnorm(beta.old,0,100,log=T))
    
    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(list(beta=beta.star,mu=mu.star))
    } else{return(list(beta=beta.old,mu=mu.old))}
}

# 2. Update lambda > 0

lambda.update = function(lambda,y,z,mu,p,rho,B=200,h=200){
    it = length(lambda); lambda.old = lambda[it]
    sd.lambda = ifelse(it <= B, .1, sd(log(lambda[(it-h+1):it])))
    lambda.star = exp(rnorm(1,log(lambda.old),sd.lambda))

    log.r = llike1(y,z,mu,lambda.star,p,rho) + 
        dgamma(lambda.star,.001,.001,log=T) - log(lambda.old) - 
        llike1(y,z,mu,lambda.old,p,rho) - 
        dgamma(lambda.old,.001,.001,log=T) + log(lambda.star)
    
    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(lambda.star)
    } else {return(lambda.old)}
}

# 3. Update eta such that log(p[i]/(1-p[i])) = x[i]'eta

eta.update = function(eta,y,x,z,mu,lambda,p.old,rho,B=200,h=200){
    n = nrow(x); d = ncol(x)
    it = if(is.vector(eta)){1} else{nrow(eta)}
    eta.old = if(is.vector(eta)){eta} else{eta[it,]}
    
    if(it <= B){
        cov.eta = diag(d)*.01
    } else{cov.eta = cov(eta[(it-h+1):it,]) + .0001*diag(d)}
    eta.star = eta.old + c(rnorm(d)%*%chol(cov.eta))
    
    p.star = c(exp(x%*%eta.star)/(1 + exp(x%*%eta.star)))
    
    log.r = llike1(y,z,mu,lambda,p.star,rho) + 
	  sum(dnorm(eta.star,0,100,log=T)) -
        llike1(y,z,mu,lambda,p.old,rho) - 
	  sum(dnorm(eta.old,0,100,log=T))
    
    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(list(eta=eta.star,p=p.star))
    } else{return(list(eta=eta.old,p=p.old))}
}

##### Update parameters associated with 
##### Y2 ~ Beta(mu*phi, (1-mu)*phi)

# 4. Update theta such that log(mu[i]/(1-mu[i])) = x[i]'theta

theta.update = function(theta,y,x,z,mu.old,phi,rho,B=200,h=200){
    n = nrow(x); d = ncol(x)
    it = if(is.vector(theta)){1} else{nrow(theta)}
    theta.old = if(is.vector(theta)){theta} else{theta[it,]}
    
    if(it <= B){
        cov.theta = diag(d)*.01
    } else{cov.theta = cov(theta[(it-h+1):it,]) + .0001*diag(d)}
    theta.star = theta.old + c(rnorm(d)%*%chol(cov.theta))
    
    mu.star = c(exp(x%*%theta.star)/(1 + exp(x%*%theta.star)))
    
    log.r = llike2(y,z,mu.star,phi,rho) + 
        sum(dnorm(theta.star,0,100,log=T)) -
        llike2(y,z,mu.old,phi,rho) - 
        sum(dnorm(theta.old,0,100,log=T))
    
    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(list(theta=theta.star,mu=mu.star))
    } else{return(list(theta=theta.old,mu=mu.old))}
}

# 5. Update phi > 0

phi.update = function(phi,y,z,mu,rho,B=200,h=200){
    it = length(phi); phi.old = phi[it]
    sd.phi = ifelse(it <= B, .1, sd(log(phi[(it-h+1):it])))
    phi.star = exp(rnorm(1,log(phi.old),sd.phi))
    
    log.r = llike2(y,z,mu,phi.star,rho) + 
        dgamma(phi.star,.01,.01,log=T) - log(phi.old) -
        llike2(y,z,mu,phi.old,rho) - 
        dgamma(phi.old,.01,.01,log=T) + log(phi.star)
    
    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(phi.star)
    } else {return(phi.old)}
}

##### 6. Update (Z1, Z2)

# 6(a). Update Z1 | Y1, Z2

z1.update = function(y,z,mu,lambda,p,rho){
    id0 = which(y==0)
    id1 = which(y!=0)
    z.draw = rep(0,length(y))

    u = pgam(y[id1],mu[id1],lambda,p[id1])
    u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
    z.draw[id1] = qnorm(u)

    if(length(id0)>0){
	z.draw[id0] = rtnorm(length(id0),mean=rho[id0]*z[id0],
		sd=sqrt(1-rho[id0]^2),upper=c(qnorm(p[id0])))
    }
return(z.draw)}

# 6(b). Update Z2 | Y2

z2.update = function(y,mu,phi){
    u = pbeta(y,mu*phi,(1-mu)*phi)
    u = ifelse(abs(u-.5) < .4999,u,sign(u-.5)*.4999+.5)
    z.draw = qnorm(u)
return(z.draw)}

##### 7(a). Update tau such that log( (1+rho[i])/(1-rho[i]) ) = x[i]'tau

tau.update = function(tau,z,x,rho.old,B=200,h=200){
    n = nrow(x); d = ncol(x)
    it = if(is.vector(tau)){1} else {nrow(tau)}
    tau.old = if(is.vector(tau)){tau} else{tau[it,]}
    
    if(it<=B){
        cov.tau = diag(d)*.01
    } else{cov.tau = cov(tau[(it-h+1):it,]) + .0001*diag(d)}
    tau.star = tau.old + c(rnorm(d)%*%chol(cov.tau))
    
    rho.star = c((exp(x%*%tau.star)-1)/(exp(x%*%tau.star)+1))
    
    log.r = mvdens.z(z,rho.star) + sum(dnorm(tau.star,0,100,log=T)) - 
        mvdens.z(z,rho.old) - sum(dnorm(tau.old,0,100,log=T))
    
    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(list(tau=tau.star,rho=rho.star))
    } else{return(list(tau=tau.old,rho=rho.old))}
}

##### 7(b). Update rho if it is a scalar

rho.update = function(z,rho,B=200,h=200){
    it = length(rho); rho.old = rho[it]
    x = log((1+rho)/(1-rho)); x.old = x[it]
    sd.x = ifelse(it <= B, .1, sd(x[(it-h+1):it]))
    x.star = rnorm(1,x.old,sd.x)
    rho.star = (exp(x.star)-1)/(exp(x.star)+1)

    log.r = mvdens.z.2(z,rho.star) - mvdens.z.2(z,rho.old) + 
		2*log(exp(x.star)+1) - 
		2*log(exp(x.old)+1) + x.old - x.star

    if((!is.na(log.r)) & (log.r > log(runif(1)))){
        return(rho.star)
    } else {return(rho.old)}
}

#############################################################################
##### MCMC output

mcmc.fun = function(y,x,x.p,n.mcmc=20000,burn=5000,thin=50){
    beta = theta = tau = matrix(0,nrow=n.mcmc,ncol=ncol(x))
    eta = matrix(0,nrow=n.mcmc,ncol=ncol(x.p))
    lambda = phi = rep(1,n.mcmc)

    mu1 = mu2 = c(exp(x%*%(beta[1,])))
    p = c(exp(x.p%*%(eta[1,]))/(1+exp(x.p%*%(eta[1,]))))
    rho = c((exp(x%*%(tau[1,]))-1)/(1+exp(x%*%(tau[1,]))))
    z = matrix(rnorm(n*2),n,2)

    for(i in 2:n.mcmc){
	  # print(i);flush.console()
	  beta.mcmc = beta.update(beta[1:(i-1),],y[,1],x,z,mu1,
		lambda[i-1],p,rho)
	  beta[i,] = beta.mcmc$beta
	  mu1 = beta.mcmc$mu

	  lambda[i] = lambda.update(lambda[1:(i-1)],y[,1],z,mu1,p,rho)

	  eta.mcmc = eta.update(eta[1:(i-1),],y[,1],x.p,z,mu1,lambda[i],p,rho)
	  eta[i,] = eta.mcmc$eta
	  p = eta.mcmc$p

	  z[,1] = z1.update(y[,1],z[,2],mu1,lambda[i],p,rho)

	  theta.mcmc = theta.update(theta[1:(i-1),],y[,2],x,z,mu2,phi[i-1],rho)
	  theta[i,] = theta.mcmc$theta
	  mu2 = theta.mcmc$mu

	  phi[i] = phi.update(phi[1:(i-1)],y[,2],z,mu2,rho)

	  z[,2] = z2.update(y[,2],mu2,phi[i])
	
	  # z = z.update(y,mu1,lambda[i],p,mu2,phi[i],rho)

	  tau.mcmc = tau.update(tau[1:(i-1),],z,x,rho)
	  tau[i,] = tau.mcmc$tau
	  rho = tau.mcmc$rho
    }
    samp = seq(burn+1,n.mcmc,by=thin)
    mcmc.out = cbind(beta,eta,theta,tau,lambda,phi)[samp,]

    colnames(mcmc.out) = c(paste0("beta",0:(ncol(beta)-1)),
	paste0("eta",0:(ncol(eta)-1)),
	paste0("theta",0:(ncol(theta)-1)),paste0("tau",0:(ncol(tau)-1)),
	"lambda","phi")

return(mcmc.out)}
