# Naive-Bayessian-Regression



How to Execute these functions in R or Rstudio:

First, load the neccesary functions
```r
source("scr/br.source")
source("scr/getOverlapping.source")
```

Simulate data
```r 
n=1000 #ndatos
p=10   #covariates
x1<-floor(runif(n,0,3)) #simulate genotypes
        X<-matrix(nrow=n,ncol=p,NA)
        for(i in 1:p){X[,i]<-floor(runif(n)*3 )}

beta<-rep(0,p)

beta[1]<-rnorm(1,10,10) #allele subst effects
beta[3]<-rnorm(1,10,10) #allele subst effects
beta[7]<-rnorm(1,10,10) #allele subst effects
mu=20
y<-mu+X%*%beta+rnorm(n,0,sd=50) #generate data
```

Now we can perform the Bayesian regression of the phenotype *y* on the covariates with the *br function*, and (optional) hyperparameteres for a priori residual variance of v=10 and S=2500. (Default values are v=1 and S=1 if no hyperparameters are provided).
```r
fit<-br(y,X,niter=10000,burnin=1000,permutation = F,v0=10,S0=2500)
```

you can check the results as
```r
#Solutions
##Posterior means for intercept
fit$itcp_Hat
##Posterior means for effects
fit$beta_Hat
##Posterior distribution for intercept
fit$itcp
##Posterior distribution for effects
fit$beta_pd[,1]
##Posterior distribution for ve
fit$ve
##Posterior mean for ve
fit$ve_Hat

burnin=fit$burnin;niter=fit$niter
#Plot posterior distribution for effects
plot(density(fit$beta_pd[(burnin+1):niter,1]))
#Plot posterior distribution for ve
plot(density(fit$ve[(burnin+1):niter]))

#Plot chain for ve
plot(fit$ve,type="l")
```

You can also perform a Bayesian permutation test on the data with the same model
```r
fit.perm<-br(y,X,niter=10000,burnin=1000,permutation = T)
```

Then we obtain the overlapping proportion between the posterior distribution of all covariates from the original model and the permutated model using the *getOverlapping function*.
```r 
significancy<-getOverlapping(fit,fit.perm)
```
We can extract all covariates with an overlapping lower than a given proportion (e.g.0.05):
```r
which(significancy$overlapping<0.05)
```
These are the significant covariates from the permutation test. Then we can plot these distributions for a given factor (identified as the order of its column in X) as follows.

```r
#Plot overlapping for factor
factor=6
df2 <- data.frame(x = c(significancy$d1dens[[factor]]$x, significancy$d1dens[[factor]]$x,significancy$d1dens[[factor]]$x), 
                  y = c(significancy$d1dens[[factor]]$y, significancy$joint[[factor]],significancy$d2dens[[factor]]$y),
                  Data = rep(c("Model", "overlap", "null modell"), each = length(significancy$d1dens[[factor]]$x)))

ggplot(df2, aes(x, y, fill = Data)) + 
    geom_area(position = position_identity(), color = "black") +
    scale_fill_brewer(palette = "Pastel2") +
    theme_bw()
```

