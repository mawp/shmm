# Todo:
# DONE 1. Construct generator inside TMB
# DONE 2. Implement uniformization algorithm, input a vector that will get multiplied and output it afterward
# 3. Implement filter
# 4. Implement likelihood function
# 5. Implement smoother
require(TMB)
compile('../src/simple.cpp')
scriptname <- "../src/simple"
dyn.load(dynlib(scriptname))

## Test data
set.seed(123)
y <- rep(1900:2010,each=2)
year <- factor(y)
quarter <- factor(rep(1:4,length.out=length(year)))
period <- factor((y > mean(y))+1)
## Random year+quarter effect, fixed period effect:
B <- model.matrix(~year+quarter-1)
A <- model.matrix(~period-1)
B <- as(B,"dgTMatrix")
A <- as(A,"dgTMatrix")
u <- rnorm(ncol(B)) ## logsdu=0
beta <- rnorm(ncol(A))*100
eps <- rnorm(nrow(B),sd=1) ## logsd0=0
x <- as.numeric(A%*%beta+B%*%u+eps)

# SHMM stuff
require(Matrix)
mu <- 1
n <- 800
S <- Matrix(0, n, n)
doS <- Diagonal(n-1, 1)
S[1:(n-1), 2:n] <- doS
S[2:n, 1:(n-1)] <- S[2:n, 1:(n-1)] + doS
diag(S) <- -apply(S, 1, sum)

dt <- 1
mu <- 8
AA <- S*mu
F <- max(abs(AA))
#I <- Diagonal(n)
I <- Matrix(0, n, n)
diag(I) <- rep(1, n)
P <- AA/F + I
m <- 37
m <- ceiling(F*dt + 4*sqrt(F*dt) + 5)

pvec <- Matrix(0, 1, n)
pvec[1, round(n/2)] <- 1

P <- as(P, 'dgTMatrix')
I <- as(I, 'dgTMatrix')
pvec <- as(pvec, 'dgTMatrix')

## Fit model
map <- NULL
get.pout <- function(){
obj <- MakeADFun(data=list(x=x, B=B, A=A, P=P, I=I, pvec=pvec, mu=mu, F=F, dt=dt, m=m),
                 parameters=list(u=u*0+.1,beta=beta*0+.1,logsdu=.1,logsd0=0.1),
                 random="u",
                 DLL='simple'
                 )
    return(obj$report()$pout2)
}



#newtonOption(smartsearch=FALSE)
#system.time(obj$fn(obj$par))
#obj$gr(obj$par)
#obj$control <- list(parscale=obj$par*0+1e-1,trace=10)
#obj$hessian <- TRUE
#opt <- do.call("optim",obj)

#Pout <- obj$report()$P
#pout <- obj$report()$pout

#Sout <- obj$report()$S
#fact <- obj$report()$fact
#pt <- obj$report()$pt



# --- Uniformization in R ---
uniformization <- function(A,dt){
  ## Uniformization is an efficient way to compute the matrix exponential for a generator matrix
  ## See Grassmann 1977 - Transient solutions in Markovian queueing systems
  
  N <- dim(A)[1]
  # Find the numerical largest rate
  F <- max(abs(diag(A)))
  # Calculate number of iterations based expression in Grassmann 1977, eq 10
  m <- ceiling(F*dt + 4*sqrt(F*dt) + 5)
  # Insert warning if m>140 ??
  I <- Diagonal(N)
  P <- A/F + I # Create sub-stochastic matrix (eq 8)

  S <- I
  pt <- I
  FPdt <- F*P*dt
  for(i in 1:m){
    S <- S %*% FPdt
    fact <- exp(lgamma(i+1))
    pt <- pt + S/fact
  }
  # Multiply by the constant
  pt <- pt * exp(-F*dt)
  pt
}

G <- S*mu

system.time(pout2 <- get.pout())
#system.time(obj$fn(obj$par))
system.time(pout2b <- pvec %*% uniformization(G, dt))
system.time(pout2c <- pvec %*% expm(G*dt))
sum(pout2 - pout2b)
sum(pout2 - pout2c)
