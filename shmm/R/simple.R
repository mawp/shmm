# Todo:
# DONE 1. Construct generator inside TMB
# DONE 2. Implement uniformization algorithm, input a vector that will get multiplied and output it afterward
# DONE 3. Make 2D version
# DONE 4. Implement filter
# DONE 5. Implement simulator and data likelihood calculator for simple example
# DONE 6. Implement likelihood function and estimate
#  7. Implement smoother. How: use a flag as a parameter fixed using map
#  8. Exploit symmetry of G somehow? see http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html
#     Symmetry is broken if advection is used
#  9. Land in simulator and optimised code that doesn't include in P the inaccessible land cells
# 10. Think about how to choose m, this depends on the expected range of D and has a large influence on speed
#     What happens if the used m is smaller than the required m? will estimation work?
rm(list=ls())
require(TMB)
source('shmmfuns.R')
compile('../src/simple.cpp')
scriptname <- "../src/simple"
dyn.load(dynlib(scriptname))
inp <- list()

# Parameters
inp$ini$sdx <- 6
inp$ini$sdy <- 5
osd <- 1 # Observation error sd
inp$dt <- 1
nt <- 30

# Simulate movements
set.seed(654)
Xx <- cumsum(rnorm(nt, 0, inp$ini$sdx*sqrt(inp$dt)))
Xy <- cumsum(rnorm(nt, 0, inp$ini$sdy*sqrt(inp$dt)))
Yx <- Xx + rnorm(nt, 0, osd)
Yy <- Xy + rnorm(nt, 0, osd)

# SHMM stuff
require(Matrix)
# X
xmin <- min(Xx) - 0.1*diff(range(Xx))
xmax <- max(Xx) + 0.1*diff(range(Xx))
inp$dx <- 1
xx <- seq(xmin, xmax, by=inp$dx)
nx <- length(xx)
# Y
ymin <- min(Xy) - 0.1*diff(range(Xy))
ymax <- max(Xy) + 0.1*diff(range(Xy))
inp$dy <- 1
yy <- seq(ymin, ymax, by=inp$dy)
ny <- length(yy)
# Land
indsy <- which(yy > 6 & yy < 8)
indsx <- which(xx > 1 & xx < 4)
dummy <- matrix(0, nx, ny)
#dummy[indsx, indsy] <- 1
land <- which(dummy==1)
cat('nx: ', nx, ' ny: ', ny, ' nt: ', nt, ' no states: ', nx*ny, '\n')

tic <- Sys.time()

# Data likelihood
n <- nx*ny
inp$datlik <- array(0, dim=c(nt, nx, ny))
for (i in 1:nt) inp$datlik[i, , ] <- outer(dnorm(xx, Yx[i], osd), dnorm(yy, Yy[i], osd))

# Fit SHMM
rep <- fit.shmm(inp)


# Get parameter estimates
ests <- opt$par
sds <- sqrt(diag(rep$cov.fixed))
sdxest <- sqrt(2*exp(ests[1]))*dx
sdxll <- sqrt(2*exp(ests[1] - 2*sds[1]))*dx
sdxul <- sqrt(2*exp(ests[1] + 2*sds[1]))*dx
sdyest <- sqrt(2*exp(ests[2]))*dy
sdyll <- sqrt(2*exp(ests[2] - 2*sds[2]))*dy
sdyul <- sqrt(2*exp(ests[2] + 2*sds[2]))*dy

pars <- rbind(c(sdxest, sdxll, sdxul, sdx), c(sdyest, sdyll, sdyul, sdy))
colnames(pars) <- c('est', 'll', 'ul', 'true')
rownames(pars) <- c('sdx', 'sdy')

# Get distributions
phi <- obj$report()$phi
psi <- obj$report()$psi
pred <- obj$report()$pred

# Calculate tracks
XX <- as.vector(matrix(xx, nx, ny))
YY <- as.vector(matrix(yy, nx, ny, byrow=TRUE))
xest <- numeric(nt)
yest <- numeric(nt)
for(t in 1:nt){
    xest[t] <- sum(phi[t, ] * XX)
    yest[t] <- sum(phi[t, ] * YY)
}

toc <- Sys.time()

par(mfrow=c(2, 1))
plot(xest, typ='l')
lines(Xx, col=3)
plot(yest, typ='l')
lines(Xy, col=3)

plot(xest, yest, typ='l')
lines(Xx, Xy, typ='l', col=3)
points(Yx, Yy, typ='p', col=2)


for(t in 1:nt){
    matpred <- matrix(pred[t, ], nx, ny)
    matphi <- matrix(phi[t, ], nx, ny)
    par(mfrow=c(2, 2))
    image(xx, yy, matpred, main='pred')
    lines(Xx, Xy, typ='b')
    lines(Yx, Yy, typ='b', col=3)
    image(xx, yy, matphi, main='phi')
    lines(Xx, Xy, typ='b')
    lines(Yx, Yy, typ='b', col=3)
    Sys.sleep(0.2)
}
