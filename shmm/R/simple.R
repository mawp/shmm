# Todo:
# DONE 1. Construct generator inside TMB
# DONE 2. Implement uniformization algorithm, input a vector that will get multiplied and output it afterward
# DONE 3. Make 2D version
# DONE 4. Implement filter
# DONE 5. Implement simulator and data likelihood calculator for simple example
# DONE 6. Implement likelihood function and estimate
# 7. Implement smoother. How: use a flag as a parameter fixed using map
# 8. Exploit symmetry of G somehow? see http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html
#    Symmetry is broken if advection is used
# 9. Land in simulator and optimised code that doesn't include in P the inaccessible land cells
rm(list=ls())
require(TMB)
source('shmmfuns.R')
compile('../src/simple.cpp')
scriptname <- "../src/simple"
dyn.load(dynlib(scriptname))

# Parameters
sdx <- 3
sdy <- 4
osd <- 1
dt <- 1
nt <- 30

# Simulate movements
set.seed(654)
Xx <- cumsum(rnorm(nt, 0, sdx*sqrt(dt)))
Xy <- cumsum(rnorm(nt, 0, sdy*sqrt(dt)))
Yx <- Xx + rnorm(nt, 0, osd)
Yy <- Xy + rnorm(nt, 0, osd)

# SHMM stuff
require(Matrix)
# X
xmin <- min(Xx) - 0.1*diff(range(Xx))
xmax <- max(Xx) + 0.1*diff(range(Xx))
dx <- 1
xx <- seq(xmin, xmax, by=dx)
nx <- length(xx)
# Y
ymin <- min(Xy) - 0.1*diff(range(Xy))
ymax <- max(Xy) + 0.1*diff(range(Xy))
dy <- 1
yy <- seq(ymin, ymax, by=dy)
ny <- length(yy)
# Land
indsy <- which(yy > 6 & yy < 8)
indsx <- which(xx > 1 & xx < 4)
dummy <- matrix(0, nx, ny)
#dummy[indsx, indsy] <- 1
land <- which(dummy==1)


tic <- Sys.time()

# Data likelihood
n <- nx*ny
datlik <- matrix(0, nt, n)
for(i in 1:nt) datlik[i, ] <- as.vector(outer(dnorm(xx, Yx[i], osd), dnorm(yy, Yy[i], osd)))

# TMB
# Close diagonals (jump north-south)
Sns <- make.ns(n, nx, land)
# Far diagonals (jump east-west)
Sew <- make.ew(n, nx, land)
logDx <- log(0.5 * (sdx/dx)^2)  # East west move rate (diffusion)
logDy <- log(0.5 * (sdy/dy)^2)  # North south move rate (diffusion)

# Generator
AA <- make.generator(nx, ny, dx, dy, sdx, sdy, land=land)
F <- max(abs(AA))
I <- Matrix(0, n, n)
diag(I) <- rep(1, n)
P <- AA/F + I
m <- ceiling(F*dt + 4*sqrt(F*dt) + 5) # Expression from Grassmann
#G <- uniformization(AA, dt)

P <- as(P, 'dgTMatrix')
I <- as(I, 'dgTMatrix')
Sew <- as(Sew, 'dgTMatrix')
Sns <- as(Sns, 'dgTMatrix')

# Create objective function
lgam <- lgamma(2:(m+2))
data <- list(datlik=datlik, I=I, dt=dt, m=m, Sns=Sns, Sew=Sew, lgam=lgam)
pars <- list(logDx=logDx, logDy=logDy)
obj <- MakeADFun(data=data, parameters=pars, random=NULL, DLL='simple') # This one requires memory

# Estimate
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
system.time(rep <- sdreport(obj))

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
