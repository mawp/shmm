# Todo:
# DONE 1. Construct generator inside TMB
# DONE 2. Implement uniformization algorithm, input a vector that will get multiplied and output it afterward
# DONE 3. Make 2D version
# DONE 4. Implement filter
# DONE 5. Implement simulator and data likelihood calculator for simple example
# DONE 6. Implement likelihood function and estimate
# 7. Implement smoother
# 8. Exploit symmetry of G somehow? see http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html
#    Symmetry is broken if advection is used
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
nt <- 80

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

# Generator
AA <- make.generator(nx, ny, dx, dy, sdx, sdy, land=land)

# Data likelihood
n <- dim(AA)[1]
datlik <- matrix(0, nt, n)
for(i in 1:nt) datlik[i, ] <- as.vector(outer(dnorm(xx, Yx[i], osd), dnorm(yy, Yy[i], osd)))
i <- 3
image(xx, yy, matrix(datlik[i, ], nx, ny))

ini <- matrix(0, nx, ny)
ini[round(0.2*nx), round(0.5*ny)] <- 1
inivec <- as.vector(ini)
inivec <- Matrix(inivec, 1, n)
inivec[land] <- 0
if(sum(inivec)==0) stop('Not initialised correctly!')

# TMB
# Close diagonals (jump north-south)
Sns <- make.ns(n, nx, land)
# Far diagonals (jump east-west)
Sew <- make.ew(n, nx, land)
logDx <- log(0.5 * (sdx/dx)^2)  # East west move rate (diffusion)
logDy <- log(0.5 * (sdy/dy)^2)  # North south move rate (diffusion)

mu <- 3
F <- max(abs(AA))
#I <- Diagonal(n)
I <- Matrix(0, n, n)
diag(I) <- rep(1, n)
P <- AA/F + I
m <- ceiling(F*dt + 4*sqrt(F*dt) + 5) # Expression from Grassmann

#pvec <- Matrix(0, 1, n)
#pvec[1, round(n/2)] <- 1
pvec <- inivec

P <- as(P, 'dgTMatrix')
I <- as(I, 'dgTMatrix')
Sew <- as(Sew, 'dgTMatrix')
Sns <- as(Sns, 'dgTMatrix')
pvec <- as(pvec, 'matrix')


#datlik <- matrix(1, nt, n)
#datlik <- as(datlik, 'dgTMatrix')

## Fit model
obj <- MakeADFun(data=list(datlik=datlik, I=I, pvec=pvec, dt=dt, m=m, Sns=Sns, Sew=Sew),
                 parameters=list(logDx=logDx, logDy=logDy, u=1.0),
                 random=NULL,
                 DLL='simple'
                 )

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
rep <- sdreport(obj)
ests <- opt$par
sds <- sqrt(diag(rep$cov.fixed))

sdxest <- sqrt(2*exp(ests[1]))*dx
sdxll <- sqrt(2*exp(ests[1] - 2*sds[1]))*dx
sdxul <- sqrt(2*exp(ests[1] + 2*sds[1]))*dx
sdyest <- sqrt(2*exp(ests[2]))*dy
sdyll <- sqrt(2*exp(ests[2] - 2*sds[2]))*dy
sdyul <- sqrt(2*exp(ests[2] + 2*sds[2]))*dy


obj$report()$nt
Stmb <- obj$report()$G
phi <- obj$report()$phi
psi <- obj$report()$psi
pred <- obj$report()$pred
Ftmb <- obj$report()$F
Ptmb <- obj$report()$P

XX <- as.vector(matrix(xx, nx, ny))
YY <- as.vector(matrix(yy, nx, ny, byrow=TRUE))

xest <- numeric(nt)
yest <- numeric(nt)
for(t in 1:nt){
    xest[t] <- sum(phi[t, ] * XX)
    yest[t] <- sum(phi[t, ] * YY)
}

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

# Compare results
system.time(pout2 <- obj$report()$phi[2, ])
#system.time(obj$fn(obj$par))
system.time(pout2b <- pvec %*% uniformization(AA, dt))
system.time(pout2c <- pvec %*% expm(AA*dt))
sum(pout2 - pout2b)
sum(pout2 - pout2c)

mat <- matrix(pout2, nx, ny)

par(mfrow=c(2, 2))
image(xx, yy, ini)
image(xx, yy, mat)
