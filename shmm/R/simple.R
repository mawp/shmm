# Todo:
# DONE 1. Construct generator inside TMB
# DONE 2. Implement uniformization algorithm, input a vector that will get multiplied and output it afterward
# DONE 3. Make 2D version
# DONE 4. Implement filter
# 5. Implement simulator and data likelihood calculator for simple example
# 6. Implement likelihood function and estimate
# 7. Implement smoother
# 8. Exploit symmetry of G somehow? see http://eigen.tuxfamily.org/dox-devel/group__TutorialSparse.html
#    Symmetry is broken if advection is used
rm(list=ls())
require(TMB)
source('shmmfuns.R')
compile('../src/simple.cpp')
scriptname <- "../src/simple"
dyn.load(dynlib(scriptname))

# SHMM stuff
require(Matrix)
# X
sdx <- 2
rx <- 10
nx <- 11
xx <- seq(0, rx, length=nx)
dx <- xx[2]-xx[1]
# Y
sdy <- 1.5
ry <- 10
ny <- 21
yy <- seq(0, ry, length=ny)
dy <- yy[2]-yy[1]
# Land
indsy <- which(yy > 6 & yy < 8)
indsx <- which(xx > 1 & xx < 4)
dummy <- matrix(0, nx, ny)
dummy[indsx, indsy] <- 1
land <- which(dummy==1)
# Generator
AA <- make.generator(nx, ny, dx, dy, sdx, sdy, land=land)

n <- dim(AA)[1]
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

dt <- 1
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

nt <- 10
datlik <- matrix(1, nt, n)
#datlik <- as(datlik, 'dgTMatrix')

## Fit model
obj <- MakeADFun(data=list(datlik=datlik, I=I, pvec=pvec, mu=mu, dt=dt, m=m, logDx=logDx, logDy=logDy, Sns=Sns, Sew=Sew),
                 parameters=list(u=1.0),
                 random=NULL,
                 DLL='simple'
                 )

obj$report()$nt
Stmb <- obj$report()$G
phi <- obj$report()$phi
pred <- obj$report()$pred
Ftmb <- obj$report()$F
Ptmb <- obj$report()$P

t <- 3
matpred <- matrix(pred[t, ], nx, ny)
matphi <- matrix(phi[t, ], nx, ny)
par(mfrow=c(2, 2))
image(xx, yy, matpred)
image(xx, yy, matphi)

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
