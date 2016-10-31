rm(list=ls())
library(shmm)

set.seed(123)

inp <- list()
inp$datatypes <- 'xy'
inp$do.estimation <- TRUE
inp$do.sd.report <- FALSE

# Parameters used for simulation
inp$ini$logsdx <- log(5)
inp$ini$logsdy <- log(6)
inp$osd <- 1 # Observation error sd
#inp$dt <- 1
inp$maxm <- 25

inp$grid$dx <- 1.2 / 3
inp$grid$dy <- 1.2 / 3

# Simulate
inp <- sim.shmm(inp, nobs=10)

# Initial values
inp$ini$logsdx <- log(8)
inp$ini$logsdy <- log(8)

# Make grid
# X
xmin <- min(inp$true$X) - 0.1*diff(range(inp$true$X))
xmax <- max(inp$true$X) + 0.1*diff(range(inp$true$X))
inp$grid$xx <- seq(xmin, xmax, by=inp$grid$dx)
inp$grid$nx <- length(inp$grid$xx)
# Y
ymin <- min(inp$true$Y) - 0.1*diff(range(inp$true$Y))
ymax <- max(inp$true$Y) + 0.1*diff(range(inp$true$Y))
inp$grid$yy <- seq(ymin, ymax, by=inp$grid$dy)
inp$grid$ny <- length(inp$grid$yy)
# Land
indsy <- which(inp$grid$yy > 6 & inp$grid$yy < 8)
indsx <- which(inp$grid$xx > 1 & inp$grid$xx < 4)
dummy <- matrix(0, inp$grid$nx, inp$grid$ny)
#dummy[indsx, indsy] <- 1
inp$land <- dummy == 1
cat('nx: ', inp$grid$nx, ' ny: ', inp$grid$ny, ' no states: ', inp$grid$nx*inp$grid$ny, '\n')

# Calculate data likelihood
                                        # inp$solvetype <- 'uniformisation'
inp$solvetype <- 'implicit'
inp <- calc.data.likelihood(inp)
print(inp$dt)
print(length(inp$time))

# Setup object
datlist <- list(datlik=inp$datlik$all,
                solvetype=inp$solvetypein,
                ns=inp$ns,
                iobs=inp$iobs,
                I=inp$gen$I,
                dt=inp$dt,
                m=inp$gen$m,
                Sns=inp$gen$Sns,
                Sew=inp$gen$Sew,
                lgam=inp$gen$lgam)

parlist <- list(logDx=log(inp$gen$Dx),
                logDy=log(inp$gen$Dy),
                dosmoo=1)

compile('tinyad.cpp')
dyn.load(dynlib("tinyad"))
obj <- TMB::MakeADFun(data=datlist,
                      parameters=parlist,
                      random=NULL,
                      map=list(dosmoo=factor(NA)),
                      DLL='tinyad',
                      checkParameterOrder=FALSE)
# Make TMB quiet
TMB:::config(trace.optimize=0, DLL=inp$scriptname)
verbose <- FALSE
obj$env$tracemgc <- verbose
obj$env$inner.control$trace <- verbose
obj$env$silent <- ! verbose
obj$fn(obj$par)
obj$gr(obj$par)




