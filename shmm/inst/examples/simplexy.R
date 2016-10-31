rm(list=ls())
library(shmm)

set.seed(123)

inp <- list()
inp$datatypes <- 'xy'
inp$do.estimation <- TRUE
inp$do.sd.report <- FALSE
inp$dosmoo <- 0
inp$doviterbi <- 0

# Parameters
inp$ini$logsdx <- log(5)
inp$ini$logsdy <- log(6)
inp$osd <- 1 # Observation error sd
#inp$dt <- 1
inp$maxm <- 25

inp$grid$dx <- 1.2
inp$grid$dy <- 1.2

# Simulate data
inp <- sim.shmm(inp, nobs=10)

# Initial values
inp$ini$logsdx <- log(8)
inp$ini$logsdy <- log(8)

# Make grid
fracext <- 0.2
# X
xmin <- min(inp$true$X) - fracext*diff(range(inp$true$X))
xmax <- max(inp$true$X) + fracext*diff(range(inp$true$X))
inp$grid$xx <- seq(xmin, xmax, by=inp$grid$dx)
inp$grid$nx <- length(inp$grid$xx)
# Y
ymin <- min(inp$true$Y) - fracext*diff(range(inp$true$Y))
ymax <- max(inp$true$Y) + fracext*diff(range(inp$true$Y))
inp$grid$yy <- seq(ymin, ymax, by=inp$grid$dy)
inp$grid$ny <- length(inp$grid$yy)
# Land
indsy <- which(inp$grid$yy > 6 & inp$grid$yy < 8)
indsx <- which(inp$grid$xx > 1 & inp$grid$xx < 4)
dummy <- matrix(0, inp$grid$nx, inp$grid$ny)
inp$land <- dummy == 1

# Calculate data likelihood
inp$solvetype <- 'uniformisation'
inp <- calc.data.likelihood(inp)

cat('nx: ', inp$grid$nx, ' ny: ', inp$grid$ny, ' no states: ', inp$grid$nx*inp$grid$ny, ' ns: ', inp$ns, '\n')


# Fit shmm
res <- fit.shmm(inp)

plotshmm.distr(res, name='phi', sleep=0.05, add.map=FALSE)
