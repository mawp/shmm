# Spatial hidden Markov model (SHMM)
#    Copyright (C) 2015-2016  Martin Waever Pedersen, mawp@dtu.dk or wpsgodd@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' @name fit.shmm
#' @title Fit a continuous-time surplus production model to data.
#' @details Fit a continuous-time surplus production model to data.
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return List containing results.
#' @export
#' @examples
#' rep <- fit.shmm(inp)
#' @import TMB
fit.shmm <- function(inp, dbg=0){
    rep <- NULL
    # Check input list
    inp <- check.inp(inp)
    nt <- dim(inp$datlik)[1]
    nx <- inp$grid$nx
    ny <- inp$grid$ny
    n <- nx*ny

    datlikin <- matrix(0, nt, n)
    for (i in 1:nt) datlikin[i, ] <- as.vector(inp$datlik[i, , ])
    # TMB
    # Close diagonals (jump north-south)
    Sns <- make.ns(n, nx, inp$land)
    # Far diagonals (jump east-west)
    Sew <- make.ew(n, nx, inp$land)
    Dx <- 0.5 * (inp$ini$sdx/inp$grid$dx)^2  # East west move rate (diffusion)
    Dy <- 0.5 * (inp$ini$sdy/inp$grid$dy)^2  # North south move rate (diffusion)

    # Generator
    #AA <- make.generator(nx, ny, dx, dy, sdx, sdy, land=land)
    AA <- make.generator(nx, ny, Dx, Dy, land=inp$land)
    F <- max(abs(AA))
    I <- Matrix::Matrix(0, n, n)
    diag(I) <- rep(1, n)
    #P <- AA/F + I
    m <- ceiling(F*inp$dt + 4*sqrt(F*inp$dt) + 5) # Expression from Grassmann
    #G <- uniformization(AA, dt)

    #P <- as(P, 'dgTMatrix')
    I <- as(I, 'dgTMatrix')
    Sew <- as(Sew, 'dgTMatrix')
    Sns <- as(Sns, 'dgTMatrix')

    # Create objective function
    lgam <- lgamma(2:(m+2))
    data <- list(datlik=datlikin, I=I, dt=inp$dt, m=m, Sns=Sns, Sew=Sew, lgam=lgam)
    pars <- list(logDx=log(Dx), logDy=log(Dy))
    obj <- TMB::MakeADFun(data=data, parameters=pars, random=NULL, DLL=inp$scriptname) # This one requires memory

    # Estimate
    system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
    system.time(rep <- TMB::sdreport(obj))

    if(!is.null(rep)) class(rep) <- "shmmcls"
    return(rep)
}


#' @name make.datin
#' @title Create data list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. 
#' @return List to be used as data input to TMB::MakeADFun.
#' @export
make.datin <- function(inp, dbg=0){
    datin <- list()
    return(datin)
}


#' @name make.obj
#' @title Create TMB obj using TMB::MakeADFun and squelch screen printing.
#' @param datin Data list.
#' @param pl Parameter list.
#' @param inp List of input variables as output by check.inp.
#' @param phase Estimation phase, integer.
#' @return List to be used as data input to TMB.
#' @export
#' @import TMB
make.obj <- function(datin, pl, inp, phase=1){
    obj <- TMB::MakeADFun(data=datin, parameters=pl, random=inp$RE, DLL=inp$scriptname, hessian=TRUE, map=inp$map[[phase]])
    TMB:::config(trace.optimize=0, DLL=inp$scriptname)
    verbose <- FALSE
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    obj$fn(obj$par)
    return(obj)
}
