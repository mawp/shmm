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
#' @param dosmoo Do smoothing?
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return List containing results.
#' @export
#' @examples
#' rep <- fit.shmm(inp)
#' @import TMB
fit.shmm <- function(inp, dosmoo=1, dbg=0){
    # Check input list
    inp <- check.inp(inp)
    rep <- list(inp=inp)
    
    tic <- Sys.time()
    cat('\nCreate objective function...')
    obj <- make.obj(inp, dosmoo=dosmoo) # Smoothing not performed here

    if (inp$do.estimation){
        cat('\nEstimating parameters...')
        rep$opt <- try(nlminb(obj$par, obj$fn, obj$gr))
    } else {
        cat('\nEvaluating for fixed parameters...')
        nouse <- obj$fn() # Evaluate fn for initial parameters
    }
    rep$dosmoo <- dosmoo
    
    cat('\nSave report...')
    rep$report <- obj$report()
    for (nm in names(rep$report)){
        if (nm != 'psi'){
            nmpl <- sub('out', '', nm)
            vecdistr <- rep$report[[nm]]
            if (length(dim(vecdistr)) == 2){
                rep$report[[nmpl]] <- format.distr(vecdistr, rep)
                #rep$report[[nm]] <- NULL
            }
        }
    }

    #cat('\nSmoothing...')
    obj2 <- obj
    obj2$env$data$dosmoo <- 1 # Smooth now
    #nouse <- obj2$fn() # Evaluate fn to do smoothing
    #inp2 <- inp
    #obj2 <- make.obj(inp2) # Smoothing not performed here
    #obj2$env$data$dosmoo <- 1 # Smooth now
    nouse <- obj2$fn(obj2$env$last.par.best) # Evaluate fn to do smoothing

    
    cat('\nSave report2...')
    rep$report2 <- obj2$report()
    for (nm in names(rep$report2)){
        if (nm != 'psi'){
            nmpl <- sub('out', '', nm)
            vecdistr <- rep$report2[[nm]]
            if (length(dim(vecdistr)) == 2){
                rep$report2[[nmpl]] <- format.distr(vecdistr, rep)
                #rep$report2[[nm]] <- NULL
            }
        }
    }
    
    if (inp$do.sd.report){
        cat('\nCalculating uncertainties...')
        rep <- try(TMB::sdreport(obj))
    }

    rep$obj <- obj
    cat('\nCalculating track...')
    rep <- get.mean.track(rep)

    rep$computing.time <- difftime(Sys.time(), tic, unit='secs')
    
    if(!is.null(rep)){
        class(rep) <- "shmmcls"
    }
    cat('\n')
    return(rep)
}


#' @name make.obj
#' @title Create TMB obj using TMB::MakeADFun and squelch screen printing.
#' @param inp List of input variables as output by check.inp.
#' @param phase Estimation phase, integer.
#' @param dosmoo Do smoothing?
#' @return List to be used as data input to TMB.
#' @export
#' @import TMB
make.obj <- function(inp, phase=1, dosmoo=1){
    datlist <- make.datlist(inp, dosmoo=dosmoo)
    parlist <- make.parlist(inp, dosmoo=dosmoo)
    # This one requires memory if not using atomic
    obj <- TMB::MakeADFun(data=datlist,
                          parameters=parlist,
                          map=list(dosmoo=factor(NA)),
                          random=NULL,
                          DLL=inp$scriptname,
                          checkParameterOrder=FALSE)
    # Make TMB quiet
    TMB:::config(trace.optimize=0, DLL=inp$scriptname)
    verbose <- FALSE
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    obj$fn(obj$par)
    return(obj)
}


#' @name make.datlist
#' @title Create data list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @param dosmoo Do smoothing?
#' @return List to be used as data input to TMB::MakeADFun.
#' @export
make.datlist <- function(inp, dosmoo=1){
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
    return(datlist)
}


#' @name make.parlist
#' @title Create parameter list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @param dosmoo Do smoothing?
#' @return List to be used as parameter input to TMB::MakeADFun.
#' @export
make.parlist <- function(inp, dosmoo=1){
    parlist <- list(logDx=log(inp$gen$Dx),
                    logDy=log(inp$gen$Dy),
                    dosmoo=dosmoo)
    return(parlist)
}

