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


#' @name sim.shmm
#' @title Simulate movement data
#' @details Simulates data
#' @param inp List with parameters specified in the ini key and time step in the dt key
#' @param nobs Optional specification of the number of simulated observations.
#' @param dt Time step between observations.
#' @return A list containing the simulated data.
#' @examples
#' sim <- sim.shmm(inp)
#' @export
sim.shmm <- function(inp, nobs=40, dt=NULL){
    sim <- inp
    sim$true <- inp$ini
    if (is.null(dt)){
        sim <- set.default(inp, 'simdt', 1)
    } else {
        sim$simdt <- dt
    }
    
    if (inp$datatype == 'xy'){
        Xx <- cumsum(rnorm(nobs, 0, exp(inp$ini$logsdx)*sqrt(sim$simdt)))
        Xy <- cumsum(rnorm(nobs, 0, exp(inp$ini$logsdy)*sqrt(sim$simdt)))
        Yx <- Xx + rnorm(nobs, 0, inp$osd)
        Yy <- Xy + rnorm(nobs, 0, inp$osd)
        # Store simulation
        sim$true$X <- Xx
        sim$true$Y <- Xy
        sim$obs$X <- Yx
        sim$obs$Y <- Yy
        sim$obstime$xy <- seq(1, nobs, by=sim$simdt) * time.fac()
        #sim$obstime$Y <- seq(1, nobs, by=sim$simdt)
    }

    #sim <- check.inp(sim)
    
    return(sim)
}
