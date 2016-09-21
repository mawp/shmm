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
#' @return A list containing the simulated data.
#' @examples
#' sim <- sim.shmm(inp)
#' @export
sim.shmm <- function(inp, nobs=40){
    sim <- inp
    sim$true <- inp$ini
    
    if (inp$datatype == 'xy'){
        Xx <- cumsum(rnorm(nobs, 0, exp(inp$ini$logsdx)*sqrt(inp$dt)))
        Xy <- cumsum(rnorm(nobs, 0, exp(inp$ini$logsdy)*sqrt(inp$dt)))
        Yx <- Xx + rnorm(nobs, 0, inp$osd)
        Yy <- Xy + rnorm(nobs, 0, inp$osd)
        # Store simulation
        sim$true$X <- Xx
        sim$true$Y <- Xy
        sim$obs$X <- Yx
        sim$obs$Y <- Yy
        sim$time$X <- 1:nobs
        sim$time$Y <- 1:nobs
    }

    #sim <- check.inp(sim)
    
    return(sim)
}
