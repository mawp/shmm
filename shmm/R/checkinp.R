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


#' @name check.inp
#' @title Check list of input variables
#' @details Check list of input variables
#' @param inp List of input variables, see details for required variables.
#' @return An updated list of input variables checked for consistency and with defaults added.
#' @examples
#' inp <- check.inp(inp)
#' @export
check.inp <- function(inp){
    set.default <- function(inpin, key, val){
        if (!key %in% names(inpin)){
            inpin[[key]] <- val
        }
        return(inpin)
    }

    # Check grid
    if (!'grid' %in% names(inp)){
        stop('grid not defined!')
    } else {
        if (!'xx' %in% names(inp$grid)){
            stop('xx not specified in inp$grid!')
        }
        if (!'yy' %in% names(inp$grid)){
            stop('yy not specified in inp$grid!')
        }
        if (!'nx' %in% names(inp$grid)){
            stop('nx not specified in inp$grid!')
        }
        if (!'ny' %in% names(inp$grid)){
            stop('ny not specified in inp$grid!')
        }
        if (!'dx' %in% names(inp$grid)){
            stop('dx not specified in inp$grid!')
        }
        if (!'dy' %in% names(inp$grid)){
            stop('dy not specified in inp$grid!')
        }
        inp$grid$n <- inp$grid$nx * inp$grid$ny
    }

    
    inp$scriptname <- 'shmm'
    # Set defaults that were not manually defined
    inp <- set.default(inp, 'dt', 1)
    inp <- set.default(inp, 'do.sd.report', TRUE)
    inp <- set.default(inp, 'dosmoo', 1)

    # Calculate time vector
    inp$timerange <- range(unlist(inp$time))
    inp$times <- seq(inp$timerange[1], inp$timerange[2], by=inp$dt)
    inp$ns <- length(inp$times)
    
    return(inp)
}
