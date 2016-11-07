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


#' @name plotshmm.distr
#' @title Plots the distribution indicated by name
#' @details Distribution is plotted as an "animation"
#' @param rep Result of running fit.shmm.
#' @param name Name of the distribution to plot. Either 'phi', 'pred' or 'smoo'.
#' @param sleep Number of seconds to pause between each frame of the animation.
#' @param add.map Add map?
#' @return Noting.
#' @export
plotshmm.distr <- function(rep, name='smoo', sleep=0.05, add.map=TRUE){
    distr <- rep$report[[name]]
    nt <- dim(distr)[1]
    for (i in 1:nt){
        main <- paste('i =', i, ' Date:', inp$date[i])
        image(rep$inp$grid$xx, rep$inp$grid$yy, -distr[i, , ], xlab='X',
              ylab='Y', main=main)
        if ('tracks' %in% names(rep)){
            if ('mean' %in% names(rep$tracks)){
                lines(rep$tracks$mean$X[1:i], rep$tracks$mean$Y[1:i], col='blue')
                points(rep$tracks$mean$X[1:i], rep$tracks$mean$Y[1:i], col='blue',
                       pch=20, cex=0.7)
            }
            if ('viterbi' %in% names(rep$tracks)){
                lines(rep$tracks$viterbi$X[1:i], rep$tracks$viterbi$Y[1:i], col='red')
                points(rep$tracks$viterbi$X[1:i], rep$tracks$viterbi$Y[1:i], col='red',
                       pch=20, cex=0.7)
            }
        }
        if ('true' %in% names(rep$inp)){
            lines(rep$inp$true$X[1:i], rep$inp$true$Y[1:i], col='green')
        }
        if (add.map){
            require(maps)
            maps::map(add=TRUE)
        }
        if ('true' %in% names(rep$inp)){
            legend('topright', legend=c('True', 'Mean'), lty=1,
                   col=c('green', 'blue'), bg='white')
        } else {
            legend('topright', legend='Mean', lty=1, col='blue', bg='white')
        }
        Sys.sleep(sleep)
    }
}
