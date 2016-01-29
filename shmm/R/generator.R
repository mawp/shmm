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


#' @name make.generator
#' @title Create generator
#' @param nx Number of grid cells in x-direction
#' @param ny Number of grid cells in y-direction
#' @param Dx Diffusivity in x-direction
#' @param Dy Diffusivity in y-direction
#' @param land Vector containing linear indices of land cells.
#' @return Generator as sparse matrix.
#' @export
make.generator <- function(nx, ny, Dx, Dy, land=NULL){
    #Dx <- 0.5 * (sdx/dx)^2;  # East west move rate (diffusion)
    #Dy <- 0.5 * (sdy/dy)^2;  # North south move rate (diffusion)
    n <- nx*ny
    # Close diagonals (jump north-south)
    Sns <- make.ns(n, nx, land)
    # Far diagonals (jump east-west)
    Sew <- make.ew(n, nx, land)
    # Join NS and EW
    S <- Dx*Sew + Dy*Sns
    return(S)
}


#' @name make.ew
#' @title Make skeleton of generator matrix for east-west movement
#' @param n Total number of grid cells (nx*ny)
#' @param nx Number of grid cells in x-direction
#' @param land Vector containing linear indices of land cells.
#' @return Skeleton of generator matrix for east-west movement as sparse matrix
#' @export
make.ew <- function(n, nx, land){
    # Skeleton matrix for generator
    S <- Matrix::Matrix(0, n, n)
    doS <- Matrix::Diagonal(n, 1)
    # Close diagonals (jump north-south)
    inds1 <- 1:(n-1)
    inds2 <- 2:n
    sub <- doS[inds2, inds2]
    S[inds1, inds2] <- sub # Upper
    S[inds2, inds1] <- S[inds2, inds1] + sub # Lower
    # Remove top, bottom, land
    S <- remove.no.access(S, n, nx, land)
    # Main diagonal
    diag(S) <- -apply(S, 1, sum)
    return(S)
}


#' @name make.ns
#' @title Make skeleton of generator matrix for north-south movement
#' @param n Total number of grid cells (nx*ny)
#' @param nx Number of grid cells in x-direction
#' @param land Vector containing linear indices of land cells.
#' @return Skeleton of generator matrix for north-south movement as sparse matrix
#' @export
make.ns <- function(n, nx, land){
    # Skeleton matrix for generator
    S <- Matrix::Matrix(0, n, n)
    doS <- Matrix::Diagonal(n, 1)
    # Far diagonals (jump east-west)
    inds1 <- 1:(n-nx)
    inds2 <- (nx+1):n
    sub <- doS[inds2, inds2]
    S[inds1, inds2] <- S[inds1, inds2] + sub # Upper
    S[inds2, inds1] <- S[inds2, inds1] + sub # Lower
    # Remove top, bottom, land
    S <- remove.no.access(S, n, nx, land)
    # Main diagonal
    diag(S) <- -apply(S, 1, sum)
    return(S)
}


#' @name remove.no.access
#' @title Remove cells from generator that are inaccessible
#' @param S Generator matrix
#' @param n Total number of grid cells (nx*ny)
#' @param nx Number of grid cells in x-direction
#' @param land Vector containing linear indices of land cells.
#' @return Generator as sparse matrix with inaccessible cells removed.
#' @export
remove.no.access <- function(S, n, nx, land){
    # Take care of top of domain
    top <- seq(nx+1, n, by=nx)
    S[top, top-1] <- 0
    # Take care of bottom of domain
    bot <- seq(nx, n-nx, by=nx)
    S[bot, bot+1] <- 0
    # Remove land
    if(length(land) != 0){
        S[land, ] <- 0
        S[, land] <- 0
    }
    return(S)
}


#' @name uniformization
#' @title Calculate the matrix exponential of a generator matrix using the uniformization algorithm
#' @param A Generator matrix
#' @param dt Time step
#' @return Transition probability matrix as sparse matrix
#' @export
uniformization <- function(A, dt){
  ## Uniformization is an efficient way to compute the matrix exponential for a generator matrix
  ## See Grassmann 1977 - Transient solutions in Markovian queueing systems
  N <- dim(A)[1]
  # Find the numerical largest rate
  F <- max(abs(diag(A)))
  # Calculate number of iterations based expression in Grassmann 1977, eq 10
  m <- ceiling(F*dt + 4*sqrt(F*dt) + 5)
  # Insert warning if m>140 ??
  I <- Matrix::Diagonal(N)
  P <- A/F + I # Create sub-stochastic matrix (eq 8)
  S <- I
  pt <- I
  FPdt <- F*P*dt
  for(i in 1:m){
    S <- S %*% FPdt
    fact <- exp(lgamma(i+1))
    pt <- pt + S/fact
  }
  # Multiply by the constant
  pt <- pt * exp(-F*dt)
  pt
}
