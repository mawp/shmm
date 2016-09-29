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


#' @name add.gen
#' @title Add gen key to inp
#' @param inp Input list
#' @param ny Number of grid cells in y-direction
#' @param Dx Diffusivity in x-direction
#' @param Dy Diffusivity in y-direction
#' @param land Vector containing linear indices of land cells.
#' @return Generator as sparse matrix.
#' @export
add.gen <- function(inp){
    nx <- inp$grid$nx
    ny <- inp$grid$ny
    n <- nx*ny

    # Close diagonals (jump north-south)
    Sns <- make.ns(n, nx, inp$land)
    # Far diagonals (jump east-west)
    Sew <- make.ew(n, nx, inp$land)
    # Diffusivity
    # East west move rate (diffusion)
    inp$gen$Dx <- sd2D(exp(inp$ini$logsdx), inp$grid$dx)
    # North south move rate (diffusion)
    inp$gen$Dy <- sd2D(exp(inp$ini$logsdy), inp$grid$dy)

    # Calculate m (number of uniformization iterations) if not specified
    if (!'m' %in% names(inp$gen)){
        inp$gen$m <- calc.m(inp)
    }

    # Convert to dgTMatrix (required format by TMB)
    inp$gen$I <- as(diag.sparse(n), 'dgTMatrix')
    inp$gen$Sew <- as(Sew, 'dgTMatrix')
    inp$gen$Sns <- as(Sns, 'dgTMatrix')

    inp$gen$lgam <- lgamma(2:(inp$gen$m+2)) # Factorial
    
    return(inp)
}


#' @name diag.sparse
#' @title Create sparse identity matrix
#' @param n Dimension of matrix
#' @return Sparse identity matrix
diag.sparse <- function(n){
    I <- Matrix::Matrix(0, n, n)
    diag(I) <- rep(1, n)
    return(I)
}


#' @name calc.m
#' @title Calculate m (number of uniformization iterations)
#' @details
#' This uses the expression stated in Grassmann (1977) eq. 10.
#' Using this m guarantees that the truncation error of elements in a transition
#' probability matrix is less than 1e-4.
#' @param grid Spatial grid details
#' @param gen Generator details
#' @param land Specification of land cells
#' @return The number of uniformization iterations
calc.m <- function(inp){
    n <- inp$grid$nx * inp$grid$ny
    AA <- make.generator(inp$grid$nx, inp$grid$ny, inp$gen$Dx, inp$gen$Dy, inp$land)
    F <- max(abs(AA))
    m <- ceiling(F*inp$dt + 4*sqrt(F*inp$dt) + 5) # Expression from Grassmann
    return(m)
}


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
    S <- Matrix::bandSparse(n, k = 1, diag = list(rep(1, n-1)), symm=TRUE)
    #S <- Matrix::Matrix(0, n, n)
    #doS <- Matrix::Diagonal(n, 1)
    # Close diagonals (jump north-south)
    #inds1 <- 1:(n-1)
    #inds2 <- 2:n
    #sub <- doS[inds2, inds2]
    #S[inds1, inds2] <- sub # Upper
    #S[inds2, inds1] <- S[inds2, inds1] + sub # Lower
    # Remove top, bottom, land
    S <- remove.no.access(S, n, nx, land)
    # Main diagonal
    #S <- S - Matrix::Diagonal(x=Matrix::rowSums(S))
    diag(S) <- -Matrix::rowSums(S)
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
    S <- Matrix::bandSparse(n, k = nx, diag = list(rep(1, n-nx)), symm=TRUE)
    #S <- Matrix::Matrix(0, n, n)
    #doS <- Matrix::Diagonal(n, 1)
    # Far diagonals (jump east-west)
    #inds1 <- 1:(n-nx)
    #inds2 <- (nx+1):n
    #sub <- doS[inds2, inds2]
    #S[inds1, inds2] <- S[inds1, inds2] + sub # Upper
    #S[inds2, inds1] <- S[inds2, inds1] + sub # Lower
    # Remove top, bottom, land
    S <- remove.no.access(S, n, nx, land)
    # Main diagonal
    #S <- S - Matrix::Diagonal(x=Matrix::rowSums(S))
    diag(S) <- -Matrix::rowSums(S)
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
