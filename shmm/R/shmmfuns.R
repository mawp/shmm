
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

make.ns <- function(n, nx, land){
    # Skeleton matrix for generator
    S <- Matrix(0, n, n)
    doS <- Diagonal(n, 1)
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

make.ew <- function(n, nx, land){
    # Skeleton matrix for generator
    S <- Matrix(0, n, n)
    doS <- Diagonal(n, 1)
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

make.generator <- function(nx, ny, dx, dy, sdx, sdy, land=NULL){
    Dx <- 0.5 * (sdx/dx)^2;  # East west move rate (diffusion)
    Dy <- 0.5 * (sdy/dy)^2;  # North south move rate (diffusion)

    n <- nx*ny
    # Close diagonals (jump north-south)
    Sns <- make.ns(n, nx, land)
    # Far diagonals (jump east-west)
    Sew <- make.ew(n, nx, land)
    # Join NS and EW
    S <- Dx*Sew + Dy*Sns
    return(S)
}


# --- Uniformization in R ---
uniformization <- function(A,dt){
  ## Uniformization is an efficient way to compute the matrix exponential for a generator matrix
  ## See Grassmann 1977 - Transient solutions in Markovian queueing systems
  N <- dim(A)[1]
  # Find the numerical largest rate
  F <- max(abs(diag(A)))
  # Calculate number of iterations based expression in Grassmann 1977, eq 10
  m <- ceiling(F*dt + 4*sqrt(F*dt) + 5)
  # Insert warning if m>140 ??
  I <- Diagonal(N)
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




