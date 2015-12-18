
make.generator <- function(nx, ny, dx, dy, sdx, sdy, land=NULL){
    Dx <- 0.5 * (sdx/dx)^2;  # Calculate diffusion parameter
    Dy <- 0.5 * (sdy/dy)^2;  # Calculate diffusion parameter

    n <- nx*ny
    # Skeleton matrix for generator
    S <- Matrix(0, n, n)
    doS <- Diagonal(n, 1)
    # Close diagonals (jump north-south)
    inds1 <- 1:(n-1)
    inds2 <- 2:n
    sub <- doS[inds2, inds2]
    S[inds1, inds2] <- Dx*sub # Upper
    S[inds2, inds1] <- S[inds2, inds1] + Dx*sub # Lower
    # Far diagonals (jump east-west)
    inds1 <- 1:(n-nx)
    inds2 <- (nx+1):n
    sub <- doS[inds2, inds2]
    S[inds1, inds2] <- S[inds1, inds2] + Dy*sub # Upper
    S[inds2, inds1] <- S[inds2, inds1] + Dy*sub # Lower
    # Take care of top of domain
    top <- seq(nx+1, n, by=nx)
    S[top, top-1] <- 0
    # Take care of bottom of domain
    bot <- seq(nx, n-nx, by=nx)
    S[bot, bot+1] <- 0
    # Remove land
    if(!is.null(land)){
        S[land, ] <- 0
        S[, land] <- 0
    }
    # Main diagonal
    diag(S) <- -apply(S, 1, sum)
    return(S)
}
