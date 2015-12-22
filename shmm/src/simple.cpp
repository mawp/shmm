// Spatial hidden Markov model
// 21.12.2015
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(datlik);     // Data likelihood
  DATA_SPARSE_MATRIX(I);   // Identity matrix
  DATA_SCALAR(dt);         // Time step
  DATA_INTEGER(m);         // Number of iterations of uniformization
  DATA_SPARSE_MATRIX(Sns); // North-south generator skeleton
  DATA_SPARSE_MATRIX(Sew); // East-west generator skeleton
  DATA_VECTOR(lgam);       // Factorial
  PARAMETER(logDx);        // log diffusion in east-west (x) direction
  PARAMETER(logDy);        // log diffusion in north-south (y) direction

  int nt = datlik.rows();
  int n = datlik.cols();

  // Calculate components for uniformization
  Type Dx = exp(logDx);
  Type Dy = exp(logDy);
  Type F = 2 * (Dx + Dy); // Absolute largest jump rate, max(abs(diag(G)))
  Eigen::SparseMatrix<Type> G = Dx*Sew + Dy*Sns; // Make generator
  Eigen::SparseMatrix<Type> P = G/F + I; // Sub-stochastic matrix
  Eigen::SparseMatrix<Type> FPdt = F*P*dt;
  //Eigen::SparseMatrix<Type> FPdt = F*dt*((Dx*Sew + Dy*Sns)/F + I);

  // Initialise HMM grids
  matrix<Type> pred(nt, n);
  matrix<Type> phi(nt, n);
  vector<Type> psi(nt-1);
  pred.row(0) = datlik.row(0) / datlik.row(0).sum();
  phi.row(0) = pred.row(0);

  // Filter loop
  for(int t=1; t<nt; t++){
    // Time update using uniformization algorithm
    matrix<Type> svec = phi.row(t-1);
    matrix<Type> predtmp = svec;
    for(int i=0; i<m; i++){
      svec = svec * FPdt; // Vector x Matrix (this can be optimised?)
      //predtmp = predtmp + svec/exp(lgamma(Type(i+2))); // exp(lgamma()) is factorial
      predtmp = predtmp + svec/exp(lgam(i)); // exp(lgamma()) is factorial
    }
    predtmp = predtmp * exp(-F*dt);
    predtmp = predtmp / predtmp.sum(); // Ensure total probability mass is 1, should be a minor correction
    pred.row(t) = predtmp; // Store prediction
    
    // Data update
    matrix<Type> post = pred.row(t).cwiseProduct(datlik.row(t)); // Element-wise product
    psi(t-1) = post.sum(); 
    phi.row(t) = post / (psi(t-1) + 1e-20); // Add small value to avoid division by zero
  }

  // Negative log likelihood
  Type ans = -sum(log(psi));

  // Reports
  REPORT(pred);
  REPORT(phi);
  REPORT(psi);

  return ans;
}

