// Normal linear mixed model specified through sparse design matrices.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //DATA_SPARSE_MATRIX(datlik);     // 
  DATA_MATRIX(datlik);     // 
  //DATA_SPARSE_MATRIX(Pin);   // Sub-stochastic matrix
  DATA_SPARSE_MATRIX(I);   // Identity matrix
  //DATA_SPARSE_MATRIX(pvec);// Vector to be multiplied by P
  DATA_MATRIX(pvec);// Vector to be multiplied by P
  //DATA_SCALAR(mu);         // Movement parameter
  //DATA_SCALAR(Fin);          // Numerical largest rate
  DATA_SCALAR(dt);         // Time step
  DATA_INTEGER(m);         // Number of iterations of uniformization
  DATA_SPARSE_MATRIX(Sns); // 
  DATA_SPARSE_MATRIX(Sew); // 
  //DATA_SCALAR(logDx);      // 
  //DATA_SCALAR(logDy);      // 
  PARAMETER(logDx);      // 
  PARAMETER(logDy);      // 
  PARAMETER(u);            // Random effects vector

  // Do dummy optimisation
  Type ans=0;
  ans -= dnorm(u, Type(0.0), Type(1.0), 1);

  // === SHMM stuff ===
  Type Dx = exp(logDx);
  Type Dy = exp(logDy);
  Eigen::SparseMatrix<Type> G = Dx*Sew + Dy*Sns; // Make generator
  Type F = 2 * (Dx + Dy); // Absolute largest jump rate, max(abs(diag(G)))
  Eigen::SparseMatrix<Type> P = G/F + I; // Sub-stochastic matrix

  int nt = datlik.rows();
  int n = datlik.cols();
  matrix<Type> pred(nt, n);
  matrix<Type> phi(nt, n);
  vector<Type> psi(nt-1);
  //pred.row(0) += pred.row(0); // Doesn't work with sparse, works with dense
  //pred.col(0) += pred.col(0); // Works with both sparse and dense
  pred.row(0) = pvec; // Doesn't work with sparse, works with dense
  phi.row(0) = pvec; // Doesn't work with sparse, works with dense
  Eigen::SparseMatrix<Type> FPdt = F*P*dt;
  // Filter loop
  for(int t=1; t<nt; t++){
    // Time update

    // --- Uniformization begin --- 
    matrix<Type> svec = phi.row(t-1);
    matrix<Type> pout2 = svec;
    Type ind = 1.0;
    for(int i=0; i<m; i++){
      Type fact = exp(lgamma(ind+1.0));    
      //fact = exp(lgamma(CppAD::<double>(i+1)));
      // Update state vector
      svec = svec * FPdt; // Vector x Matrix
      pout2 = pout2 + svec/fact;
      // Update index
      ind += 1.0; // This should be changed to use i converted to Type
    }
    pout2 = pout2 * exp(-F*dt);
    // --- uniformization end ---

    pout2 = pout2 / pout2.sum(); // Ensure total probability mass is 1, should be a minor correction
    pred.row(t) = pout2; // Store prediction
    
    // Data update
    matrix<Type> post = pred.row(t).cwiseProduct(datlik.row(t)); // Element-wise product
    psi(t-1) = post.sum(); // Add small value to avoid division by zero
    phi.row(t) = post / (psi(t-1) + 1e-20);

  }

  ans -= sum(log(psi));

  REPORT(nt);
  REPORT(pred);
  REPORT(phi);
  REPORT(psi);
  REPORT(G);
  REPORT(F);
  REPORT(P);

  return ans;
}

