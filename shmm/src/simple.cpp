// Normal linear mixed model specified through sparse design matrices.
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);          // Observations
  DATA_SPARSE_MATRIX(B);   // Random effect design matrix
  DATA_SPARSE_MATRIX(A);   // Fixed effect design matrix
  //DATA_SPARSE_MATRIX(S);   // Skeleton of generator matrix
  DATA_SPARSE_MATRIX(P);   // Sub-stochastic matrix
  DATA_SPARSE_MATRIX(I);   // Identity matrix
  DATA_SPARSE_MATRIX(pvec);// Vector to be multiplied by P
  DATA_SCALAR(mu);         // Movement parameter
  DATA_SCALAR(F);          // Numerical largest rate
  DATA_SCALAR(dt);         // Time step
  DATA_INTEGER(m);         // Number of iterations of uniformization
  PARAMETER_VECTOR(u);     // Random effects vector
  PARAMETER_VECTOR(beta);  // Fixed effects vector
  PARAMETER(logsdu);       // Random effect standard deviations
  PARAMETER(logsd0);       // Measurement standard deviation

  //int n = P.rows();

  // Uniformization

  //Eigen::SparseMatrix<Type> S = I;
  //Eigen::SparseMatrix<Type> pt = I;
  Eigen::SparseMatrix<Type> FPdt = F*P*dt;
  Eigen::SparseMatrix<Type> svec = pvec;
  Eigen::SparseMatrix<Type> pout2 = pvec;
  Type ind = 1.0;
  Type fact;
  for(int i=0; i<m; i++){
    fact = exp(lgamma(ind+1.0));    
    //fact = exp(lgamma(CppAD::<double>(i+1)));
    // Create transition matrix (actually not needed)
    //S = S * FPdt; // Matrix x Matrix
    //pt = pt + S/fact;
    // Update state vector
    svec = svec * FPdt; // Vector x Matrix
    pout2 = pout2 + svec/fact;
    // Update index
    ind += 1.0; // This should be changed to use i converted to Type
  }
  // Multiply by the constant
  //pt = pt * exp(-F*dt);
  //Eigen::SparseMatrix<Type> pout = pvec * pt;
  pout2 = pout2 * exp(-F*dt);


  // Distribution of random effect (u):
  Type ans=0;
  ans-=dnorm(u, Type(0) /*mean*/, exp(logsdu) /*sd*/, 1 /*log?*/).sum();

  // Distribution of obs given random effects (x|u):
  vector<Type> y=A*beta+B*u;
  ans-=dnorm(x,y,exp(logsd0),1).sum();

  //REPORT(P);
  //REPORT(pout);
  REPORT(pout2);
  //REPORT(S);
  //REPORT(fact);
  //REPORT(pt);

  return ans;
}

