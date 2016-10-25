// Spatial hidden Markov model
// 21.12.2015
#include <TMB.hpp>

/* Class to build generator and project state one step forward  */
//struct {
  // Data for atomic function (note: double types)
  // Initialize at first evaluation with double types
namespace shmm {
  // Input that is constant (do not depend on parameters)
  int m;
  Eigen::SparseMatrix<double> I;
  Eigen::SparseMatrix<double> Sns;
  Eigen::SparseMatrix<double> Sew;
  double dt;
  vector<double> lgam;

  template<class Type>
  struct shmm_parms {
    // Input that is *not* constant:
    matrix<Type> svec;
    Type Dx;
    Type Dy;
    // Method: vector -> shmm_parms
    void operator=(vector<Type> x){
      int n = x.size() - 2;
      svec.resize(1, n);
      svec << x.head(n).transpose();
      Dx = x(n);
      Dy = x(n+1);
    }
    // Method: shmm_parms -> vector
    operator vector<Type>(){
      int n = svec.size();
      vector<Type> x(n + 2);
      x.head(n) << svec.vec();
      x(n) = Dx;
      x(n+1) = Dy;
      return x;
    }
  };


  // Build generator and project state one step forward
  template<class Type>
  matrix<Type> forwardProject(matrix<Type> svec, Type Dx, Type Dy){
    // Cast required 'double' objects to 'Type'
    Eigen::SparseMatrix<Type>   I = shmm::  I.cast<Type>();
    Eigen::SparseMatrix<Type> Sns = shmm::Sns.cast<Type>();
    Eigen::SparseMatrix<Type> Sew = shmm::Sew.cast<Type>();
    Type                       dt = shmm::dt;
    vector<Type>             lgam = shmm::lgam.cast<Type>();

    // Build generator
    Type F = 2 * (Dx + Dy); // Absolute largest jump rate, max(abs(diag(G)))
    Eigen::SparseMatrix<Type> G = Dx*Sew + Dy*Sns; // Make generator
    Eigen::SparseMatrix<Type> P = G/F + I; // Sub-stochastic matrix
    Eigen::SparseMatrix<Type> FPdt = F*P*Type(dt);

    // One-step forward
    matrix<Type> predtmp = svec;
    for(int i=0; i<m; i++){
      svec = svec * FPdt; // Vector x Matrix (this can be optimised?)
      predtmp = predtmp + svec/exp(lgam(i)); // exp(lgamma()) is factorial
    }
    predtmp = predtmp * exp(-F*dt);
    predtmp = predtmp / predtmp.sum(); // Ensure total probability mass is 1, should be a minor correction
    return predtmp;
  }
//};

  // Wrapper: vector input -> vector output
  template<class Type>
  vector<Type> forwardProject(vector<Type> input){
    shmm_parms<Type> parms;
    parms = input;
    vector<Type> output(parms.svec.size());
    output << forwardProject(parms.svec, parms.Dx, parms.Dy).vec();
    return output;
  }
  REGISTER_ATOMIC(forwardProject)

  // User version
  template<class Type>
  matrix<Type> ForwardProject(matrix<Type> svec, Type Dx, Type Dy){
    shmm_parms<Type> parms = {svec, Dx, Dy};
    vector<Type> x = parms;
    vector<Type> y = forwardProject(x);
    matrix<Type> out(1, y.size());
    out << y.transpose();
    return out;
  }
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(datlik);     // Data likelihood
  DATA_INTEGER(dosmoo);    // If 1 smoothing is done
  DATA_INTEGER(ns);        // Number of time steps of solver
  DATA_VECTOR(iobs);       // Indices to which observations correspond
  PARAMETER(logDx);        // log diffusion in east-west (x) direction
  PARAMETER(logDy);        // log diffusion in north-south (y) direction

  // Transfer all constant data to namespace 'shmm'
  if(isDouble<Type>::value){
#define Type double
    DATA_SPARSE_MATRIX(I);   // Identity matrix
    DATA_SCALAR(dt);         // Time step
    DATA_INTEGER(m);         // Number of iterations of uniformization
    DATA_SPARSE_MATRIX(Sns); // North-south generator skeleton
    DATA_SPARSE_MATRIX(Sew); // East-west generator skeleton
#undef Type
    shmm::m = m;
    shmm::I = I;
    shmm::Sns = Sns;
    shmm::Sew = Sew;
    shmm::dt = dt;
      // Dirty trick didn't work for DATA_VECTOR:
    // DATA_VECTOR(lgam);       // Factorial
    shmm::lgam = asVector<double>(getListElement(objective_function::data,"lgam",&isNumeric));
  }

  int nobs = datlik.rows();
  int n = datlik.cols();

  // // Calculate components for uniformization
  Type Dx = exp(logDx);
  Type Dy = exp(logDy);

  // Initialise HMM grids
  matrix<Type> pred(ns, n);
  matrix<Type> phi(ns, n);
  vector<Type> psi(nobs - 1);
  // First state is at time of first observation
  pred.row(0) = datlik.row(0) / datlik.row(0).sum();
  phi.row(0) = pred.row(0);

  // Filter loop
  for(int t=1; t<ns; t++){
    // Time update using uniformization algorithm
    matrix<Type> svec = phi.row(t-1);
    matrix<Type> predtmp = shmm::ForwardProject(svec, Dx, Dy);
    pred.row(t) = predtmp; // Store prediction
    
    if (iobs(t) > 0){
      // Data update
      int ind = CppAD::Integer(iobs(t)-1);
      matrix<Type> post = pred.row(t).cwiseProduct(datlik.row(ind)); // Element-wise product
      psi(ind - 1) = post.sum(); 
      phi.row(t) = post / (psi(ind - 1) + 1e-20); // Add small value to avoid division by zero
    } else {
      // No data update
      phi.row(t) = pred.row(t);
    }
  }

  // Negative log likelihood
  Type ans = -sum(log(psi));

  // Smoothing
  // TODO: only run smoothing once after estimation is completed
  matrix<Type> smoo(ns, n);
  if (dosmoo == 1){
    // Smoothing loop
    smoo.row(ns-1) = phi.row(ns-1);
    for(int t=1; t < ns; t++){
      int tt = ns - t;
      // Time update using uniformization algorithm
      matrix<Type> predrow = pred.row(tt);
      for (int i=0; i < n; i++){
	predrow(0, i) += 1e-10;
      }
      matrix<Type> ratio = smoo.row(tt).cwiseQuotient(predrow);
      //matrix<double> asd(1, n);
      //for (int i=0; i < n; i++){
      //asd(0, i) = std::isnan(asDouble(ratio(0, i)));
	//cout 
      //}
      //matrix<double> asd = std::isnan(asDouble(ratio));
      matrix<Type> ratiotmp = shmm::ForwardProject(ratio, Dx, Dy);
      matrix<Type> post = phi.row(tt-1).cwiseProduct(ratiotmp);
      post = post / (post.sum() + 1e-20);
      smoo.row(tt-1) = post;
    }
  }

  // Store a subset of distribution for output
  matrix<Type> smooout(nobs, n);
  matrix<Type> phiout(nobs, n);
  matrix<Type> predout(nobs, n);
  for(int t=0; t<ns; t++){
    if (iobs(t) > 0){
      int ind = CppAD::Integer(iobs(t)-1);
      smooout.row(ind) = smoo.row(t);
      phiout.row(ind) = phi.row(t);
      predout.row(ind) = pred.row(t);
    }
  }

  // Reports
  REPORT(predout);
  REPORT(phiout);
  REPORT(psi);
  REPORT(smooout);

  return ans;
}

