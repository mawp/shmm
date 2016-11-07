// Spatial hidden Markov model
// 21.12.2015
#include <TMB.hpp>

/* Class to build generator and project state one step forward  */
namespace shmm {
  // Input that is constant (do not depend on parameters)
  int m;
  Eigen::SparseMatrix<double> I;
  Eigen::SparseMatrix<double> Sns;
  Eigen::SparseMatrix<double> Sew;
  double dt;
  vector<double> lgam;
  matrix<double> datlik;     // Data likelihood
  int solvetype;             // Type of solver to use
  int ns;                    // Number of time steps of solver
  vector<int> iobs;          // Indices to which observations correspond

  /* 'filter' class calculates the entire likelihood (not just one
     step). We call the type 'Float' to emphasize that this is a
     tiny_ad type */
  template<class Float>
  struct filter_t {
    /* Objects we only calculate once during initialization */
    //    Eigen::SparseMatrix<Float> G;
    Eigen::SparseMatrix<Float> FPdt; // For uniformization
    Float F;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<Float> > solver; // For implicit

    // HMM grids (used by both 'eval' and 'smoothing')
    matrix<Float> pred;
    matrix<Float> phi;
    vector<Float> psi;
    matrix<Float> smoo; // output from smoothing

    /* Initialize the generator */
    void initialize(Float Dx, Float Dy) {
      // Cast required 'double' objects to 'Float'
      Eigen::SparseMatrix<Float>   I = shmm::  I.cast<Float>();
      Eigen::SparseMatrix<Float> Sns = shmm::Sns.cast<Float>();
      Eigen::SparseMatrix<Float> Sew = shmm::Sew.cast<Float>();
      Float                       dt = shmm::dt;
      vector<Float>             lgam = shmm::lgam.cast<Float>();
      // Build generator
      Eigen::SparseMatrix<Float> G = Dx*Sew + Dy*Sns;  // Make generator
      if (solvetype == 1) { // Uniformisation
	F = 2 * (Dx + Dy); // Absolute largest jump rate, max(abs(diag(G)))
	Eigen::SparseMatrix<Float> P = G/F + I; // Sub-stochastic matrix
	FPdt = F * P * Float(dt);
      }
      else if (solvetype == 2) { // Implicit solving
	Eigen::SparseMatrix<Float> A = I - G*dt;
	solver.analyzePattern(A);
	solver.factorize(A);
      }
      else error("'solvetype' can by 1 or 2");
    }

    // 1. Uniformization: project state one step forward
    matrix<Float> forwardProject(matrix<Float> svec){
      matrix<Float> predtmp = svec; // Initialise
      for(int i=0; i<m; i++){
	svec = svec * FPdt; // Vector x Matrix (this can be optimised?)
	predtmp = predtmp + svec/exp(lgam(i)); // exp(lgamma()) is factorial
      }
      predtmp = predtmp * exp(-F*dt);
      predtmp = predtmp / predtmp.sum(); // Ensure total probability mass is 1, should be a minor correction
      return predtmp;
    }

    // 2. Implicit: project state one step forward
    matrix<Float> forwardProjecti(matrix<Float> svec){
      matrix<Float> predtmp = solver.solve(svec.transpose());
      predtmp = predtmp / predtmp.sum(); // Ensure total probability mass is 1, should be a minor correction
      return predtmp.transpose();
    }

    /* Run filter loop */
    Float eval(Float Dx, Float Dy, bool do_smoothing = false) {
      matrix<Float> datlik = shmm::datlik.cast<Float>();
      int nobs = datlik.rows();
      int n = datlik.cols();
      // Initialize solver or FPdt
      initialize(Dx, Dy); // FIXME: Put in constructor (?)
      // Initialise HMM grids
      pred.resize(ns, n);
      phi.resize(ns, n);
      psi.resize(nobs - 1);
      // First state is at time of first observation
      pred.row(0) = datlik.row(0) / datlik.row(0).sum();
      phi.row(0) = pred.row(0);
      // Filter loop
      for(int t=1; t<ns; t++) {
      	// Time update using uniformization algorithm
      	matrix<Float> svec = phi.row(t-1);
      	matrix<Float> predtmp = svec;
        if (solvetype == 1) {
      	  predtmp = forwardProject(svec);
        }
        else if (solvetype == 2) {
      	  predtmp = forwardProjecti(svec);
        }
        else error("'solvetype' can by 1 or 2");
      	pred.row(t) = predtmp; // Store prediction
        if (iobs(t) > 0) {
      	  // Data update
      	  int ind = (iobs(t)-1);
      	  matrix<Float> post = pred.row(t).cwiseProduct(datlik.row(ind)); // Element-wise product
      	  psi(ind - 1) = post.sum(); 
      	  phi.row(t) = post / (psi(ind - 1) + 1e-20); // Add small value to avoid division by zero
      	} else {
	  // No data update
      	  phi.row(t) = pred.row(t);
      	}
      }
      // Negative log likelihood
      Float ans = -sum(log(psi));

      if(do_smoothing) {
	smoo.resize(ns, n);
	// Smoothing loop
	smoo.row(ns-1) = phi.row(ns-1);
	for(int t=1; t < ns; t++) {
	  int tt = ns - t;
	  // Time update using uniformization algorithm
	  matrix<Float> predrow = pred.row(tt);
	  for (int i=0; i < n; i++) {
	    predrow(0, i) += 1e-10;
	  }
	  matrix<Float> ratio = smoo.row(tt).cwiseQuotient(predrow);
	  //matrix<double> asd(1, n);
	  //for (int i=0; i < n; i++){
	  //asd(0, i) = std::isnan(asDouble(ratio(0, i)));
	  //cout 
	  //}
	  //matrix<double> asd = std::isnan(asDouble(ratio));
	  matrix<Float> ratiotmp = ratio;
	  if (solvetype == 1){
	    ratiotmp = forwardProject(ratio);
	  } else if (solvetype == 2) {
	    ratiotmp = forwardProjecti(ratio);
	  } else error("'solvetype' can by 1 or 2");
	  matrix<Float> post = phi.row(tt-1).cwiseProduct(ratiotmp);
	  post = post / (post.sum() + 1e-20);
	  smoo.row(tt-1) = post;
	}	
      }
      
      return ans;
    }
    
  };

  
  // ****** How to use it in TMB:
  // 1. Create an evaluator 'eval' for previous class
  template<class Float>
  Float eval(Float Dx, Float Dy) {
    filter_t<Float> f;
    return f.eval(Dx, Dy);
  }
  // 2. Run 'eval' through tiny_ad and obtain an atomic function
  //    'func'.  The '11' tells tiny_ad that we need
  //    derivatives wrt. both Dx and Dy.
  TMB_BIND_ATOMIC(func, 11, eval(x[0], x[1]))
  // 3. Create a more user-friendly version ('func' takes vector
  //    arguments and there's a final invisible argument that
  //    corresponds to the derivative order)
  template<class Type>
  Type hmm_nll(Type Dx, Type Dy) {
    vector<Type> args(3); // Last index reserved for derivative order
    args << Dx, Dy, 0;
    return shmm::func(CppAD::vector<Type>(args))[0];
  }
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(logDx);        // log diffusion in east-west (x) direction
  PARAMETER(logDy);        // log diffusion in north-south (y) direction

  // Transfer all constant data to namespace 'shmm'
  if(isDouble<Type>::value){
#define Type double

    DATA_MATRIX(datlik);     // Data likelihood
    DATA_INTEGER(solvetype); // Type of solver to use
    DATA_INTEGER(ns);        // Number of time steps of solver
    DATA_IVECTOR(iobs);      // Indices to which observations correspond

    DATA_SPARSE_MATRIX(I);   // Identity matrix
    DATA_SCALAR(dt);         // Time step
    DATA_INTEGER(m);         // Number of iterations of uniformization
    DATA_SPARSE_MATRIX(Sns); // North-south generator skeleton
    DATA_SPARSE_MATRIX(Sew); // East-west generator skeleton
#undef Type
    shmm::datlik = datlik;
    shmm::solvetype = solvetype;
    shmm::ns = ns;
    shmm::iobs = iobs;
    
    shmm::m = m;
    shmm::I = I;
    shmm::Sns = Sns;
    shmm::Sew = Sew;
    shmm::dt = dt;
    // Dirty trick didn't work for DATA_VECTOR:
    // DATA_VECTOR(lgam);       // Factorial
    shmm::lgam = asVector<double>(getListElement(objective_function::data,"lgam",&isNumeric));
  }

  Type Dx = exp(logDx);
  Type Dy = exp(logDy);

  Type ans = shmm::hmm_nll(Dx, Dy);

  /* With this construct 'dosmoo' flag becomes redundant. The branch
     is only entered by obj$report(). */
  if(isDouble<Type>::value) {
    shmm::filter_t<double> filter;
    filter.eval(asDouble(Dx), asDouble(Dy), true /* do_smoothing */ );
    // Store a subset of distribution for output
    int nobs = shmm::datlik.rows();
    int n = shmm::datlik.cols();
    int ns = shmm::ns;
    vector<int> iobs = shmm::iobs;
    matrix<double> smooout(nobs, n);
    matrix<double> phiout(nobs, n);
    matrix<double> predout(nobs, n);
    for(int t=0; t<ns; t++){
      if (iobs(t) > 0){
	int ind = (iobs(t)-1);
	smooout.row(ind) = filter.smoo.row(t);
	phiout.row(ind) = filter.phi.row(t);
	predout.row(ind) = filter.pred.row(t);
      }
    }
    vector<double> psi = filter.psi;
    // Reports
    REPORT(predout);
    REPORT(phiout);
    REPORT(psi);
    REPORT(smooout); 
  }
  
  return ans;

  /*
  // Viterbi
  // Forward sweep
  vector<Type> ksi1 = phi.row(0);
  vector<Type> ksi2 = ksi1;
  matrix<Type> zerovec(1, n);
  for(int i=0; i < n; i++){
    zerovec(0, i) = 0.0;
  }
  Type out = max(ksi1); 
  REPORT(out);
  for(int t=0; t<ns; t++){
    for(int i=0; i < n; i++){
      matrix<Type> vec = zerovec;
      vec(0, i) = phi(t, i);
      matrix<Type> newvec = shmm::ForwardProject(vec, Dx, Dy);
    }
  }
  // Backard sweep
  */


}

