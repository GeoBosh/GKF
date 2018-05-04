#include <iostream>
#include <map>
#include <RcppArmadillo.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string> // this should be already included in <sstream>

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::interfaces(cpp)]]

mat getWt_(uvec ind, unsigned nNA,unsigned d){
  mat res(nNA,d,fill::zeros);
  for ( unsigned i=0; i < nNA; i++){
    res(i,ind(i))=1.0;
  }
  return(res);
}

mat ProjectSpd(mat B, double eig_tol){
  mat D,Q ;
  vec d ;
  
  eig_sym(d, Q, B);
  unsigned nB= B.n_rows;
  double eps = d(nB-1)*eig_tol; 
  uvec ind = find(d > eps);
  unsigned nN0 = ind.n_elem ;
  
  if (nN0==0){
    Rf_warning("Matrix seems negative semi-definite");
    mat Empty ;
    return Empty ;
  } 
  
  mat Wt = getWt_(ind,nN0,nB);
  mat Q_tilda = Q*Wt.t();
  D.zeros(nN0,nN0);
  D.diag() = d(ind);
  return(Q_tilda*D*Q_tilda.t()); 
}

mat ProjectUnitDiag(mat X, bool doSym, bool corr, bool keepDiag ,vec diagX0,
		    unsigned n){
  mat Y=X;
  
  if(doSym)
    Y = 0.5*(X + X.t());
  if(corr) 
    Y.diag()= ones(n)  ;
  else if(keepDiag) 
    Y.diag() = diagX0 ;
  
  return(Y);
}

void EnforceDefinitePositive(mat &X, double posd_tol, unsigned n){
  mat D,Q ;
  vec d ;
  // begin from posdefify(sfsmisc)
  eig_sym(d, Q, X);
  double Eps = posd_tol * abs(d(d.n_elem -1)); //here we use the posd.tol sensibility
  if (d(0) < Eps) { // if smallest eigenvalue smaller than eps
    uvec indSmall = find(d < Eps) ; 
    d.elem(indSmall)= Eps*ones<vec>(indSmall.n_elem); // replace the small eig val by Eps
    vec odiag = diagvec(X) ;
    D.zeros(n,n);
    D.diag()=d;
    X= Q*D*Q.t(); // rebuild X
    vec VEps = Eps*ones<vec>(n) ;
    for (unsigned i=0; i<n;i++) {
      if (VEps(i) < odiag(i))
	VEps(i) = odiag(i);
	}
    D.diag() = sqrt(VEps / X.diag());
    X= D*X*D;
  }
}

mat ComputeChol (mat X, vec &d, unsigned n){
  mat Q;
  mat D_sqrt(n,n,fill::zeros);
  
  eig_sym(d, Q, X);
   
  D_sqrt.diag() = sqrt(d);
  return(Q*D_sqrt);
}
 
// [[Rcpp::export]]
Rcpp::List nearPD(arma::mat x_, 
		  arma::vec w_, //  is a vector defining a diagonal weight matrix diag(w).
		  bool corr = false, bool keepDiag = false,
		  bool EnforcePosDef = true,  // if TRUE, we insure the matrix is df
		  bool doSym = false, // symmetrize after tcrossprod()
		  bool ensureSymmetry = false, 
		  double eig_tol = 1e-6, // defines relative positiveness of eigenvalues compared to largest
		  double conv_tol  = 1e-7, // convergence tolerance for algorithm
		  double posd_tol  = 1e-3, // tolerance for enforcing positive definiteness
		  unsigned maxit    = 100, // maximum number of iterations allowed
		  bool chol = true // approximate cholesky decomposition
		  ){
  mat X,Xold,Y,Yold,R,B;
    
  if(ensureSymmetry) { // only if needed/wanted ...
    x_ = 0.5*(x_ + x_.t());
   }
  
  unsigned n = x_.n_rows;
  vec diagX0 ;
  if(keepDiag) diagX0 = diagvec(x_) ;

  // init matrices
  // --- for Dykstra's correction,  allow memory
  mat D_S(x_.n_rows,x_.n_cols,fill::zeros);
  X=x_;
  Y=X;

  // init vars
  unsigned iter =0 ;
  bool converged = false;
  double conv = datum::inf;
  vec rel_diff(3,fill::ones);

  // weight matrix
  mat Whalf=sqrt(w_*w_.t());

   while (iter < maxit && !converged) {
     Xold=X ;
     R = Y - D_S ;
     // project onto PSD matrices  X_k  =  P_S (R_k)
     mat B = Whalf % R ; 
     X = ProjectSpd(B,eig_tol);
     X = X / Whalf; 
     // ------------------------------------------------------- 
     D_S = X-R; //update Dykstra's correction D_S = \Delta S_k
      
     Yold=Y;
     //project onto symmetric and possibly 'given diag' matrices:
     Y= ProjectUnitDiag(X,doSym,corr,keepDiag,diagX0,n);
      
     rel_diff(0) = norm(X-Xold, "fro") / norm(X, "fro"); // Xdiff
     rel_diff(1) = norm(Y-Yold, "fro") / norm(Y, "fro"); //Ydiff
     rel_diff(2) = norm(Y-X, "fro") / norm(Y, "fro"); // XYdiff
     conv = max(rel_diff); 
     iter += 1 ;
     
     converged = (conv <= conv_tol) ;
     X=Y;
   }
   
   string msg;
   if(!converged){
     msg = "nearPD()' did not converge in " + 
       static_cast<ostringstream*>( &(ostringstream() << iter) )->str() + 
       " iterations" ;

     Rf_warning(msg.c_str());   
   }

   // if do2eigen, we remove small and negative eigenvalues
   if(EnforcePosDef) 
     EnforceDefinitePositive(X,posd_tol,n) ;   
   
   Y= ProjectUnitDiag(X,doSym,corr,keepDiag,diagX0,n);
   X=Y;
   
   // add the approximate Cholesky decomposition here if needed
   mat sqX ;
   vec d;
   if (chol) 
     sqX = ComputeChol(X,d,n);
   
   // compute the diff norm between original and approx matrix
   double diffNorm = norm(x_-X,"fro");
   
   return Rcpp::List::create(Rcpp::Named("mat")=X,
			     Rcpp::Named("eigenvalues")=d,
			     Rcpp::Named("corr")=corr,
			     Rcpp::Named("normF")=diffNorm,
			     Rcpp::Named("iterations")=iter,
			     Rcpp::Named("rel.tol")= conv,
			     Rcpp::Named("converged")= converged,
			     Rcpp::Named("CholDec")= sqX
			     );
}
 
