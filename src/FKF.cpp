/*
 --------- References --------------------------------------------------------
 @article{jungbacker2007monte,
  title={Monte Carlo estimation for nonlinear non-Gaussian state space models},
  author={Jungbacker, Borus and Koopman, Siem Jan},
  journal={Biometrika},
  volume={94},
  number={4},
  pages={827--839},
  year={2007},
  publisher={Biometrika Trust}
}

@book{durbin2012time,
  title={Time series analysis by state space methods},
  author={Durbin, James and Koopman, Siem Jan},
  number={38},
  year={2012},
  publisher={Oxford University Press}
}
*/

#include <iostream>
#include <map>
#include <RcppArmadillo.h>
#include "../inst/include/GKF.h"

using namespace arma;
using namespace std;
using namespace Rcpp;
using namespace GKF;

// [[Rcpp::depends(RcppArmadillo)]]

// Smoothed matrix inverse
// when Condition number is larger than 1e12, the diagonal is inflated with
// epsilon =1e-6

mat InvMat(mat input){
  double condNbr = cond(input);
  
  if (input.n_rows != input.n_cols){
    return(pinv(input));
  }
  else if (condNbr > 1e12){
    mat eps ;
    eps.eye(input.n_rows,input.n_cols);
    return(inv(input + 1e-6*eps));
  }
  else {
    return(inv(input));
  }
}

// [[Rcpp::export]]
NumericMatrix InvMat_Rcpp(NumericMatrix input){
  mat input_ = as < mat > (input);
  mat res = InvMat(input_) ;
  return(wrap(res));
}


// Handle the conversion from R array to arma cube
// if a matrix is inputed it will be considered as 
// a cube with only one slice.

cube ConvertToCube(NumericVector myArray){
  NumericVector vecArray(myArray);
  IntegerVector arrayDims = vecArray.attr("dim");
  unsigned ThirdDim = 1; // matrix case
  
  if (arrayDims.size()==3){
    ThirdDim=arrayDims[2];
  }
  
  cube cubeArray(vecArray.begin(),arrayDims[0], arrayDims[1], ThirdDim, false);
  return(cubeArray);  
}

//  Sample from the multivariate normal distribution N(0,sigma)
//  sigma_sqrt is supposed to be the result of the chol decomposition
//  i.e sigma = sigma_sqrt*sigma_sqrt  where sigma_sqrt=Robust_chol(sigma)
// each col of the returned matrix is a sample from N(0,sigma)

mat mvnSample(unsigned n, mat sigma_sqrt){
  mat Y = sigma_sqrt*randn(sigma_sqrt.n_cols,n);
  return(Y);
}

// Wt is used to get red of missing values 
// ind contains the indexes of the non NA values
// nNA is the number of non NA and d is the original
// size of the observation vector
// one use of Wt is the following:
// yt*=Wt*yt where yt is a vector containing NA values
// yt* is a nNA vector with the the observed values only
// See Koopman(2012) p111 for more details

mat getWt(uvec ind, unsigned nNA,unsigned d){
  mat res(nNA,d,fill::zeros);
  for ( unsigned i=0; i < nNA; i++){
    res(i,ind(i))=1.0;
  }
  return(res);
}

// [[Rcpp::export]]
Rcpp::List nearPD_Rcpp(NumericMatrix x, 
		       NumericVector w,
		       bool corr, bool keepDiag,
		       bool EnforcePosDef,  // if TRUE, we insure the matrix is df
		       bool doSym, // symmetrize after tcrossprod()
		       bool ensureSymmetry, 
		       double eig_tol, // defines relative positiveness of eigenvalues compared to largest
		       double conv_tol, // convergence tolerance for algorithm
		       double posd_tol, // tolerance for enforcing positive definiteness
		       unsigned maxit, // maximum number of iterations allowed
		       bool chol // compute the Sqrt (Cholesky) decompoisition on the new matrix
		       ){
  mat x_ = as < mat > (x);
  vec w_ = as < vec > (w);  
  return(nearPD(x_,w_,corr,keepDiag,EnforcePosDef,doSym,ensureSymmetry, 
		eig_tol,conv_tol,posd_tol,maxit,chol)); 
}

// Rework of the Sqrt (Cholesky ?) decomposition
// we use eigenVal decomposition (we expect only symmetric matrices)
// to compute the 'square-root' of a matrix in the following way:
// if B = U*D*U.t(), where the matrix of eigenVectors(U=(U1|U2|...|Un))
// and D=diag(lambda_1,..,lambda_n), the 'square-root' of B is :
// sqrt_B = U_tilda*diag(sqrt(lambda_1),..,sqrt(lambda_k)); k<=n
// if the eigenvalues of B are found to be small (compared to eig_tol) or
// negative, the nearPD function is called to 'clean' it.
// the returned matrix is defined by R*R.t()=B (opposite of chol in Armadillo)

mat Robust_chol(mat B, double eig_tol = 1e-3){
  // get the eigen decomposition : B = eigvec*diag(eigval)*eigvec.t()
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, B);
  
  // apply the eps cut-off
  unsigned nB= B.n_rows;
  double eps = eigval(nB-1)*eig_tol; 
  
  uvec ind = find(eigval > eps);
  unsigned nN0 = ind.n_elem ;
  mat SqrtDec;
  
  if ( nN0 != B.n_rows ){
    vec w = ones < vec > (B.n_rows);
    Rcpp::List lS = nearPD(B,w,false,false,true,false,
			   false,1e-6, 1e-7,eig_tol);
    SqrtDec = as < mat > (lS["CholDec"]);
  }
  else {
    mat Wt = getWt(ind,nN0,nB);
    mat U_tilda = eigvec*Wt.t();
    mat D_sqrt(nN0,nN0,fill::zeros);
    D_sqrt.diag() = sqrt(eigval(ind));
    
    SqrtDec=U_tilda*D_sqrt;
  }
  return(SqrtDec);
}

// [[Rcpp::export]]
NumericMatrix Robust_chol_Rcpp(NumericMatrix B, double eig_tol=0.001){
  mat B_ = as < mat > (B);
  return(wrap(Robust_chol(B_,eig_tol)));
}

bool custom_isnan(double var)
{
  volatile double d = var;
  return d != d;
}

// Loglik(theta) computed iteratively
Rcpp::List ThetaLogLikUncond(vec ct, vec dt, vec thetaVal,vec xtildaPrev, 
			     mat Tt, mat Qt, mat Zt, unsigned iter){
  double qt=0;
  double lt=0;

  mat Xterm,xDx;
  
  vec mu=ct + Zt*dt;
  vec xtilda=InvMat(Zt)*(thetaVal-mu);
  
  if (iter==0){
    Xterm = InvMat(Tt)*xtilda;
  }
  else{
    Xterm = xtilda-Tt*xtildaPrev;
  }

  xDx=Xterm.t()*InvMat(Qt)*Xterm;
  if (!custom_isnan(xDx[0])) qt= xDx[0];
  
  lt= 2*log(det(Zt)) + log(det(Qt));
  if (custom_isnan(lt)) lt=0;
  
return Rcpp::List::create(Rcpp::Named("Qt")=qt,
				Rcpp::Named("lt")=lt,
				Rcpp::Named("xtilda")=xtilda);

}


/*
==============================================================================================================
--------------------------------- Fast Kalman Filter ---------------------------------------------------------
==============================================================================================================
 Kalman filter Implementation
 Same notation as fkf package except that Ht (GGt in fkf obs.var) and Qt (HHt in fkf state.var)
 are not supposed to be definite positive (only non-singular)
 The detail of the implementation follows Eq A1 for Jungbacker-Koopman(2007) with slightly modified 
 notations xt <-> yt and At <-> Ht
 the Inversion of Ft is made by a smoothed solve so LLt decomposition is not required anymore
 If ComputeThetaLik is true, we compute the loglikelihood log p(theta) for the signal theta
 In this the Rcpp list ThetaVal should contains an array (under the name theta), each slice of this array
 contains the signal sample where we compute the (log)likelihood. This (log)likelihood is computed iteratively. 
 Dimension:
 --- Inputs
 ** yt is dxn matrix ; yt can be entirely (one col) or partially (some elements of the col) NA
 ** at is mxn matrix (and hence a0) 
 ** n is the number of time observation
 ** P0 is mxm matrix
 ** dt is mxn matrix
 ** ct is dxn matrix
 ** Tt is mxmxn cube
 ** Zt is dxmxn cube
 ** Ht is dxdxn cube
 ** Qt is mxmxn cube
 ** ThetaVal a list containing an dxnxM array named theta 
 ** ComputeThetaLik is a bool
 --- Output
 an Rcpp list with the following components:
      ** at is mxn+1 matrix
      ** att is mxn matrix
      ** Pt is mxmxn+1 cube
      ** Ptt mxmxn cube
      ** vt is dxn matrix
      ** Ft is dxdxn cube
      ** Fti is dxdxn cube;  inverse of Ft
      ** Kt is mxdxn cube; definition different from FKF (See Koopman (2007) A1)
      ** logLik is n vec containing the contribution of each time point in the log-likelihood
 In order to obtain the entire log-likelihood, one has to sum the loglik vector

 ************* Validated *****************
*/

// [[Rcpp::export]]
Rcpp::List FKF(NumericVector a0_, NumericMatrix P0_, NumericMatrix dt_, NumericMatrix ct_,
	       NumericVector Tt_, NumericVector Zt_, NumericVector Ht_, NumericVector Qt_, // arrays
	       NumericMatrix yt_, bool ComputeThetaLik =false, Rcpp::List ThetaVal=0){

  //Conversion
  vec a0 = as < vec > (a0_);
  mat P0 = as < mat > (P0_);
  mat dt = as < mat > (dt_);
  mat ct = as < mat > (ct_);
  cube Tt = ConvertToCube(Tt_);
  cube Zt = ConvertToCube(Zt_);
  cube Ht = ConvertToCube(Ht_);
  cube Qt = ConvertToCube(Qt_);
  mat yt = as < mat > (yt_);

  // Extract dimension
  uword d,m,n;
  d = yt.n_rows;
  n = yt.n_cols;
  m = a0.n_elem;

  // Allocate memory for the results object
  mat at(m,n+1,fill::zeros);
  mat att(m,n,fill::zeros);
  cube Pt(m,m,n+1,fill::zeros);
  cube Ptt(m,m,n,fill::zeros);
  mat vt(d,n,fill::zeros);
  cube Ft(d,d,n,fill::zeros);
  cube Fti(d,d,n,fill::zeros); // inverse of Ft
  cube Kt(m,d,n,fill::zeros); // definition different from FKF (See Koopman Paper)
  vec logLik(n,fill::zeros);
  // loglik theta
  cube theta ;
  uword M = 0 ;
  mat thetalogLikD,xtilda;
  vec thetalogLik;

  if (ComputeThetaLik) {
    theta= ConvertToCube(ThetaVal["theta"]);
    M = theta.n_slices;
    thetalogLikD.zeros(n,M);
    thetalogLik.zeros(M);
    xtilda.zeros(d,M);
  }

  //initialisation
  at.col(0)=a0;
  Pt.slice(0)=P0;

  for (unsigned i=0; i<n; i++){
    // check Missing value
    uvec ind = find_finite(yt.col(i));
    unsigned nNA = ind.n_elem; // number of non na
    double logCte = -0.5*d*log( 2.0*datum::pi ); // was -0.5*nNA*log( 2.0*datum::pi )

    if (nNA == 0){ //===================== The entire observation is missing =========================
      // Updating Equations
      Ft.slice(i)= Ht.slice(i);
      Fti.slice(i)= InvMat(Ft.slice(i));
      Kt.slice(i).zeros();
      att.col(i)=at.col(i);
      Ptt.slice(i)=Pt.slice(i);

      // Prediction Equations
      at.col(i+1)= dt.col(i) + Tt.slice(i)*at.col(i); 
      Pt.slice(i+1)= Qt.slice(i) +  Tt.slice(i)*Pt.slice(i)*Tt.slice(i).t();
      // Likelihood
      logLik(i)=0; // just to comfirm
      
    }
    else if (nNA < d ){ //===================== Some observations are missing =========================
      // init modified size object
      vec yti0 = yt.col(i);
      vec yti = yti0.elem(ind);
      vec cti0 = ct.col(i);
      vec cti = cti0.elem(ind);
      mat Tti=Tt.slice(i);
      
      mat Wt = getWt(ind,nNA,d);
      mat Zti = Wt*Zt.slice(i) ;
      mat Ztit = Zti.t() ;
      mat Hti = Wt*Ht.slice(i)*Wt.t();
      mat Pti = Pt.slice(i);

      // Updating Equations
      vec vti=yti-cti-Zti*at.col(i);
      mat Ftin= Hti + Zti*Pti*Ztit;
      mat Ftiin = InvMat(Ftin);
      mat Kti= Pti*Ztit*Ftiin;
      att.col(i)=at.col(i) + Kti*vti;
      Ptt.slice(i) = Pti - Pti*Ztit*Kti.t();
      // Prediction Equations
      at.col(i+1)= dt.col(i) + Tt.slice(i)*att.col(i) ;
      Pt.slice(i+1)= Qt.slice(i) + Tt.slice(i)*Ptt.slice(i)*Tt.slice(i).t();
      
      // Likelihood
      mat vFv = vti.t()*Ftiin*vti;
      logLik(i)= logCte  - 0.5 *log(det(Ftin)) - 0.5*vFv[0];
      // theta Likelihood
      if (ComputeThetaLik){
	Rcpp::List ltheta ;
	
	for (unsigned s =0; s < M ; s++){
	  mat QQi = Qt.slice(i);
	  if (i==0) QQi=P0;
	  vec Vtheta = theta(span(),span(i,i),span(s,s));
	  ltheta = ThetaLogLikUncond(ct.col(i),dt.col(i),Vtheta,xtilda.col(s),
				     Tti,QQi,Zt.slice(i),i);
	  
	  double qt = as < double > (ltheta["Qt"]);
	  double llt = as < double > (ltheta["lt"]);
	  thetalogLikD(i,s)= logCte  - 0.5 *llt - 0.5*qt;
	  thetalogLik(s) += thetalogLikD(i,s) ;
	  xtilda.col(s)= as < vec > (ltheta["xtilda"]);
	}
      }

      // store object with nan elsewhere
      vt.col(i)= vt.col(i)*datum::nan ;
      Ft.slice(i)= Ft.slice(i)*datum::nan;
      Fti.slice(i)= Fti.slice(i)*datum::nan;
      Kt.slice(i)= Kt.slice(i)*datum::nan;

      vt(span(0,nNA-1),i)=vti;
      Ft(span(0,nNA-1),span(0,nNA-1),span(i,i))=Ftin;
      Fti(span(0,nNA-1),span(0,nNA-1),span(i,i))=Ftiin;
      Kt(span(),span(0,nNA-1),span(i,i))=Tt.slice(i)*Kti; // update Kt to meet koopman def
    }
    else { //===================== No missing values ==================================================
      mat Zti=Zt.slice(i);
      mat Ztit= Zti.t();
      mat Pti=Pt.slice(i);
      mat Tti=Tt.slice(i);
      vec cti =ct.col(i);
      vec dti= dt.col(i);
      // Updating Equations
      vec vti = yt.col(i)-ct.col(i)-Zti*at.col(i);
      vt.col(i)= vti ;
      Ft.slice(i)= Ht.slice(i) + Zti*Pti*Ztit;
      mat Ftii = InvMat(Ft.slice(i));
      Fti.slice(i) = Ftii;
      mat Kti= Pti*Ztit*Ftii;
      att.col(i)= at.col(i) + Kti*vti;
      Ptt.slice(i)=Pti - Pti*Ztit*Kti.t();
      // Prediction Equations
      at.col(i+1)= dt.col(i) + Tti*att.col(i);
      Pt.slice(i+1)= Qt.slice(i) +  Tti*Ptt.slice(i)*Tti.t();
      // update Kt to meet koopman def
      Kt.slice(i)= Tti*Kti;
      // Likelihood
      mat vFv = vti.t()*Ftii*vti;
      logLik(i)= logCte  - 0.5 *log(det(Ft.slice(i))) - 0.5*vFv[0];
      // theta Likelihood
      if (ComputeThetaLik){
	Rcpp::List ltheta ;
	
	for (unsigned s =0; s < M ; s++){
	  mat QQi = Qt.slice(i);
	  if (i==0) QQi=P0;
	  vec Vtheta = theta(span(),span(i,i),span(s,s));
	  ltheta = ThetaLogLikUncond(cti,dti,Vtheta,xtilda.col(s),
				     Tti,QQi,Zti,i);
	  double qt = as < double > (ltheta["Qt"]);
	  double llt = as < double > (ltheta["lt"]);
	  thetalogLikD(i,s)= logCte  - 0.5 *llt - 0.5*qt;
	  thetalogLik(s) += thetalogLikD(i,s) ;
	  xtilda.col(s)= as < vec > (ltheta["xtilda"]);
	}
      }
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("at")=at,
			    Rcpp::Named("Pt")=Pt,
			    Rcpp::Named("att")=att,
			    Rcpp::Named("Ptt")=Ptt,
			    Rcpp::Named("vt")=vt,
			    Rcpp::Named("Ft")=Ft,
			    Rcpp::Named("Fti")=Fti,
			    Rcpp::Named("Kt")=Kt,
			    Rcpp::Named("logLik")=sum(logLik),
			    Rcpp::Named("logLik_details")=logLik,
			    Rcpp::Named("thetalogLik")=thetalogLik,
			    Rcpp::Named("thetalogLik_details")=thetalogLikD);
}

/*
==============================================================================================================
-------------------------Newton-Raphson Updating Step ---------------------------------------------------------
==============================================================================================================
 Eq A10-A11 for Jungbacker-Koopman(2007) with slightly modified notation
 xt <-> yt and At <-> Ht
 It requires the implementation of the KF first
 The recursion is backward from n,..,1
 returns the matrix g+=[g1+|g2+|....|gn+]
 Only partial missing values in the observation vector are expected
  ************* Validated *****************
*/

// [[Rcpp::export]]
NumericMatrix NRUpdatingStep(NumericVector a0_, NumericMatrix P0_, NumericMatrix dt_, NumericMatrix ct_,
			     NumericVector Tt_, NumericVector Zt_, NumericVector Ht_, NumericVector Qt_, // arrays
			     NumericMatrix yt_){
  
  // Convert to Arma Objects
  cube Tt = ConvertToCube(Tt_);
  cube Zt = ConvertToCube(Zt_);
  cube Ht = ConvertToCube(Ht_);
  mat yt = as < mat > (yt_);

  // Call Kalman Filter and Extract needed objects
  Rcpp::List FKFRes = FKF(a0_,P0_,dt_,ct_,Tt_,Zt_,Ht_,Qt_,yt_);
  cube Kt = ConvertToCube(FKFRes["Kt"]);
  cube Fti= ConvertToCube(FKFRes["Fti"]);
  mat vt = as < mat > (FKFRes["vt"]);

  // Allocate Memory to the new objects
  // Extract dimension
  unsigned d,m,n,j;
  d = yt.n_rows;
  n = yt.n_cols;
  m = a0_.size();

  // Allocate memory for the results object
  vec et(d,fill::zeros);
  vec st(m,fill::zeros); // already initialised at time n
  mat gtp(d,n,fill::zeros); // same size as yt
  
  for (unsigned i=0; i<n; i++){
    // the backward index
    j=n-1-i;
    
    // check Missing value
    uvec ind = find_finite(yt.col(j));
    unsigned nNA = ind.n_elem; // number of not na
    
    if (nNA < d) {
      // init modified size object
      vec ytj0 = yt.col(j);
      vec ytj = ytj0.elem(ind);
      vec vtj = vt(span(0,nNA-1),j); // due to our KF convention. NA stored at the tail
      
      mat Wt_ = getWt(ind,nNA, d);
      mat Ztj = Wt_*Zt.slice(j) ;
      mat Htj = Wt_*Ht.slice(j)*Wt_.t();
      mat Ftij = Fti(span(0,nNA-1),span(0,nNA-1),span(j,j));
      mat Ktj = Kt(span(),span(0,nNA-1),span(j,j));
            
      gtp.col(j)=gtp.col(j)*datum::nan ;

      et(span(0,nNA-1))= Ftij*vtj - Ktj.t()*st ;
      gtp(span(0,nNA-1),span(j,j))= ytj - Htj*et(span(0,nNA-1));
      st=Ztj.t()*et(span(0,nNA-1))+ Tt.slice(j).t()*st;
    }
    else{
      // Updating equation, start with s
      et= Fti.slice(j)*vt.col(j) - Kt.slice(j).t()*st;
      gtp.col(j)= yt.col(j) - Ht.slice(j)*et;
      st=Zt.slice(j).t()*et + Tt.slice(j).t()*st;
    }
  }
 
  return(wrap(gtp));
}


/*
==============================================================================================================
------------------------- Signal Sampling ---------------------------------------------------------
==============================================================================================================
 Theorem2 (Eq 17) in Jungbacker-Koopman(2007) with slightly modified notation xt <-> yt and At <-> Ht
 It requires the implementation of the KF first
 The recursion is backward from n,..,1
 On top of the usual KF inputs, user needs to specify:
 *** thetaHat : dxn matrix containing the posterior (signal) mode (each col is a an obs time)
 *** M        : the number of sample size required
 returns a cube (3D-array) where each slice is a sample theta^i ~f(theta|y) = [theta_hat_1+u1|...|theta_hat_n+un]
         an M-vector of (log)likelihood
*/
// [[Rcpp::export]]
Rcpp::List thetaSampling(NumericVector a0_, NumericMatrix P0_, NumericMatrix dt_, NumericMatrix ct_,
			 NumericVector Tt_, NumericVector Zt_, NumericVector Ht_, NumericVector Qt_, // arrays
			 NumericMatrix yt_, NumericMatrix thetaHat_, unsigned M){
  
  // Convert to Arma Objects
  cube Tt = ConvertToCube(Tt_);
  cube Zt = ConvertToCube(Zt_);
  cube Ht = ConvertToCube(Ht_);
  mat P0 = as < mat > (P0_);
  mat dt = as < mat > (dt_);
  mat ct = as < mat > (ct_);
  mat yt = as < mat > (yt_);
  mat thetaHat = as < mat > (thetaHat_);

  // Call Kalman Filter and Extract needed objects
  Rcpp::List FKFRes = FKF(a0_,P0_,dt_,ct_,Tt_,Zt_,Ht_,Qt_,yt_);
  cube Kt = ConvertToCube(FKFRes["Kt"]);
  cube Fti= ConvertToCube(FKFRes["Fti"]);
  mat vt = as < mat > (FKFRes["vt"]);

  // Allocate Memory to the new objects
  // Extract dimension
  unsigned d,m,n,j,s;
  mat tmp(1,1,fill::zeros);
  d = yt.n_rows;
  n = yt.n_cols;
  m = a0_.size();

  // Allocate memory for the results object
  vec wt(d,fill::zeros);
  mat rt(m,M,fill::zeros);// already initialised at time n
  mat WtM(d,M,fill::zeros); // M vector same size as yt
  mat Rt(d,m,fill::zeros);
  mat Nt(m,m,fill::zeros); // already initialised at time n
  mat Ct(d,d,fill::zeros);
  // Likelihood specific
  mat Bt(d,d,fill::zeros);
  mat Bti(d,d,fill::zeros);
  vec b(d,fill::zeros);

  cube theta(d,n,M,fill::zeros); // result
  vec logLik(M,fill::ones); //result
  // init the log likelihood vector with the cte mn/2*log(2*Pi)
  double cte = -0.5*d*n*log( 2.0*datum::pi);
  logLik =logLik*cte;

  for (unsigned i=0; i<n; i++){
    // the backward index
    j=n-1-i;
    
    // check Missing value
    uvec ind = find_finite(yt.col(j));
    unsigned nNA = ind.n_elem; // number of not na
    
    if (nNA < d) {
      vec vtj = vt(span(0,nNA-1),j);
      mat Wt_ = getWt(ind,nNA, d);
      mat Ztj = Wt_*Zt.slice(j) ;
      mat Htj = Wt_*Ht.slice(j)*Wt_.t();
      mat Ftij = Fti(span(0,nNA-1),span(0,nNA-1),span(j,j));
      mat Ktj = Kt(span(),span(0,nNA-1),span(j,j));
      mat Htnj = InvMat(Htj);
      mat Ttj = Tt.slice(j);

      //updating equations
      mat Ctj  = Htnj - Ftij - Ktj.t()*Nt*Ktj;
      mat Rtj  = InvMat(Ctj)*(Htnj*Ztj -Ktj.t()*Nt*Ttj);
      mat Btj  = Robust_chol(Ctj);
      mat Btji = InvMat(Btji); 
      vec btj(nNA,fill::zeros);
    
      // Simulate from Mvn using Ct
      WtM(span(0,nNA-1),span()) = mvnSample(M,Btj);
      // for each simulated sample compute u_j and add it to thetaHat_j
      for (s=0;s<M;s++){
	vec wtj=WtM(span(0,nNA-1),s);
	vec ut = Htj*(wtj + Ftij*vtj - Ktj.t()*rt.col(s));
	
	theta(span(),span(j,j),span(s,s))= theta(span(),span(j,j),span(s,s))*datum::nan; 
	theta(span(0,nNA-1),span(j,j),span(s,s)) = thetaHat(span(0,nNA-1),j) + ut ;
	// update rt
	rt.col(s)= Ztj.t()*Htnj*ut-Rtj.t()*wtj+Ttj.t()*rt.col(s) ;
	// ------------ update the Likelihood
	// f(theta;y)
	btj= Btji*wtj;
	tmp = -0.5*(btj.t()*btj) ;
	tmp(0,0) = tmp(0,0) -log(abs(det(Htj))) - log(det(Btj))  ;
	logLik(s)= logLik(s) + tmp(0,0);
      }
      // update Nt
      Nt=Rtj.t()*Ctj*Rtj - Ztj.t()*Htnj*Ztj + Ttj.t()*Nt*Ttj;
    }
    else{
      // Compute Ct and Rt
      mat Htj  = Ht.slice(j);
      mat Htnj = InvMat(Htj);
      mat Ftij = Fti.slice(j);
      mat Ktj = Kt.slice(j);
      mat Ztj = Zt.slice(j);
      mat Ttj = Tt.slice(j);
      
      Ct  = Htnj  - Ftij - Ktj.t()*Nt*Ktj;
      Rt  = InvMat(Ct)*(Htnj*Ztj - Ktj.t()*Nt*Ttj);
      Bt  = Robust_chol(Ct);
      
      Bti = InvMat(Bt) ;
      // Simulate from Mvn using Ct
      WtM = mvnSample(M,Bt);
      // for each simulated sample compute u_j and add it to thetaHat_j
      for (s=0;s<M;s++){
	wt=WtM.col(s);
	vec ut = Htj*(wt + Ftij*vt.col(j) - Ktj.t()*rt.col(s));
	theta(span(),span(j,j),span(s,s))= thetaHat.col(j) + ut;
	// update rt
	rt.col(s)=Ztj.t()*Htnj*ut-Rt.t()*wt+Ttj.t()*rt.col(s);
	// update the Likelihood
	b=Bti*wt;
	tmp = -0.5*(b.t()*b) ;
	// the abs for Bt is just to prevent numerical issue with negative det
	tmp(0,0) = tmp(0,0) -log(abs(det(Htj))) - log(abs(det(Bt)))  ;
	logLik(s)= logLik(s) + tmp(0,0);
    }
      // update Nt
      Nt=Rt.t()*Ct*Rt - Ztj.t()*Htnj*Ztj + Ttj.t()*Nt*Ttj;
    }
  } 
  return Rcpp::List::create(Rcpp::Named("theta")=theta,
			    Rcpp::Named("logLik")=logLik);
}


/*
==============================================================================================================
------------------------- Gaussian Signal Smoothing ----------------------------------------------------------
==============================================================================================================
 Signal smoothing as described in Koopman(2012, chapter 4.5.3)
 We use 2nd version of (4.69) to compute eta_hat and then use the following definition of theta_hat
 theta_hat= yt - eta_hat ;; eta_hat=E[eta|Yn]
 We only return theta_Hat as well as its variance in a Rcpp::List
 Partly missing observation in yt are allowed
 ************* Validated *****************
*/
// [[Rcpp::export]]
Rcpp::List GaussianSignalSmoothing(NumericVector a0_,NumericMatrix P0_,NumericMatrix dt_,NumericMatrix ct_,
				   NumericVector Tt_,NumericVector Zt_,NumericVector Ht_,NumericVector Qt_,//arrays
				   NumericMatrix yt_){

  cube Tt = ConvertToCube(Tt_);
  cube Zt = ConvertToCube(Zt_);
  cube Ht = ConvertToCube(Ht_);
  mat yt = as < mat > (yt_);
  
  // Call Kalman Filter and Extract needed objects
  Rcpp::List FKFRes = FKF(a0_,P0_,dt_,ct_,Tt_,Zt_,Ht_,Qt_,yt_);
  cube Kt = ConvertToCube(FKFRes["Kt"]);
  cube Fti= ConvertToCube(FKFRes["Fti"]);
  mat vt = as < mat > (FKFRes["vt"]);

  // Allocate Memory to the new objects
  // Extract dimension
  unsigned d,m,n,j;
  d = yt.n_rows;
  n = yt.n_cols;
  m = a0_.size();

  // Allocate memory for the results object
  // temp
  vec rt(m,fill::zeros); //already initialised
  mat Nt(m,m,fill::zeros); // already initialised at time n
  vec ut(m,fill::zeros); 
  vec eta_hat(d,fill::zeros);
  mat Dt(d,d,fill::zeros);
  // result
  mat theta_hat(d,n,fill::zeros); // same size as yt
  cube Vtheta_hat(d,d,n,fill::zeros); // variance

  for (unsigned i=0; i<n; i++){
    // the backward index
    j=n-1-i;
    uvec ind = find_finite(yt.col(j));
    unsigned nNA = ind.n_elem; // number of non na
    if (nNA < d) { 
      // init modified size object
      vec ytj0 = yt.col(j);
      vec ytj = ytj0.elem(ind);
      vec vtj = vt(span(0,nNA-1),j); // KF convention
      
      mat Wt_ = getWt(ind,nNA, d);
      mat Ztj = Wt_*Zt.slice(j) ;
      mat Htj = Wt_*Ht.slice(j)*Wt_.t();
      mat Ftij = Fti(span(0,nNA-1),span(0,nNA-1),span(j,j));
      mat Ktj = Kt(span(),span(0,nNA-1),span(j,j));
      
      // init theta_hat
      theta_hat.col(j)=theta_hat.col(j)*datum::nan ;
      Vtheta_hat.slice(j)=Vtheta_hat.slice(j)*datum::nan ;
      // updating temp equations
      vec utj= Ftij*vtj - Ktj.t()*rt;
      mat Dt = Ftij + Ktj.t()*Nt*Ktj;
      // compute qtities of interest
      vec etaj_hat = Htj*utj;
      theta_hat(span(0,nNA),j)=ytj-etaj_hat;
      Vtheta_hat(span(0,nNA),span(0,nNA),span(j,j)) = Htj - Htj.t()*Dt*Htj;
      // update temp
      rt = Ztj.t()*utj + Tt.slice(j).t()*rt;
      Nt = Ztj.t()*Dt*Ztj + Tt.slice(j).t()*Nt*Tt.slice(j) 
	-Ztj.t()*Ktj.t()*Nt*Tt.slice(j) - Tt.slice(j)*Nt*Ktj*Ztj;
    }
    else{
      // compute temp
      ut = Fti.slice(j)*vt.col(j) - Kt.slice(j).t()*rt;
      Dt=  Fti.slice(j) + Kt.slice(j).t()*Nt*Kt.slice(j) ;
      
      // compute qtities of interest
      eta_hat = Ht.slice(j)*ut;
      theta_hat.col(j)=yt.col(j)-eta_hat;
      Vtheta_hat.slice(j) = Ht.slice(j) - Ht.slice(j)*Dt*Ht.slice(j);
      // update temp
      rt = Zt.slice(j).t()*ut + Tt.slice(j).t()*rt ;
      Nt = Zt.slice(j).t()*Dt*Zt.slice(j) + Tt.slice(j).t()*Nt*Tt.slice(j).t() -
	Zt.slice(j).t()*Kt.slice(j).t()*Nt*Tt.slice(j) - Tt.slice(j)*Nt*Kt.slice(j)*Zt.slice(j); 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("theta_hat")=theta_hat,
			    Rcpp::Named("Vtheta_hat")=Vtheta_hat);
}


/*
==============================================================================================================
------------------------- Gaussian Signal Sampling ---------------------------------------------------------
==============================================================================================================
 Gaussin (Signal) simulation smoother by mean correction as described in Koopman(2012) 4.92
 theta_tilda ~p(theta|Yn) is obtained as follow:
 **** Draw eta+_t ~ N(0,Ht), nu+_t ~ N(0,Qt)
 **** Draw alpha+_t using the state recursion starting with alpha+_1 ~N(a0,P0)
 **** Obtain theta+_t = ct +Zt*alpha+_t ;; y+_t = theta+_t + eta+_t
 **** theta_hat = E[theta|Yn] , theta+_hat = E[theta|Y+] obtained from the Kalman filter smoother
 **** theta_tilda = theta+ - theta+_hat + theta_hat

 result is a Rcpp::List with the following
 *** theta_tilda: cube with M slice, each slice is a dxn matrix representing the sample
 *** logLik     : M vector containing the log-likelihood of each sample
 Note that this log-likelihood is given by: 
 lg = -dn/2 log(2*Pi) -0.5*sum_t=1^n log|Vt|-sum_t=1^n 0.5*(theta_tilda_t-thetaHat_t)'*Vt^-1*(theta_tilda_t-thetaHat_t)
 where thetaHat and Vt are given by the GaussianSignalSmoothing and theta_tilda is the actual sample
  ************* Validated *****************
*/

// [[Rcpp::export]]
Rcpp::List GaussianthetaSampling(NumericVector a0_, NumericMatrix P0_, NumericMatrix dt_, NumericMatrix ct_,
				 NumericVector Tt_, NumericVector Zt_, NumericVector Ht_, NumericVector Qt_,//arrays
				 NumericMatrix yt_, unsigned M){
  
  // Convert to Arma Objects
  vec a0 = as < vec > (a0_);
  mat P0 = as < mat > (P0_);
  mat dt = as < mat > (dt_);
  mat ct = as < mat > (ct_);
  cube Tt = ConvertToCube(Tt_);
  cube Zt = ConvertToCube(Zt_);
  cube Ht = ConvertToCube(Ht_);
  cube Qt = ConvertToCube(Qt_);
  mat yt = as < mat > (yt_);

  // Allocate Memory to the new objects
  // Extract dimension
  unsigned d,m,n,j,s;
  d = yt.n_rows;
  n = yt.n_cols;
  m = a0_.size();

  // Allocate memory for the results object
  map<unsigned,mat> Vti; //inverse of the smoothed variance matrix
  mat etaPlus(d,M,fill::zeros);
  mat nuPlus(m,M,fill::zeros);
  cube thetaPlus(d,n,M,fill::zeros); // used to get theta_tilda
  cube ytPlus(d,n,M,fill::zeros); // used for obtaining thetaHatPlus

  // get the Smoothed values with the observation Yn
  Rcpp::List Smoothed =  GaussianSignalSmoothing(a0_,P0_,dt_,ct_,Tt_,Zt_,Ht_,Qt_,yt_);
  mat thetaHat = as < mat > (Smoothed["theta_hat"]);
  cube Vtheta = ConvertToCube(Smoothed["Vtheta_hat"]);
  
  cube theta_tilda(d,n,M,fill::zeros); // result
  vec logLik(M,fill::zeros); // result
  logLik = logLik -0.5*d*n*log(2*datum::pi); // likelihood cte
  // init Simulated vals
  mat alphaPlus = mvnSample(M,Robust_chol(P0)); //mxM
  alphaPlus.each_col() += a0;
  
  // ------ construct thetaPlus ---------------------
  for (j=0; j<n; j++){
    // check Missing value
    uvec ind = find_finite(yt.col(j));
    unsigned nNA = ind.n_elem; // number of not na
    nuPlus = mvnSample(M,Robust_chol(Qt.slice(j)));
  
    if (nNA < d) {
      // size adapted necessary objects
      mat Wt_ = getWt(ind,nNA, d);
      mat Ztj = Wt_*Zt.slice(j) ;
      mat Htj = Wt_*Ht.slice(j)*Wt_.t();
      etaPlus(span(0,nNA),span()) = mvnSample(M,Robust_chol(Htj));

      for (s=0;s<M;s++){
	// init the particular column with nan
	thetaPlus(span(),span(j,j),span(s,s)) = thetaPlus(span(),span(j,j),span(s,s)) * datum::nan ;
	// update the non-nan part
	thetaPlus(span(0,nNA),span(j,j),span(s,s)) = ct(span(0,nNA),j) + Ztj*alphaPlus.col(s);
	mat ttp = thetaPlus(span(0,nNA),span(j,j),span(s,s));
	ytPlus(span(0,nNA),span(j,j),span(s,s)) = ttp + etaPlus(span(0,nNA),s);
	// update alpha for each sample (at time j) 
	alphaPlus.col(s) = dt.col(j) + Tt.slice(j)*alphaPlus.col(s) + nuPlus.col(s);
      }
      mat Vth = Vtheta(span(0,nNA),span(0,nNA),span(j,j));
      // update the likelihood with the det term
      logLik = logLik - 0.5*log(det(Vth));
      // store the inverse of the Vtheta matrix
      Vti[j]=InvMat(Vth);
    }
    else{
      etaPlus= mvnSample(M,Robust_chol(Ht.slice(j)));  
       for (s=0;s<M;s++){
	 thetaPlus(span(),span(j,j),span(s,s)) = ct.col(j) + Zt.slice(j)*alphaPlus.col(s);
	 mat tp = thetaPlus(span(),span(j,j),span(s,s));
	 ytPlus(span(),span(j,j),span(s,s)) = tp +  etaPlus.col(s);
	 alphaPlus.col(s) = dt.col(j) + Tt.slice(j)*alphaPlus.col(s) + nuPlus.col(s);
       }
       // update the likelihood with the det term
       logLik = logLik - 0.5*log(det(Vtheta.slice(j)));
       // store the inverse of the Vtheta matrix
       Vti[j]=InvMat(Vtheta.slice(j));
    }
  }

  
  // ---------------------- Construct the actual sample  -----------------------
  // We had to do it outside the first loop because we need the entire ytPlus 
  // to compute the simulated sample (ytPlus gives ThetaHatPlus
  for (s=0;s<M;s++){
    // construct actual sample theta_tilda
    Rcpp::List SmoothedPlus = GaussianSignalSmoothing(a0_,P0_,dt_,ct_,Tt_,Zt_,Ht_,Qt_,wrap(ytPlus.slice(s)));
    mat thetaHatPlus = as < mat > (SmoothedPlus["theta_hat"]);
    theta_tilda.slice(s) = thetaPlus.slice(s) - thetaHatPlus + thetaHat ;
    
    // update the likelihood : for each sample, we update the likelihood with the quadratic form
    // which depends on the sampled theta_tilda constructed just before.
    for (j=0; j<n; j++){
      // check Missing value
      uvec ind = find_finite(yt.col(j));
      unsigned nNA = ind.n_elem; // number of non na
      
      if (nNA < d) {
	vec th_tilda = theta_tilda(span(0,nNA),span(j,j),span(s,s));
	vec th_Hat = thetaHat(span(0,nNA),j) ;
	vec DiffTheta = th_tilda  - th_Hat;
	mat Vtit = Vti[j];
	mat tVt = DiffTheta.t()*Vtit*DiffTheta;
	logLik(s) = logLik(s) - 0.5*tVt(0,0);	
      }
      else{
	vec th_tilda = theta_tilda(span(),span(j,j),span(s,s));
	vec DiffTheta =  th_tilda - thetaHat.col(j);
	mat Vtit = Vti[j];
	mat tVt = DiffTheta.t()*Vtit*DiffTheta ;
	logLik(s)=logLik(s) - 0.5*tVt(0,0);
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("theta")=theta_tilda,
			    Rcpp::Named("logLik")=logLik);
			    
}

