#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include "gaston.h"
#ifndef GASTONRELIK
#define GASTONRELIK
//#define ANY(_X_) (std::any_of(_X_.begin(), _X_.end(), [](bool x) {return x;})) 

using namespace Rcpp;
using namespace Eigen;

template<typename T1, typename T2, typename T3, typename A>
double re_likelihood(const Eigen::MatrixBase<T1> & y, const Eigen::MatrixBase<T3> & x, const std::vector<T2,A> & K, 
               NumericVector theta) {
  int n(y.rows()), p(x.cols());
  int s(K.size());

  MatrixXd V(n,n), P(n,n), Vi(n,n);
  MatrixXd XViX(p,p), XViX_i(p,p);
  MatrixXd ViX(n,p);
  VectorXd Py(n);

  double log_detV, detV, d1, log_d1;

  // calcul de V 
  V = theta(0)*MatrixXd::Identity(n,n);
  for(int j = 0; j < s; j++) V.noalias() += theta[j+1]*K[j];

  // Calcul de Vi = inverse(V)
  sym_inverse(V,Vi,log_detV,detV,1e-7); // detruit V

  // Calcul de X' Vi X et de P
  ViX.noalias() = Vi * x;
  XViX.noalias() = x.transpose() * ViX;
  sym_inverse(XViX, XViX_i, log_d1, d1, 1e-5);
  P.noalias() = Vi - ViX * XViX_i * ViX.transpose();

  // Py = P * y (en tenant ompte de la symmétrie de P)
  Py.noalias()   =  P.selfadjointView<Lower>() * y;
  double logL = -0.5*(log_detV + log_d1 + Py.dot(y.col(0)));
  return logL;    
}

template<typename T1, typename T2, typename A>
double re_likelihood_nofix(const Eigen::MatrixBase<T1> & y, const std::vector<T2,A> & K, NumericVector theta) {
  int n(y.rows());
  int s(K.size());

  MatrixXd V(n,n), P(n,n), Vi(n,n);
  VectorXd Py(n);

  double log_detV, detV;

  // calcul de V 
  V = theta(0)*MatrixXd::Identity(n,n);
  for(int j = 0; j < s; j++) V.noalias() += theta[j+1]*K[j];

  // Calcul de Vi = inverse(V)
  sym_inverse(V,Vi,log_detV,detV,1e-7); // detruit V

  // Py = P * y = Vi * y 
  Py.noalias()   =  Vi.selfadjointView<Lower>() * y;
  double logL = -0.5*(log_detV + Py.dot(y.col(0)));
  return logL;    
}
#endif

