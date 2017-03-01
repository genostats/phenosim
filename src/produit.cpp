#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston.h"

using namespace Rcpp;
using namespace RcppParallel;


struct capucine : public Worker {
  // input 
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  double * p;
  double * v;
  const bool dominance;
  //output
  double * vA;

  //le constructeur 
  capucine(uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, double * p, double * v, bool dominance) 
            : data(data), ncol(ncol), true_ncol(true_ncol), nrow(nrow), p(p), v(v), dominance(dominance) {
    vA = new double[ncol];
    std::fill(vA, vA+ncol, 0); 
  }

  //constructeur pour le split
  capucine(capucine & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), nrow(Q.nrow), p(Q.p), v(Q.v), dominance(Q.dominance) {
    vA = new double[ncol];
    std::fill(vA, vA+ncol, 0);
  }

  // destructeur
  ~capucine() {
    delete [] vA;
  }

  //worker
  void operator()(size_t beg, size_t end) {
    double gg[4];
    gg[3] = 0;
    for(size_t i = beg; i < end; i++) {
      if(p[i] == 0 || p[i] == 1) continue; // SNP monomorphe [centré -> 0, pas réduit...]
      if(dominance) {
        gg[0] = p[i]/(1.-p[i]);
        gg[1] = -1.;
        gg[2] = (1.-p[i])/p[i];
      } else {
        double mu_ = 2*p[i];
        double sd_ = sqrt(2*p[i]*(1-p[i]));
        gg[0] = -mu_/sd_;
        gg[1] = (1-mu_)/sd_;
        gg[2] = (2-mu_)/sd_;
      }
      size_t k = 0;
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t x = data[i][j];
        for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
          vA[k++] += v[i]*gg[x&3];
          x >>= 2;
        }
      }
    }
  }

  // join
  void join(const capucine & Q) {
    std::transform(vA, vA + ncol, Q.vA, vA, std::plus<double>());
    // autrement dit : vA += vA.K;
  }

};


//[[Rcpp::export]]
NumericVector vector_product(XPtr<matrix4> p_A, NumericVector p, NumericVector v, bool sparse, bool dominance) {
  int n = p_A->nrow; // nb snps
  int m = p_A->ncol; // nb inds
  if(n != v.size()) stop("Dimensions mismatch");

  if(sparse) {
    LogicalVector w = (v != 0);
    NumericVector v1 = v[w];
    NumericVector p1 = p[w];
    int nb_snps = v1.length();
    uint8_t ** data = new uint8_t * [nb_snps];
    size_t k = 0;
    for(size_t i = 0; i < n; i++) {
      if(v[i] != 0.) {
        data[k++] = p_A->data[i];
      }
    }
    capucine X(data, p_A->ncol, p_A->true_ncol, p_A->nrow, p1.begin(), v1.begin(), dominance);
    parallelReduce(0, nb_snps, X);
    delete [] data; 

    NumericVector R(m);
    std::copy(X.vA, X.vA+m, R.begin());
    return R;
  } else {
    capucine X(p_A->data, p_A->ncol, p_A->true_ncol, p_A->nrow, p.begin(), v.begin(), dominance);
    parallelReduce(0, n, X);

    NumericVector R(m);
    std::copy(X.vA, X.vA+m, R.begin());
    return R;
  }
}

RcppExport SEXP ps_vector_product(SEXP p_ASEXP, SEXP pSEXP, SEXP vSEXP, SEXP sparseS, SEXP domiS) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseS);
    Rcpp::traits::input_parameter< bool >::type domi(domiS);
    rcpp_result_gen = Rcpp::wrap(vector_product(p_A, p, v, sparse, domi));
    return rcpp_result_gen;
END_RCPP
}

