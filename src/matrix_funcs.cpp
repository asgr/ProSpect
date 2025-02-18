#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// the below seems slower than base R crossprod(matrix, vector), but worth a try
// [[Rcpp::export(".colSums_wt_cpp")]]
NumericVector colSums_wt_cpp(NumericMatrix mat, NumericVector vec_wt=1) {
  int nrow = mat.nrow(), ncol = mat.ncol();
  if(vec_wt.length() != nrow && vec_wt.length() != 1){
    throw std::range_error("Weight length error!");
  }
  NumericVector output(ncol);

  if(vec_wt.length() == 1){
    if(vec_wt(0) == 1){
      output = colSums(mat);
    }else{
      for (int j = 0; j < ncol; j++) {
        double sum = 0.0;
        for (int i = 0; i < nrow; i++) {
          sum += mat(i, j)*vec_wt(i);
        }
        output(j) = sum;
      }
    }
  }else{
    for (int j = 0; j < ncol; j++) {
      double sum = 0.0;
      for (int i = 0; i < nrow; i++) {
        sum += mat(i, j)*vec_wt(i);
      }
      output(j) = sum;
    }
  }
  return output;
}

// [[Rcpp::export(".vec_add_cpp")]]
void vec_add_cpp(NumericVector vec_base, NumericVector vec_add=0) {
  int nvec = vec_base.length();
  if(vec_add.length() != nvec && vec_add.length() != 1){
    throw std::range_error("Vectors length error!");
  }
  if(vec_add.length() == 1){
    if(vec_add(0) != 0){
      for (int i = 0; i < nvec; i++) {
        vec_base(i) += vec_add(0);
      }
    }
  }else{
    for (int i = 0; i < nvec; i++) {
      vec_base(i) += vec_add(i);
    }
  }
}

// [[Rcpp::export(".vec_mult_cpp")]]
void vec_mult_cpp(NumericVector vec_base, NumericVector vec_mult=1) {
  int nvec = vec_base.length();
  if(vec_mult.length() != nvec && vec_mult.length() != 1){
    throw std::range_error("Vectors length error!");
  }
  if(vec_mult.length() == 1){
    for (int i = 0; i < nvec; i++) {
      if(vec_mult(0) != 1){
        vec_base(i) *= vec_mult(0);
      }
    }
  }else{
    for (int i = 0; i < nvec; i++) {
      vec_base(i) *= vec_mult(i);
    }
  }
}
