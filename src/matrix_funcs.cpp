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
          sum += mat(i, j)*vec_wt(0);
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

// [[Rcpp::export(".mat_vec_mult_col_cpp")]]
NumericMatrix mat_vec_mult_col(NumericMatrix mat, NumericVector vec=1, int col_lim=0) {
  int nrow = mat.nrow(), ncol = mat.ncol();

  NumericMatrix output(nrow, ncol);

  if(vec.length() != nrow && vec.length() != 1){
    throw std::range_error("Vec length error!");
  }

  if(vec.length() == 1){
    if(vec(0) == 1){
      for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++) {
          output(i, j) = mat(i, j);
        }
      }
    }else{
      for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++) {
          if(j < col_lim){
            output(i, j) = mat(i, j) * vec(0);
          }else{
            output(i, j) = mat(i, j);
          }
        }
      }
    }
  }else{
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < nrow; i++) {
        if(j < col_lim){
          output(i, j) = mat(i, j) * vec(i);
        }else{
          output(i, j) = mat(i, j);
        }
      }
    }
  }
  return(output);
}

// NumericMatrix mat_vec_mult_row(NumericMatrix mat, NumericVector vec=1, int row_lim=0) {
//   int nrow = mat.nrow(), ncol = mat.ncol();
//
//   NumericMatrix output(nrow, ncol);
//
//   if(vec.length() != ncol && vec.length() != 1){
//     throw std::range_error("Vec length error!");
//   }
//
//   if(vec.length() == 1){
//     if(vec(0) == 1){
//       for (int j = 0; j < ncol; j++) {
//         for (int i = 0; i < nrow; i++) {
//           output(i, j) = mat(i, j);
//         }
//       }
//     } else {
//       for (int j = 0; j < ncol; j++) {
//         for (int i = 0; i < row_lim; i++) {
//           output(i, j) = mat(i, j) * vec(0);
//         }
//         for (int i = row_lim; i < nrow; i++) {
//           output(i, j) = mat(i, j);
//         }
//       }
//     }
//   } else {
//     for (int j = 0; j < ncol; j++) {
//       for (int i = 0; i < row_lim; i++) {
//         output(i, j) = mat(i, j) * vec(i);
//       }
//       for (int i = row_lim; i < nrow; i++) {
//         output(i, j) = mat(i, j);
//       }
//     }
//   }
//   return(output);
// }

// [[Rcpp::export(".mat_vec_mult_row_cpp")]]
void mat_vec_mult_row(NumericMatrix mat, NumericVector vec=1, int row_lim=0) {
  int ncol = mat.ncol();

  if(vec.length() != ncol && vec.length() != 1){
    throw std::range_error("Vec length error!");
  }

  if(vec.length() == 1){
    if(vec(0) == 1){
      // No need to modify as vec(0) == 1 means no change
    } else {
      for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < row_lim; i++) {
          mat(i, j) = mat(i, j) * vec(0);
        }
      }
    }
  } else {
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < row_lim; i++) {
        mat(i, j) = mat(i, j) * vec(i);
      }
    }
  }
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
