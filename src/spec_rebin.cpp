#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".wave_rebin_cpp")]]
void wave_rebin_cpp(NumericVector wave, NumericVector wave_bin, NumericVector wave_bin_lo,
                             NumericVector wave_bin_hi, bool logbin = false){
  int i;
  int wave_N = wave.length();

  // set up wave_bin grids

  for(i = 1; i < wave_N - 1; i++){
    if(logbin){
      wave_bin_lo(i) = sqrt(wave(i) * wave(i - 1));
      wave_bin_hi(i) = sqrt(wave(i) * wave(i + 1));
    }else{
      wave_bin_lo(i) = (wave(i) + wave(i - 1))/2;
      wave_bin_hi(i) = (wave(i) + wave(i + 1))/2;
    }
  }

  // fix the extremes to be sensible values
  if(logbin){
    wave_bin_lo(0) = pow(wave(0),2) / wave_bin_lo(1);
    wave_bin_hi(wave_N - 1) = pow(wave(wave_N - 1),2) / wave_bin_lo(wave_N - 2);
  }else{
    wave_bin_lo(0) = 2*wave(0) - wave_bin_lo(1);
    wave_bin_hi(wave_N - 1) = 2*wave(wave_N - 1) - wave_bin_lo(wave_N - 2);
  }

  wave_bin_lo(wave_N - 1) = wave_bin_hi(wave_N - 2);
  wave_bin_hi(0) = wave_bin_lo(1);

  for(i = 0; i < wave_N; i++){
    wave_bin(i) = wave_bin_hi(i) - wave_bin_lo(i);
  }
}

// [[Rcpp::export(".spec_rebin_cpp")]]
NumericVector spec_rebin_cpp(NumericVector wave_in,
                                   NumericVector flux_in,
                                   NumericVector wave_out,
                                   Nullable<NumericVector> invar_in = R_NilValue,  // <-- new input: variance
                                   bool logbin_in = false, bool logbin_out = false) {
  int wave_in_N = wave_in.length();
  int wave_out_N = wave_out.length();

  if (wave_in_N < 2 || wave_out_N < 2) stop("wave_in and wave_out must have at least 2 elements");
  if (flux_in.size() != wave_in_N) stop("flux_in must have same length as wave_in");

  // If invar_in is provided, extract it
  bool use_invar = invar_in.isNotNull();
  if (use_invar) {
    NumericVector invar_vec;
    invar_vec = as<NumericVector>(invar_in); //need to copy into the new structure so methods work
    if (invar_vec.size() != wave_in_N) stop("invar_in must have same length as wave_in");
  }


  double wave_in_min = wave_in(0);
  double wave_in_max = wave_in(wave_in_N - 1);
  double wave_out_min = wave_out(0);
  double wave_out_max = wave_out(wave_out_N - 1);

  NumericVector flux_out;

  if(use_invar){
    flux_out = NumericVector(wave_out_N * 2);
  }else{
    flux_out = NumericVector(wave_out_N);
  }

  if (wave_out_min > wave_in_max || wave_out_max < wave_in_min) {
    return flux_out; // no overlap
  }

  NumericVector wave_bin_in(wave_in_N);
  NumericVector wave_bin_in_lo(wave_in_N);
  NumericVector wave_bin_in_hi(wave_in_N);
  NumericVector wave_bin_out(wave_out_N);
  NumericVector wave_bin_out_lo(wave_out_N);
  NumericVector wave_bin_out_hi(wave_out_N);

  wave_rebin_cpp(wave_in, wave_bin_in, wave_bin_in_lo, wave_bin_in_hi, logbin_in);
  wave_rebin_cpp(wave_out, wave_bin_out, wave_bin_out_lo, wave_bin_out_hi, logbin_out);

  int jstart = 0;
  double bin_weight = 0.0;
  double total_weight = 0.0;

  for (int i = 0; i < wave_out_N; i++) {
    if (wave_bin_out_hi(i) < wave_in_min) continue;
    if (wave_bin_out_lo(i) > wave_in_max) break;
    total_weight = 0.0;

    for (int j = jstart; j < wave_in_N; j++) {
      if (wave_bin_in_hi(j) <= wave_bin_out_lo(i)) { jstart++; continue; }
      if (wave_bin_in_lo(j) >= wave_bin_out_hi(i)) break;

      // Compute overlap fraction
      if (wave_bin_in_lo(j) <= wave_bin_out_lo(i)) {
        if (wave_bin_in_hi(j) < wave_bin_out_hi(i)) {
          bin_weight = (wave_bin_in_hi(j) - wave_bin_out_lo(i)) / wave_bin_in(j);
        } else {
          bin_weight = (wave_bin_out_hi(i) - wave_bin_out_lo(i)) / wave_bin_in(j);
        }
      } else {
        if (wave_bin_in_hi(j) <= wave_bin_out_hi(i)) {
          bin_weight = 1.0;
        } else {
          bin_weight = (wave_bin_out_hi(i) - wave_bin_in_lo(j)) / wave_bin_in(j);
        }
      }

      if(use_invar){
        flux_out(i) += flux_in(j) * bin_weight * invar_vec(j);
        total_weight += bin_weight * invar_vec(j);
      }else{
        flux_out(i) += flux_in(j) * bin_weight;
        total_weight += bin_weight;
      }
    }

    if (total_weight > 0) {
      flux_out(i) /= total_weight;
    }

    if(use_invar){
      flux_out(i + wave_out_N) = total_weight;
    }
  }

  return flux_out;
}
