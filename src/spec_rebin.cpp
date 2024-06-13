#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".wave_rebin_cpp")]]
NumericVector wave_rebin_cpp(NumericVector wave, NumericVector wave_bin, NumericVector wave_bin_lo,
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

  return NULL;
}

// [[Rcpp::export(".spec_rebin_cpp")]]
NumericVector spec_rebin_cpp(NumericVector wave_in, NumericVector flux_in, NumericVector wave_out, bool logbin_in = false, bool logbin_out = false) {
  NumericVector flux_out(wave_out.length());
  NumericVector flux_out_weight(wave_out.length());

  NumericVector wave_bin_in(wave_in.length());
  NumericVector wave_bin_in_lo(wave_in.length());
  NumericVector wave_bin_in_hi(wave_in.length());
  NumericVector wave_bin_out(wave_out.length());
  NumericVector wave_bin_out_lo(wave_out.length());
  NumericVector wave_bin_out_hi(wave_out.length());

  int wave_in_N = wave_in.length();
  int wave_out_N = wave_out.length();

  int i, j;
  double temp_weight = 0;

  // set up wave_bin grids

  wave_rebin_cpp(wave_in, wave_bin_in, wave_bin_in_lo, wave_bin_in_hi, logbin_in);

  wave_rebin_cpp(wave_out, wave_bin_out, wave_bin_out_lo, wave_bin_out_hi, logbin_out);

  int jstart = 0;

  for(i = 0; i < wave_out_N ; i++){
    for(j = jstart; j < wave_in_N; j++){
      // Rcpp::Rcout << i << " " << j << " " << " "  << wave_bin_in_lo(j) << " "  << wave_bin_in_hi(j) << " "  << wave_bin_out_lo(i) << " "  << wave_bin_out_hi(i) << "\n";
      if(wave_bin_in_hi(j) <= wave_bin_out_lo(i)){
        jstart++;
      }else if(wave_bin_in_lo(j) >= wave_bin_out_hi(i)){
        j = wave_in_N;
      }else if(wave_bin_in_lo(j) <= wave_bin_out_lo(i)){ // left edge of wave_in lower than wave_out
        if(wave_bin_in_hi(j) < wave_bin_out_hi(i)){ // right edge of wave_in lower than wave_out
          // 2 wave_in bin cross the lower edge of the wave_out bin but not upper, so some flux is contributed
          temp_weight = (wave_bin_in_hi(j) - wave_bin_out_lo(i))/wave_bin_in(j);
          flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
          flux_out_weight(i) = flux_out_weight(i) + temp_weight;
        }else{ // right edge of wave_in higher than wave_out
          // 4 wave_out bin sits fully inside the wave_in bin, so some flux is contributed
          temp_weight = (wave_bin_out_hi(i) - wave_bin_out_lo(i))/wave_bin_in(j);
          flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
          flux_out_weight(i) = flux_out_weight(i) + temp_weight;
        }
      }else{ // left edge of wave_in higher than wave_out
        if(wave_bin_in_hi(j) <= wave_bin_out_hi(i)){  // right edge of wave_in lower than wave_out
          // 1 wave_in bin sits fully inside the wave_out bin, so all flux is contributed
          flux_out(i) = flux_out(i) + flux_in(j);
          flux_out_weight(i) = flux_out_weight(i) + 1;
        }else{ // right edge of wave_in higher than wave_out
          // 3 wave_in bin cross the higher edge of the wave_out bin, so some flux is contributed
          temp_weight = (wave_bin_out_hi(i) - wave_bin_in_lo(j))/wave_bin_in(j);
          flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
          flux_out_weight(i) = flux_out_weight(i) + temp_weight;
        }
      }

      // if(temp_weight < 0 || temp_weight > 1){
      //   Rcpp::Rcout << i << " " << j << " " << temp_weight << " "  << wave_bin_in_lo(j) << " "  << wave_bin_in_hi(j) << " "  << wave_bin_out_lo(i) << " "  << wave_bin_out_hi(i) << "\n";
      // }
    }
    // Rcpp::Rcout << temp_weight << "\n";
    if(flux_out_weight(i) > 0){
      flux_out(i) = flux_out(i) / flux_out_weight(i);
    }
  }

  return flux_out;
}



