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
  
  // set up wave_bin_in grids
  
  for(i = 1; i < wave_in_N - 1; i++){
    if(logbin_in){
      wave_bin_in_lo(i) = sqrt(wave_in(i) * wave_in(i - 1));
      wave_bin_in_hi(i) = sqrt(wave_in(i) * wave_in(i + 1));
    }else{
      wave_bin_in_lo(i) = (wave_in(i) + wave_in(i - 1))/2;
      wave_bin_in_hi(i) = (wave_in(i) + wave_in(i + 1))/2;
    }
    wave_bin_in(i) = wave_bin_in_hi(i) - wave_bin_in_lo(i);
  }
  // fix the extremes to be sensible values
  if(logbin_in){
    wave_bin_in(0) = pow(wave_in(0),2) / wave_in(1);
    wave_bin_in(wave_in_N - 1) = pow(wave_in(wave_in_N - 1),2) / wave_in(wave_in_N - 2);
  }else{
    wave_bin_in(0) = wave_in(1) - wave_in(0);
    wave_bin_in(wave_in_N - 1) = wave_in(wave_in_N - 1) - wave_in(wave_in_N - 2);
  }
  
  wave_bin_in_lo(0) = wave_bin_in_lo(1) - wave_bin_in(0);
  wave_bin_in_lo(wave_in_N - 1) = wave_bin_in_lo(wave_in_N - 2) + wave_bin_in(wave_in_N - 1);
  
  wave_bin_in_hi(0) = wave_bin_in_hi(1) - wave_bin_in(0);
  wave_bin_in_hi(wave_in_N - 1) = wave_bin_in_hi(wave_in_N - 2) + wave_bin_in(wave_in_N - 1);
  
  // set up wave_bin_out grids
  
  for(i = 1; i < wave_out_N - 1; i++){
    if(logbin_out){
      wave_bin_out_lo(i) = sqrt(wave_out(i) * wave_out(i - 1));
      wave_bin_out_hi(i) = sqrt(wave_out(i) * wave_out(i + 1));
    }else{
      wave_bin_out_lo(i) = (wave_out(i) + wave_out(i - 1))/2;
      wave_bin_out_hi(i) = (wave_out(i) + wave_out(i + 1))/2;
    }
    wave_bin_out(i) = wave_bin_out_hi(i) - wave_bin_out_lo(i);
  }
  // fix the extremes to be sensible values
  if(logbin_out){
    wave_bin_out(0) = pow(wave_out(0),2) / wave_out(1);
    wave_bin_out(wave_out_N - 1) = pow(wave_out(wave_out_N - 1),2) / wave_out(wave_out_N - 2);
  }else{
    wave_bin_out(0) = wave_out(1) - wave_out(0);
    wave_bin_out(wave_out_N - 1) = wave_out(wave_out_N - 1) - wave_out(wave_out_N - 2);
  }
  
  
  wave_bin_out_lo(0) = wave_bin_out_lo(1) - wave_bin_out(0);
  wave_bin_out_lo(wave_out_N - 1) = wave_bin_out_lo(wave_out_N - 2) + wave_bin_out(wave_out_N - 1);
  
  wave_bin_out_hi(0) = wave_bin_out_hi(1) - wave_bin_out(0);
  wave_bin_out_hi(wave_out_N - 1) = wave_bin_out_hi(wave_out_N - 2) + wave_bin_out(wave_out_N - 1);
  
  int jstart = 0;
  
  for(i = 0; i < wave_out_N ; i++){
    // Rcpp::Rcout << "i  " << i << "\n";
    for(j = jstart; j < wave_in_N; j++){
      // Rcpp::Rcout << "j  " << j << "\n";
      if(wave_bin_in_hi(j) <= wave_bin_out_lo(i)){
        jstart++;
      }else if(wave_bin_in_lo(j) >= wave_bin_out_hi(i)){
        j = wave_in_N;
      }else if(wave_bin_in_lo(j) >= wave_bin_out_lo(i) && wave_bin_in_hi(j) <= wave_bin_out_hi(i)){
        // 1 wave_in bin sits fully inside the wave_out bin, so all flux is contributed
        flux_out(i) = flux_out(i) + flux_in(j);
        flux_out_weight(i) = flux_out_weight(i) + 1;
      }else if(wave_bin_in_lo(j) <= wave_bin_out_lo(i) && wave_bin_in_hi(j) <= wave_bin_out_hi(i)){
        // 2 wave_in bin cross the lower edge of the wave_out bin, so some flux is contributed
        temp_weight = (wave_bin_in_hi(j) - wave_bin_out_lo(i))/wave_bin_in(j);
        flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
        flux_out_weight(i) = flux_out_weight(i) + temp_weight;
      }else if(wave_bin_in_lo(j) <= wave_bin_out_hi(i) && wave_bin_in_hi(j) >= wave_bin_out_hi(i)){
        // 3 wave_in bin cross the higher edge of the wave_out bin, so some flux is contributed
        temp_weight = (wave_bin_out_hi(i) - wave_bin_in_lo(j))/wave_bin_in(j);
        flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
        flux_out_weight(i) = flux_out_weight(i) + temp_weight;
      }else if(wave_bin_in_lo(j) <= wave_bin_out_lo(i) && wave_bin_in_hi(j) >= wave_bin_out_hi(i)){
        // 4 wave_out bin sits fully inside the wave_in bin, so some flux is contributed
        temp_weight = (wave_bin_out_hi(i) - wave_bin_out_lo(i))/wave_bin_in(j);
        flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
        flux_out_weight(i) = flux_out_weight(i) + temp_weight;
      }else if(wave_bin_in_lo(j) >= wave_bin_out_lo(i) && wave_bin_in_lo(j) <= wave_bin_out_hi(i)){
        // 5 wave_out bin cross the lower edge of the wave_in bin, so some flux is contributed
        temp_weight = (wave_bin_out_hi(i) - wave_bin_in_lo(j))/wave_bin_in(j);
        flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
        flux_out_weight(i) = flux_out_weight(i) + temp_weight;
      }else if(wave_bin_in_lo(j) >= wave_bin_out_lo(i) && wave_bin_in_lo(j) <= wave_bin_out_hi(i)){
        // 6 wave_out bin cross the higher edge of the wave_in bin, so some flux is contributed
        temp_weight = (wave_bin_in_hi(j) - wave_bin_out_lo(i))/wave_bin_in(j);
        flux_out(i) = flux_out(i) + flux_in(j)*temp_weight;
        flux_out_weight(i) = flux_out_weight(i) + temp_weight;
      }else{
        Rcpp::Rcout << i << " " << j << "\n";
      }
      if(temp_weight < 0 || temp_weight > 1){
        Rcpp::Rcout << i << " " << j << " " << temp_weight << " "  << wave_bin_in_lo(j) << " "  << wave_bin_in_hi(j) << " "  << wave_bin_out_lo(i) << " "  << wave_bin_out_hi(i) << "\n";
      }
    }
    flux_out(i) = flux_out(i) / flux_out_weight(i);
  }
  
  return flux_out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

# /*** R
# wave_in = seq(1e3,9e3, by=5)
# flux_in = rnorm(length(wave_in), mean=seq(100,200,len=length(wave_in)))
# wave_out = seq(2e3,6e3, by=20)
# flux_out = spec_rebin_cpp(wave_in, flux_in, wave_out)
# 
# #yay
# magplot(wave_in, flux_in, type='l')
# lines(wave_out, flux_out, col='red')
# 
# wave_out2 = seq(2e3,6e3, by=2)
# flux_out2 = spec_rebin_cpp(wave_in, flux_in, wave_out2)
# 
# #yay
# magplot(wave_in, flux_in, type='l')
# lines(wave_out2, flux_out2, col='red')
# */
