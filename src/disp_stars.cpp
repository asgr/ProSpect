#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// disp_stars_cpp: Rcpp implementation of the stellar dispersion convolution.
//
// Equivalent to the R loop in SFHfunc()/SFHburst() when disp_stars=TRUE:
//   for(i in seq_along(grid)){
//     z_seq_log = log10(1 + grid[i]*z_disp)
//     new_wave_log = wave_log + z_seq_log
//     new_lum = lum_log - z_seq_log
//     new_lum = 10^approx(x=new_wave_log, y=new_lum, xout=wave_log,
//                         rule=2, yleft=new_lum[1], yright=new_lum[n])$y
//     lum_conv = lum_conv + weights[i] * new_lum
//   }
//   return(lum_conv * res)
//
// Inputs:
//   wave_log  - log10 wavelengths (monotone increasing, length n)
//   lum_log   - log10 luminosity (length n)
//   z_disp    - wavelength-dependent velocity dispersion in z units (length n)
//   grid      - Gaussian quadrature grid points, e.g. seq(-range, range, by=res)
//   weights   - dnorm(grid)
//   res       - grid spacing (multiplied into final result)
//
// [[Rcpp::export(".disp_stars_cpp")]]
NumericVector disp_stars_cpp(NumericVector wave_log,
                              NumericVector lum_log,
                              NumericVector z_disp,
                              NumericVector grid,
                              NumericVector weights,
                              double res) {
  int n  = wave_log.size();
  int ng = grid.size();

  if (lum_log.size() != n)  stop("lum_log must have same length as wave_log");
  if (z_disp.size()  != n)  stop("z_disp must have same length as wave_log");
  if (weights.size() != ng) stop("weights must have same length as grid");

  NumericVector out(n, 0.0);
  // Reusable buffers for the shifted wavelength / luminosity grids
  NumericVector x_src(n);
  NumericVector y_src(n);

  for (int gi = 0; gi < ng; gi++) {
    double g = grid[gi];
    double w = weights[gi];

    // Build shifted grids: x_src = wave_log + log10(1 + g*z_disp)
    //                      y_src = lum_log  - log10(1 + g*z_disp)
    for (int j = 0; j < n; j++) {
      double val = 1.0 + g * z_disp[j];
      if (val <= 0.0) {
        stop("1 + grid[i]*z_disp <= 0: cannot take log10. "
             "Reduce veldisp or the grid range.");
      }
      double lv = std::log10(val);
      x_src[j] = wave_log[j] + lv;
      y_src[j] = lum_log[j]  - lv;
    }

    // Linear interpolation of y_src(x_src) onto wave_log.
    // Uses rule=2: extrapolate with boundary values y_src[0] / y_src[n-1].
    // Two-pointer sweep works because wave_log is monotone increasing and
    // x_src remains monotone for the small shifts typical of veldisp.
    // lo is intentionally reset to 0 for each grid point gi: since wave_log
    // is scanned left-to-right each iteration, lo must start fresh.
    int lo = 0;
    for (int j = 0; j < n; j++) {
      double xq = wave_log[j];
      double y_interp;

      if (xq <= x_src[0]) {
        y_interp = y_src[0];
      } else if (xq >= x_src[n - 1]) {
        y_interp = y_src[n - 1];
      } else {
        // Advance lo so that x_src[lo] <= xq < x_src[lo+1]
        while (lo < n - 2 && x_src[lo + 1] <= xq) {
          lo++;
        }
        double t = (xq - x_src[lo]) / (x_src[lo + 1] - x_src[lo]);
        y_interp = y_src[lo] + t * (y_src[lo + 1] - y_src[lo]);
      }

      out[j] += w * std::pow(10.0, y_interp);
    }
  }

  // Scale by the grid spacing
  for (int j = 0; j < n; j++) {
    out[j] *= res;
  }

  return out;
}
