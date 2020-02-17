# owl (development version)

## Major changes

* The `standardize_features` argument in `owl()` has been deprecated in favor
  of two new arguments: `scale` and `center` in order to provide
  fine-grained control of normalization. 
* Just-in-time centering (for sparse predictor matrices) is no longer supported.
  As a result, `center` cannot be set to false when the predictor matrix
  is sparse in the call to `owl()`. 
* The FISTA solver has been replaced with an ADMM solver for
  OLS (`family = "gaussian"`). Two new arguments were added to 
  control stopping criterion for the ADMM solver: `tol_rel` and `tol_abs`.
  
## Minor changes

* `print.Owl()` no longer prints the regularization path when called.

## Bug fixes

* Lambda sequences when `lambda = "gaussian"` are now computed properly
  when the number of features is larger than the number of observations.
  
# owl 0.1.1

## Minor changes

* A reference to Bogdan et al. (2015) has been added to the description of the 
  package.

## Bug fixes

* The URL in the readme pointing to the code of conduct is now absolute.

# owl 0.1.0

* First release of owl.
