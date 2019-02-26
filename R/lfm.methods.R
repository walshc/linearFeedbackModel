# This hides the large matrices from printing when the function is called:
print.lfm <- function(x) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

# This shows the table when summary(lfm) is called (not finished yet):
summary.lfm <- function(x) {
  cat("Call:\n")
  x$call
  cat("\nCoefficients:\n")
  lmtest::coeftest(x, vcov = x$vcov)
}
