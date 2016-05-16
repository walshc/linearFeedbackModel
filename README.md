# linearFeedbackModel

This R package is work in progress.

### Description

`lfm.R` estimates the first-order linear feedback model in Blundell, Griffith and Windmeijer (2002) via Generalized Method of Moments.

### Usage:
Usage is very similar to the `pgmm()` function in package `plm`:

```r
lfm(formula, data, effect = c("twoways", "individual"), model = c("onestep", "twosteps"))
```

### References
 - Blundell, Richard & Griffith, Rachel & Windmeijer, Frank, 2002. "Individual effects and dynamics in count data models," *Journal of Econometrics*, Elsevier, vol. 108(1), pages 113-131, May.


