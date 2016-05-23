# linearFeedbackModel

This R package is work in progress and might not work as expected.

### Description

`lfm()` estimates the first-order linear feedback model in [Blundell, Griffith and Windmeijer (2002)](http://personalpages.manchester.ac.uk/staff/rachel.griffith/PublishedPapers/BlundellGriffithWindmeijer02.pdf) via Generalized Method of Moments.

### Installation
```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("walshc/linearFeedbackModel")
```

### Usage:
Usage is very similar to the `pgmm()` function in package `plm`:

```r
lfm(y ~ lag(y, k = 1) + x | lag(y, k = 2:4) + lag(x, k = 1:4),
    data = data, effect = "individual", model = "onestep")
```
Type `help(lfm)` for more details.

### References
 - Blundell, Richard & Griffith, Rachel & Windmeijer, Frank, 2002. ["Individual effects and dynamics in count data models"](http://personalpages.manchester.ac.uk/staff/rachel.griffith/PublishedPapers/BlundellGriffithWindmeijer02.pdf), *Journal of Econometrics*, Elsevier, vol. 108(1), pages 113-131, May.


