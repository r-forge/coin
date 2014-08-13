### Conditional Expectation and Covariance
setClass("ExpectCovar",
    representation = representation(
        expectation = "numeric",
        covariance  = "matrix",
        dimension   = "integer"
   )
)

### Expectation and Covariance of the influence function
### (+ sum of weights)
setClass("ExpectCovarInfluence",
    contains = "ExpectCovar",
    representation = representation(
        sumweights = "numeric"
    )
)
