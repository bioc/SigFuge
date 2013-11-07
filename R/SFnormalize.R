SFnormalize <- function(data, flag = 1) {

  n <- ncol(data)
  p <- nrow(data)

  # identify low expression samples
  if (length(flag) == 1) {

    # use SigFuge default cutoffs
    if (flag == 1) {
      medcov <- apply(data, 2, function(x) median(x[x > 0]))
      percov <- apply(data, 2, function(x) mean(x > 0))
      flag <- (medcov < 5) | (percov < .1)

    # no low coverage samples
    } else if (flag == 0) {
      flag <- as.logical(rep(0, n))

    # input not recognized
    } else {
      stop("flag must be 0, 1 or logical n-vector.")
    }

  # else check if user input is logical n-vector
  } else if (length(flag) != n || !is.logical(flag)) {
    stop("flag must be 0, 1 or logical n-vector.")      
  }

  if (sum(!flag) < 3) {
    stop("Too many lowly expressed samples to carry out 
          normalization procedure.
          Must have at least 3 samples above threshold.
          More samples are recommended for reasonable 
          inference from p-value.")
  }

  # normalize non-low expression samples
  data.norm <- data[, !flag]
  coverage <- apply(data.norm, 2, sum)
  data.norm <- apply(data.norm, 2, 
                  function(x) x / sum(x) * median(coverage))


  return(list(data.norm = data.norm, flag = flag))

}