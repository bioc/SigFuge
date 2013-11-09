SFpval <- function(data, normalize = 1, flag = 1) {

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> change indents to 2 spaces
  # normalize data
  if (normalize == 1) {
    data <- SFnormalize(data, flag)$data.norm
  }

  # calculate sigclust p-value on normalized data
  data <- log10(data+1)
  sc <- sigclust::sigclust(t(data), nsim = 100, nrep = 2, icovest = 2)
      # nrep=1 appears to produre inconsistent results

  return(sc)
<<<<<<< HEAD
=======
    # normalize data
    if (normalize == 1) {
        data <- SFnormalize(data, flag)$data.norm
    }

    # calculate sigclust p-value on normalized data
    data <- log10(data+1)
    sc <- sigclust::sigclust(t(data), nsim=100, nrep=2, icovest=2)
        # nrep=1 appears to produre inconsistent results
  
    return(sc)
>>>>>>> Correct orientation of data input to sigclust()
=======
>>>>>>> change indents to 2 spaces

}