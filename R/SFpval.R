SFpval <- function(data, normalize = 1, flag = 1) {
    
    ## normalize data
    if (normalize == 1) {
        data <- SFnormalize(data, flag)$data.norm
    }

    ## calculate sigclust p-value on normalized data
    data <- log10(data+1)
    sc <- sigclust::sigclust(t(data), nsim = 100, nrep = 2, icovest = 2)
    ## nrep=1 appears to produre inconsistent results

    return(sc)

}
