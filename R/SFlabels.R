SFlabels <- function(normData) {
    
    ## obtain 2-means clustering labels from transformed data
    km <- kmeans(t(log10(normData$data.norm + 1)),
                 2, nstart=100)$cluster

    ## sort cluster labels based on first entry for consistency
    if (km[1] == 2) {
        km <- 3 - km
    }

    ## combine 2-means labels with low expression flag
    labels <- rep(1, length(normData$flag))
    labels[which(!normData$flag)] <- km + 1

    return(labels)
}
