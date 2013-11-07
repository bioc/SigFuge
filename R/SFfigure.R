SFfigure <- function(data, locusname, annot = c(), flip.fig = 1, 
                     label.exon = 1, print.n = 1, data.labels = 0,
                     label.colors = c(), flag = 1, lplots = 2, 
                     log10 = 1, summary.type = "median",
                    savestr = c(), titlestr = c(), pval = 1) {


  # data dims and force data to matrix
  n <- ncol(data)
  p <- nrow(data)
  data <- as.matrix(data)

  if (is.null(colnames(data))) {
    colnames(data) <- paste("sample",1:n)
  }
    

  # parse savestr to separate name and extension (if specified)
  useAlpha <- TRUE
  if (!is.null(savestr)) {

    savestr <- unlist(strsplit(savestr,"[.]"))
    parts <- length(savestr)
    if (parts == 1) {
      savestr <- c(savestr,".pdf")
    } else {
      savestr <- c(paste(savestr[1:(parts-1)], collapse="."),
              paste0(".", tolower(savestr[parts])))
    }

    # .eps/.ps does not support semi-transparency
    if (savestr[2] == ".eps" || savestr[2] == ".ps") {
      warning(paste(".eps/.ps extensions do not support semi-transparency.", 
          "Figures will be plotted using white background.",
          sep="\n"))
      useAlpha <- FALSE
    }
  
  } # end of savestr if statement

  # whether "_plot#" should be added to savestr
  onefig <- (length(lplots) == 1)


  # parse annotations if provided
  exon.annot <- c()
  chr <- ""
  dir <- "+"
  if (class(annot) == "GRanges") {
    annot <- as.data.frame(annot)
  }
  if (is.data.frame(annot)) {
    exon.annot <- rbind(annot$start, annot$end)
    chr <- as.character(unique(annot$seqnames))
    dir <- as.character(unique(annot$strand))   
  } else if (!is.null(annot)) {
    warning("annot input ignored")
  }


  # flip data if on neg strand
  if (dir == "-") {
    if (flip.fig == 1) { 
      data <- matlab::flipud(data)
    }
  } else if (dir != "+") {
    warning(paste("dir must be '+' or '-'.",
          "assuming '+' strand",
          sep="\n"))
    dir <- "+"
  }


  # exon boundaries for plots
  if (!is.null(exon.annot)) {
    exon.start <- exon.annot[1,]
    exon.end <- exon.annot[2,]

    # vector of exon boundary labels
    if (dir == "-" && flip.fig == 1) {
      exon.SE <- rbind(rev(exon.start), rev(exon.end))
      exon.mod <- rbind(c("",rev(exon.start)), c(rev(exon.end),""))
    } else {
      exon.SE <- rbind(exon.start, exon.end)
      exon.mod <- rbind(c("",exon.end), c(exon.start,""))
    }

    # data.frame defining exon boundaries
    len <- c(0, cumsum(exon.SE[2,]-exon.SE[1,]))
    exon.df <- data.frame(p = (1:(length(len)-1))%%2,
                          s = len[-length(len)],
                          e = len[-1])

  } # end of exon boundary if statement


  # default plot title
  if (is.null(titlestr)) {
      titlestr <- paste(locusname, " locus, SigFuge analysis")
  }


  # calculate SigFuge p-value and add to title
  if (pval == 1) {
    pvalue <- SFpval(data,1)@pvalnorm
    titlestr <- paste0(titlestr,", pval = ",signif(pvalue,3))
  }

  
  # define cluster labels by SigFuge clustering or user input
  if (length(data.labels) != n) {
  
    # use default SigFuge labels
    if (!all(data.labels == 0)) {
      warning(paste("data.labels input should be '0' or n-vector.",
          "Using default SigFuge labels.",
          sep="\n"))
    }
    SFout <- SFnormalize(data, flag=flag)
    data.labels <- SFlabels(SFout)
    data.labels <- paste("Cluster", data.labels)
    Klevels <- paste("Cluster", 1:3)
    K <- 3

  # use input cluster labels
  } else {

    Klevels <- unique(data.labels)
    K <- length(Klevels)

  } # end of cluster if statement
  data.labels <- factor(data.labels, levels=Klevels)
  clustSize <- table(data.labels)


  # primary data.frame for ggplot2 input
  data <- as.data.frame(data)
  data$base <- 1:p
  data.df <- reshape::melt.data.frame(data, id.vars=c("base"))
  colnames(data.df) <- c("base", "sample", "coverage")
  data.df$cluster <- rep(data.labels, each=p)
  

  # figure boundaries
  xrng <- range(data.df$base)
  if (log10 == 1) {
    yrng <- range(log10(data.df$coverage+1))
    yr <- yrng[2]-yrng[1]
    y.cut <- 0:ceiling(yrng[2])
    y.cut <- c(y.cut, y.cut+log10(5))
  } else if (log10 == 0) {
    yrng <- range(data.df$coverage)
    yr <- yrng[2]-yrng[1]
    y.cut <- c(0:5)*ceiling(yrng[2]/5)
    y.cut <- c(y.cut, y.cut+5)  
  } else {
    warning("log10 should be either 0 or 1.\n
    Setting log10 to 1.")
    yrng <- range(log10(data.df$coverage+1))
    yr <- yrng[2]-yrng[1]
    y.cut <- 0:ceiling(yrng[2])
    y.cut <- c(y.cut, y.cut+log10(5))
    log10 <- 1
  }

  # main ggplot object
  if (log10 == 0) {
    main.plot <- ggplot(data.df, aes(base,coverage)) + 
                scale_y_continuous(limits=yrng, expand=c(.1,0)) +
                xlab(paste(chr,"genomic positions")) +
                theme(panel.background=element_rect(fill="white"),
                        panel.grid.major=element_line(color="#CCCCCC"),
                        panel.grid.minor=element_line(color="white")) +
                geom_hline(yintercept=0) +
                ylab(expression(paste("read depth")))
  } else {
    main.plot <- ggplot(data.df, aes(base,log10(coverage+1))) + 
                scale_y_continuous(limits=yrng, breaks=y.cut, labels=c(0,10^y.cut[-1]), expand=c(.1,0)) +
                xlab(paste(chr,"genomic positions")) +
                theme(panel.background=element_rect(fill="white"),
                        panel.grid.major=element_line(color="#CCCCCC"),
                        panel.grid.minor=element_line(color="white")) +
                geom_hline(yintercept=0) +
                ylab(expression(paste("read depth (",log[10]," plotting)")))
  }
  
  # add exon boundaries to ggplot object if desired
  if (!is.null(exon.annot)) {
    if (!label.exon) {
        exon.mod[] <- ""
    }   
    main.plot <- main.plot +
                    scale_x_continuous(breaks=len, 
                                 labels=apply(exon.mod,2,paste,collapse=' \n '), 
                                 expand=c(0,0)) +
                    geom_rect(aes(NULL, NULL, xmin=s, xmax=e, fill=factor(p)),
                                    ymin=yrng[1]-yr*.1, ymax=yrng[2]+yr*.3, 
                                    data=exon.df, alpha=(1/5)*useAlpha, 
                                    show_guide=F) +
                    scale_fill_manual(values=c("#FFCC99","#99FFFF")) +
                    theme(axis.text.x=element_text(angle=90,vjust=1/2))
  }



  ## create specified lplots
  if (is.null(savestr)) {
    allplots <- list()
  }
    
  # plot of all curves, no cluster information
  if (any(lplots == 1)) {
  
    tplot <- main.plot + 
                geom_line(aes(color=factor(sample),group=sample), 
                        alpha=(1/2)^useAlpha, size=.4, show_guide=F) + 
                scale_color_hue(l=50, c=100) +
                ggtitle(titlestr)
    if (!is.null(savestr)) {
      fname <- paste0(savestr[1],ifelse(onefig,"","_1"),savestr[2])
      ggsave(filename=fname, plot=tplot, width=10, height=7)
    } else {
      allplots <- c(allplots, list("plot1"=tplot))
    }

  } # end of (lplots == 1) if statement


  # if plots with clustering information are specified...
  if (any(lplots > 1)) {

    # use colors specified by user, if possible
    if (!is.null(label.colors)) {
      if (nrow(label.colors) == K && ncol(label.colors) == 3) {
        mcolors <- rgb(label.colors[,1], label.colors[,2], label.colors[,3])
        label.colors <- pmax(label.colors-.3,0)
        medcolors <- rgb(label.colors[,1], label.colors[,2], label.colors[,3])
        medcolors[clustSize==1] <- mcolors[clustSize==1]
      } else {
        warning(paste("label.colors must be Kx3 RGB matrix.",
            "Using default colors.",
            sep="\n"))
      }
    } else { # use equally spaced colors 
      hues <- seq(15, 375, length=K+1)
      mcolors <- hcl(h=hues, l=50, c=100)[1:K]    
      medcolors <- hcl(h=hues, l=20, c=100)[1:K]
    }
        
    # ggplot object with cluster information
    cluster.plot <- main.plot +
                        geom_line(aes(color=cluster, group=sample), alpha=1, size=.4) +
                        scale_color_manual(name="Clusters", values=mcolors, drop=FALSE) +
                  ggtitle(titlestr)

    # data.frame containing cluster size information
    counts.df <- data.frame(x = .75*p,
                        y1 = (1-.05*1:K)*(yrng[2]),
                        y2 = .9*(yrng[2]),
                        label = paste0("size = ", clustSize),
                        cluster = Klevels)


    # plots containing all clusters in single panel
    if (any(lplots == 2) || any(lplots == 5)) {

      med.mat <- c()
      if (log10 == 0) {
        for(k in 1:K) {
          if (summary.type == "median") {
            med.mat <- cbind(med.mat, apply(as.matrix(data[,data.labels==Klevels[k]]),1,median))
          } else if (summary.type == "mean") {
            med.mat <- cbind(med.mat, apply(as.matrix(data[,data.labels==Klevels[k]]),1,mean))                  
          } else {
            warning("summary.type needs to equal either median or mean.
                Will use default of median")
            summary.type <- "median"
            med.mat <- cbind(med.mat, apply(as.matrix(data[,data.labels==Klevels[k]]),1,median))
          }
        }
        yrng5 <- range(med.mat)
        yr5 <- yrng5[2]-yrng5[1]
        y.cut5 <- 0:ceiling(yrng5[2])
        y.cut5 <- c(y.cut5, y.cut5+5)
      } else {
        for(k in 1:K) {
          if (summary.type == "median") {
            med.mat <- cbind(med.mat, apply(log10(as.matrix(data[,data.labels==Klevels[k]])+1),1,median))
          } else if (summary.type == "mean") {
            med.mat <- cbind(med.mat, apply(log10(as.matrix(data[,data.labels==Klevels[k]])+1),1,mean))                 
          } else {
            warning("summary.type needs to equal either median or mean.
                Will use default of median")
            summary.type <- "median"
            med.mat <- cbind(med.mat, apply(log10(as.matrix(data[,data.labels==Klevels[k]])+1),1,median))
          }
        }
        yrng5 <- range(med.mat)
        yr5 <- yrng5[2]-yrng5[1]
        y.cut5 <- 0:ceiling(yrng5[2])
        y.cut5 <- c(y.cut5, y.cut5+log10(5))
      }
      
      med.df <- data.frame(base=1:p, med.mat)
      if (K==1) {
        names(med.df) <- c('base','X1')
      }

      # plot of all clusters in single panel with all curves
      if (any(lplots == 2)) {

        tplot <- cluster.plot
        if (print.n) {
          tplot <- tplot + 
                   geom_text(data=counts.df, aes(x,y1,label=paste(cluster,label)), 
                                    size=4.5, hjust=0)
        }
        for (k in 1:K) {
          tplot <- tplot +
                      geom_line(data=med.df, aes_string(x="base",y=paste0("X",k)),
                                  size=.6, linejoin="mitre", color=medcolors[k])
        }
        if (!is.null(savestr)) {
          fname <- paste0(savestr[1],ifelse(onefig,"","_2"),savestr[2])
          ggsave(filename=fname, plot=tplot, width=10, height=7)      
        } else {
          allplots <- c(allplots, list("plot2"=tplot))
        }

      } # end of (lplots == 2) if statement


      # plot of all clusters in single panel with only median curves
      if (any(lplots == 5)) {

        counts.df5 <- data.frame(x = .75*p,
                y1 = (1-.05*1:K)*(yrng5[2]),
                y2 = .9*(yrng5[2]),
                label = paste0("size = ", clustSize),
                cluster = Klevels)

        tplot <- main.plot
        if (print.n) {
          tplot <- tplot + 
                      geom_text(data=counts.df5, aes(x,y1,label=paste(cluster,label)), 
                                  size=4.5, hjust=0)
        }
        med.df <- reshape::melt.data.frame(med.df, id.vars="base", variable_name="Clusters")
        
        if (log10 == 0) {
          tplot <- tplot + 
                  geom_line(data=med.df, aes(x=base,y=value,color=Clusters), 
                              size=.6, linejoin="mitre") +
                  scale_color_manual(name="Clusters", values=medcolors, drop=FALSE, labels=Klevels) +
                  ggtitle(paste0(titlestr,", cluster ", summary.type, "s")) +
                  scale_y_continuous(limits=yrng5, expand=c(.1,0))
        } else {
          tplot <- tplot + 
                  geom_line(data=med.df, aes(x=base,y=value,color=Clusters), 
                              size=.6, linejoin="mitre") +
                  scale_color_manual(name="Clusters", values=medcolors, drop=FALSE, labels=Klevels) +
                  ggtitle(paste0(titlestr,", cluster ", summary.type, "s")) +
                  scale_y_continuous(limits=yrng5, breaks=y.cut5, labels=c(0,10^y.cut5[-1]), expand=c(.1,0))
        }
        if (!is.null(savestr)) {
          fname <- paste0(savestr[1],ifelse(onefig,"","_5"),savestr[2])
          ggsave(filename=fname, plot=tplot, width=10, height=7)      
        } else {
          allplots <- c(allplots, list("plot5"=tplot))
        }

      } # end of (lplots == 5) if statement

    } # end of (lplots == 2,5) if statement


    # plot of clusters in separate panels
    if (any(lplots == 3)) {

      tplot <- cluster.plot + 
              facet_grid(cluster~., drop=FALSE) + 
              stat_summary(fun.y=summary.type, geom="line", size=1)
      if (print.n) {
        tplot <- tplot + 
            geom_text(data=counts.df, aes(x,y2,label=paste(cluster,label)), size=3, hjust=0)
      }
      # save plot or return figures
      if (!is.null(savestr)) {
        fname <- paste0(savestr[1],ifelse(onefig,"","_3"),savestr[2])
        ggsave(filename=fname, plot=tplot, width=7, height=7)      
      } else {
        allplots <- c(allplots, list("plot3"=tplot))
      }

    } # end of (lplots == 3) if statement


    # plot each curve in a separate sheet, colored by cluster
    if (any(lplots == 4)) {

      samples <- colnames(data)

      # save plot if savestr is provided
      if (is.null(savestr)) {
        warning("Plot option 4 can only be saved to file. Must provide savestr.")
      } else {
        fname <- paste0(savestr[1],ifelse(onefig,"","_4"),savestr[2])
        pdf(file=fname, onefile=TRUE, width=10, height=7)
        for (i in 1:n) {        
          print(cluster.plot %+% subset(data.df, sample==samples[i]) + 
          ggtitle(samples[i]))
        }
        dev.off()
      }

    } # end of (lplots == 4) if statement

  } # end of (lplots > 1) if statement


  # return list of plots if no savestr specified
  if (is.null(savestr)) {
    return(allplots)
  }

}
