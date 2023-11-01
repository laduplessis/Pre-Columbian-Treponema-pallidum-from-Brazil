require(ggplot2)
require(ggsci)
require(ggpubr)
require(ggtree)
require(phytools)
require(phangorn)


rootToTipDistance <- function(tree) {
  # Memory intensive way to get the result
  # mrca  <- findMRCA(tree, tree$tip.label)
  # dists <- setNames(dist.nodes(tree)[mrca, 1:(tree$Nnode+1)], getSeqParts(tree$tip.label,1))
  return(setNames(node.depth.edgelength(tree)[1:(tree$Nnode+1)], tree$tip.label))
}


addRootToTipDistance <- function(tree, metadata) {
  
  dists  <- rootToTipDistance(tree)
  remove <- setdiff(metadata$name, tree$tip.label)
  
  if (length(remove) > 0) {
      result <- metadata[-which(metadata$name %in% remove),]
  } else {
      result <- metadata
  }
  result <- cbind(result, dists[match(result$name, names(dists))])
  colnames(result)[ncol(result)] <- "distance"
  
  return(result)  
}

sciNotation <- function(x, digits = 1) {
  if (length(x) > 1) {
    return(append(sciNotation(x[1]), sciNotation(x[-1])))
  }
  if (!x) return(0)
  exponent <- floor(log10(x))
  base <- round(x / 10^exponent, digits)
  if (base == 1) {
    as.expression(substitute(10^exponent, list(exponent = exponent))) 
  } else {
    as.expression(substitute(base %*% 10^exponent, list(base = base, exponent = exponent)))
  }
}

getTempEstPlot <- function(tempest, xvals="date", yvals="distance", cols="country", shapes=20, title="", pal=NULL, shuffle=FALSE, legend.position="bottom", ...) {
  
  if (shuffle) {
      tempest$distance <- tempest$distance[sample(nrow(tempest))]
  }
  
  title <- paste0(title,": r = ", round(cor(tempest[[xvals]], tempest[[yvals]]),3))
  
  p <- ggscatter(tempest,  x=xvals, y=yvals, title = title, 
                 fill = cols, color=cols, shape=shapes, size=2, palette=pal, # Points color, shape and size
                 #ellipse=TRUE, ellipse.border.remove=TRUE, ellipse.alpha = 0.25,
                 add = "reg.line",  # Add regression line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, cor.coef = FALSE, # Add confidence interval and correlation coefficient
                 ...)
  p <- p + xlab("Sampling year") + ylab("Root-to-tip divergence") + scale_y_continuous(labels = scales::scientific) +
           theme(axis.text=element_text(size=10), axis.title=element_text(size=11), title=element_text(size=10),
                 legend.position=legend.position, legend.text=element_text(size=11), legend.title=element_blank(), legend.box = "vertical") 
  
  return(p)
}



plotParDensity <- function(distr,  MLest, bw='sj', col=pal.dark(cblue), fill=pal.dark(cblue,0.5), add=FALSE, ...) { 
  
  par_qnt     <- quantile(distr, c(0.025, 0.5, 0.975))
  par_density <- density(distr, bw=bw, from=par_qnt[1], to=par_qnt[3])
  
  if (add==FALSE) {
    plot(c(par_density$x), c(par_density$y), type='n', col=col, ...)
  }
  polygon(c(par_density$x[1], par_density$x, par_density$x[length(par_density$x)]), c(0, par_density$y, 0), col=fill, border=NA)
  lines(par_density$x, par_density$y, col=col)
  #abline(v=par_qnt[2], col=col, lwd=1, lty=3)
  abline(v=MLest, col=col, lwd=1, lty=2)
}





checkClockSignal <- function(tempest, nreps=1000, alpha=0.05, plot=TRUE, clades=NULL, ...) {
  
  pal <- get_palette("aaas", 4)
  paltrans <- sapply(substr(get_palette("aaas", 4), 1, 7), paste0,"77")
  
  truth <- cor(tempest$distance, tempest$date)
  
  shuffled  <- vector("numeric", nreps)
  sdates    <- matrix(0, nrow=nrow(tempest), ncol=nreps)
  for (i in 1:nreps) {
      if (is.null(clades)) {
          dates <- tempest$date[sample(nrow(tempest))]
          shuffled[i] <- cor(tempest$distance, dates) 
      } else {
          dates <- vector("numeric", nrow(tempest))
          for (clade in levels(tempest[, clades])) {
              rows         <- which(tempest[,clades] == clade)
              shuffledrows <- rows[sample(length(rows))]
              dates[rows]  <- tempest$date[shuffledrows]
          }
          shuffled[i] <- cor(tempest$distance, dates)
      }
      sdates[,i] <- dates
  }
  
  # Quantiles of null distribution
  critval <- quantile(shuffled, p=c(alpha/2, 1-alpha/2))
  
  # Proportion of replicates with test statistic >= truth
  pval    <- sum(shuffled >= truth)/nreps
  
  if (plot) {
      cordens <- density(shuffled)
      
      args    <- list(...)
      if ("xlim" %in% names(args)) {
          xlim <- args$xlim
          args$xlim <- NULL
      } else {
        xlim=range(c(cordens$x, truth))
      }
      
      if ("ylim" %in% names(args)) {
          ylim <- args$ylim
      } else {
          ymax    <- max(cordens$y)*1.1
          ylim=c(0, ymax)
      }
      
      plot(cordens$x, cordens$y, type='l', col=pal[1], xlim=xlim, ylim=ylim, yaxs='i', ylab="Density", xlab="Pearson correlation coefficient", las=1)
      polygon(cordens$x, cordens$y, col=paltrans[4], border=NA)
      abline(v=critval, col=pal[2], lty=2)  
      points(truth, ymax*0.01, pch="*", col=ifelse(truth <= critval[1] || truth >= critval[2], pal[2], "#000000"), cex=2)
      text(truth, ymax*0.1, paste0("R = ",round(truth,3), "\np = ", round(pval,10)), pos=3, srt=45, col=ifelse(truth <= critval[1] || truth >= critval[2], pal[2], "#000000"), xpd=TRUE)
  }
  
  return( list(truth=truth, shuffled=shuffled, truedates=tempest$date, clades=tempest[, clades], shuffleddates=sdates, pval=pval))
}
