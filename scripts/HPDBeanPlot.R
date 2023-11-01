HPDBeanPlot <- function(mcmc.trace, names=NULL, bw = "sj", kernel = "gaussian", fill = "blue", border = "black", medcol = "red", 
                        side="both", maxwidth = 0.5, lwd = c(1,1,NA), add = FALSE, axes = TRUE, ...) {

    if (is.null(names)) {
        names <- colnames(mcmc.trace)
    }
    n <- length(names)
    
    # Actually calculate the statistics
    hpdmed      <- getHPDMedian(mcmc.trace)
    par_density <- lapply(1:ncol(mcmc.trace), function(i) density(mcmc.trace[,i], from=hpdmed[i,1], to=hpdmed[i,3]))
    
    maxdens <- max(sapply(par_density, function(x) max(x$y)))
    scaling <- maxwidth/maxdens
    
    # Draw plot axes
    if (add==FALSE) {
        plot(1,type='n', xlim=c(1-maxwidth,n+maxwidth), ylim=range(hpdmed), xaxt='n', yaxt='n', ...)
    }
    
    if (axes) {
        axis(1, at = 1:n, labels = names)
        axis(2)
    }
    

    # Plot beans
    left  <- side == "both" | side == "left"
    right <- side == "both" | side == "right"
    for (i in 1:n) {
      
        # Left bean
        if (left) {
            polygon(c(i, i-par_density[[i]]$y*scaling, i), c(par_density[[i]]$x[1], par_density[[i]]$x, par_density[[i]]$x[length(par_density[[i]]$x)]), col=fill, border=NA)
            lines(i-par_density[[i]]$y*scaling, par_density[[i]]$x, col=border, lwd=lwd[1])
        }
      
        # Right bean
        if (right) {
              polygon(c(i, i+par_density[[i]]$y*scaling, i), c(par_density[[i]]$x[1], par_density[[i]]$x, par_density[[i]]$x[length(par_density[[i]]$x)]), col=fill, border=NA)
              lines(i+par_density[[i]]$y*scaling, par_density[[i]]$x, col=border, lwd=lwd[1])
        }
    }
    
    if (left) {
        segments(1:n - maxwidth*0.75, hpdmed[,1], 1:n, hpdmed[,1], col=border, lwd=lwd[3])
        segments(1:n - maxwidth,      hpdmed[,2], 1:n, hpdmed[,2], col=medcol, lwd=lwd[2])
        segments(1:n - maxwidth*0.75, hpdmed[,3], 1:n, hpdmed[,3], col=border, lwd=lwd[3])
    }
    
    if (right) {
        segments(1:n + maxwidth*0.75, hpdmed[,1], 1:n, hpdmed[,1], col=border, lwd=lwd[3])
        segments(1:n + maxwidth,      hpdmed[,2], 1:n, hpdmed[,2], col=medcol, lwd=lwd[2])
        segments(1:n + maxwidth*0.75, hpdmed[,3], 1:n, hpdmed[,3], col=border, lwd=lwd[3])
    }
    
}
  