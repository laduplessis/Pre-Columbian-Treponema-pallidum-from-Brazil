require(ggplot2)
require(ggsci)
require(ggtree)
require(beastio)

#' Select an element from a string of tuples separated by `sep`
#' 
#' @param id The string of tuples
#' @param group The index of the tuple to select (use -1 for the last element)
#' @param sep The separator between tuples
#' 
selectFromSeq <- function(id, group, sep="\\|") {
  if (group > 0) {
      return(strsplit(id,split=sep)[[1]][group]) 
  } else {
      parts <- strsplit(id,split=sep)[[1]]
      return(parts[length(parts)])
  }
}

#' Get parts of a list of ordered tuples separated by `sep`
#' 
#' @param labels The list of tuples as strings
#' @param group The element in the tuples to return (use -1 for the last element)
#' @param sep The separator between tuples
#' 
getSeqParts <- function(labels, group=1, sep="\\|") {
  return(sapply(labels, selectFromSeq, group=group, sep=sep))
}


#' Use this function for reading in RAxML bootstrap summary file "RAxML_bipartitions.*"
#' Use read.raxml for reading in RAxML bootstrap summary file "RAxML_bipartitionsBranchLabels.*" instead <-- this is preferred!
read.raxml2 <- function(file, root = FALSE) {
  
  tree.text <- readLines(file)
  phylo <- read.tree(text=tree.text)
  if (root) {
      phylo <- phangorn::midpoint(phylo, node.labels='support')
  }
  bootstrap <- as.numeric(phylo$node.label)
  phylo$node.label <- NULL
  
  bootstrap <- data_frame(node = Ntip(phylo) + 1:phylo$Nnode, 
                          bootstrap = bootstrap)
  new("treedata", file = file, treetext = tree.text, 
      phylo = phylo, data = bootstrap)
}



#' Convert a phylo object to a treedata object (treeio/ggtree)
#' Add node ages to the treedata
phylo2treedata <- function(phylo) {
  
    if (class(phylo)[1] != "phylo") {
        # Quick hack to satisfy treeio
        if ("phylo" %in% class(phylo)) {
            class(phylo) <- "phylo"
        } else {
            stop("Error: Tree not of class phylo")
        }
    }
    intervals <- beastio::getTreeIntervals(phylo)
    
    heights <- tidytree::tibble(node = intervals$node, height=intervals$height, nodetype=intervals$nodetype)
    new("treedata", file="", treetext=write.tree(phylo), 
        phylo=phylo, data=heights)
}




#' Extract metadata from taxon labels. 
#' Function will only work for all fields if each taxon tuple has the same number of fields (i.e. is well-behaved)
#' Augment this data frame with data from a csv file
getMetadata <- function(labels, key=list(accession=1, date=-1)) {
  
  result <- list(taxon=labels)
  for (k in names(key)) {
    result[[k]] <- unname(getSeqParts(labels, key[[k]]))
  }
  return(data.frame(result, stringsAsFactors = FALSE))
}



formatDateRange <- function(x) {
  
  x$date       <- round(x$date, 2)
  x$date_upper <- round(x$date_upper, 2)
  x$date_lower <- round(x$date_lower, 2)
  
  if (x$date_upper != x$date_lower) {
      result <- paste0(x$name, " (", x$date_lower,"-", x$date_upper,")")
  } else {
      result <- paste0(x$name," (",x$date,")")
  }
  return(result)
}


#' Return ggtree object that can be plotted/saved
#' 
#' @param cols List of colours used, in order (foreground/tree, labels, highlights)
#' @param highlighted Subset of sequence ids to highlight
#' 
getTreePlot <- function(tree, metadata,  
                        dateIds=FALSE, plotIds=TRUE, plotNodes=FALSE, plotNodeBars=FALSE, 
                        plotBootstrap=TRUE, plotPosterior=FALSE, 
                        plotScale=TRUE, plotClades=TRUE, plotLegend=TRUE,
                        alignIds=FALSE, scaleWidth=0.005, nodeCols="#69913B", nodePal=scale_color_aaas(), nodeBars="height_0.95_HPD",
                        cols=c("#002147", "#69913B", "#872434", "#4891DC"), highlighted=c(), legendPosition=c(0.05,0.75), ...) {
  
  
  remove   <- setdiff(metadata$name, tree@phylo$tip.label)
  if (length(remove) > 0) {
      metadata <- metadata[-which(metadata$name %in% remove),]
  }
  
  if (dateIds) {
      # Only accession
      #metadata$name <- tree@phylo$tip.label <- as.character(metadata$name)
    
      # Accession and date
      highlighted <- sapply(match(highlighted, metadata$name), function(i) formatDateRange(metadata[i,]))
    
      tiporder <- match(metadata$name, tree@phylo$tip.label)
      metadata$name <- tree@phylo$tip.label[tiporder] <- sapply(1:nrow(metadata), function(i) formatDateRange(metadata[i,]))
  } 
  
  
  # Base tree
  p <- ggtree(tree, ladderize=TRUE, linewidth=0.25, color=cols[1]) %<+% metadata
  
  # Colour tips
  if (plotClades) {
      p <- p + geom_tippoint(aes(color=clade), size=2) + nodePal
  }
  
  if (plotNodeBars) {
      p <- p + geom_range(range=nodeBars, color=cols[4], alpha=0.5, size=0.75, center='height') 
  }
  
  # Add bootstrap on nodes
  if (plotBootstrap) {
      p <- p + geom_point2(aes(subset=!is.na(bootstrap) & bootstrap >= 50 & bootstrap < 75), color=cols[1], size=1.5) + 
               geom_text2(aes(subset=(!is.na(bootstrap) & bootstrap >= 50 & bootstrap < 75), label=bootstrap), vjust=-0.5, hjust=1.25, size=2.5, color=cols[1]) +
               geom_point2(aes(subset=!is.na(bootstrap) & bootstrap >= 75), color=cols[3], size=1.5) + 
               geom_text2(aes(subset=(!is.na(bootstrap) & bootstrap >= 75), label=bootstrap), vjust=-0.5, hjust=1.25, size=2.5, color=cols[3]) 
  }
  
  # Add posterior support on nodes
  if (plotPosterior) {
    p <- p + geom_point2(aes(subset=!is.na(posterior) & posterior >= 0.9), color=cols[3], size=0.5) + 
             geom_text2(aes(subset=(!is.na(posterior) & posterior >= 0.9), label=round(posterior,2)), vjust=-0.5, hjust=1.25, size=1.5, color=cols[3]) 
  }
  
  # Add sequence ids on tips
  if (plotIds) {
      #p <- p + geom_tiplab(aes(subset=!(label %in% highlighted)), size=3.5, hjust=-0.05, align=alignIds, color=cols[2]) 
      #p <- p + geom_tiplab(aes(subset=(label %in% highlighted)),  size=3.5, hjust=-0.05, align=alignIds, color=cols[3], fontface="bold") 
      p <- p + geom_tiplab(size=1.5, hjust=-0.05, align=alignIds, color=cols[2], show.legend=FALSE) 
      
  }
  
  
  # Add internal node labels (useful for finding clades)
  if (plotNodes) {
      p <- p + geom_text2(aes(label=node), vjust=1.5, hjust=1.25, size=2, color=cols[2]) 
  }
  
  
  # Add legend
  if (plotLegend) {
      p <- p + theme(legend.direction = "vertical", legend.position = legendPosition, 
                     legend.justification = c("left","center"), 
                     legend.text = element_text(size = 10), legend.title = element_text(size=10))
    }
  
  
  # Add scale bar
  if (plotScale) {
      p <- p + geom_treescale(y=-0.02*tree@phylo$Nnode, fontsize=3.5, linesize=0.25, offset=-0.03*tree@phylo$Nnode, width=scaleWidth) 
  }
  
  
  return(p)
}


getCladeIds <- function(cladeName, clades) as.character(clades$SequenceID[clades$Cluster == cladeName, drop=TRUE])



# Label clades
annotateClades <- function(treeplot, cladenodes, hilightClade=TRUE, labelClade=TRUE, alpha=0.4, cols="#872434", ...) {
  
  if (length(cols) < length(cladenodes)) {
    cols <- rep(cols, length(cladenodes))
  }
  
  for (i in 1:length(cladenodes)) {
    cladename <- ifelse(is.null(names(cladenodes)), paste("Clade",i), names(cladenodes)[i])
    
    if (hilightClade) {
      treeplot <- treeplot + geom_hilight(node=cladenodes[i], fill="#777777", alpha=alpha)
    }
    
    if (labelClade) {
      treeplot <- treeplot + geom_cladelabel(node=cladenodes[i], label=cladename, color=cols[i], ...)
    }
  }
  
  return(treeplot)
  
}


