## ------------------------------------------------------------------------- ##
# --------------------  read VCFtools ind stats files ----------------------- #

library(dplyr)

# read individual stats files generated using vcftools
# specify data path (dir) and name of file to be loaded (vcf, no file ending)
# will load all het, idepth and imiss files and join into one data frame

# dir <- c("results")
# vcf <- c("SNAPPER_FIL-4")

read.ind.stats <- function(dir, vcf) {
  # read depth stats
  filename <- paste(vcf, ".idepth", sep = "")
  path <- file.path(dir, filename)
  idepth <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  names(idepth) <- c("INDV",
                     "SITES", 
                     "MEAN_DEPTH")
  
  # read missing stats
  filename <- paste(vcf, ".imiss", sep = "")
  path <- file.path(dir, filename)
  imiss <- read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    select(INDV, F_MISS)
  names(imiss) <- c("INDV",
                    "MISS")
  # join stats
  temp <- left_join(imiss, idepth)
  
  # read missing stats
  filename <- paste(vcf, ".het", sep = "")
  path <- file.path(dir, filename)
  het <- read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    select(INDV, `F`)
  names(het) <- c("INDV",
                  "Fis")
  # join stats
  final <- left_join(temp, het)
}

## ------------------------------------------------------------------------- ##



## ------------------------------------------------------------------------- ##
# --------------------  read VCFtools locus stats files ----------------------- #

library(dplyr)

# read locus stats files generated using vcftools
# specify data path (dir) and name of file to be loaded (vcf, no file ending)
# will load all het, idepth and imiss files and join into one data frame

read.loc.stats <- function(dir, vcf) {
  # read depth stats
  filename <- paste(vcf, ".ldepth.mean", sep = "")
  path <- file.path(dir, filename)
  ldepth <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  names(ldepth) <- c("CHR", "POS", 
                     "MEAN_DEPTH", "VAR_DEPTH")
  
  # read missing stats
  filename <- paste(vcf, ".lmiss", sep = "")
  path <- file.path(dir, filename)
  lmiss <- read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
    select(CHR, POS, F_MISS)
  names(lmiss) <- c("CHR", "POS", "MISS")
  # join stats
  temp <- left_join(lmiss, ldepth)
}

## ------------------------------------------------------------------------- ##

# add to ggplot.R script?
# Multiple plot function (http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
