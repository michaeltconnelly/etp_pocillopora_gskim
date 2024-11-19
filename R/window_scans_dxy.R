# Window scans: Fst, Dxy
# Calculate Dxy in windows that correspond to the ANGSD Fst and Theta windows
library("readr")
library("ggplot2")

setwd("/scratch/nmnh_corals/connellym/projects/etp_pocillopora_gskim")

### ---------------- 
chrs_abbv <- list(
  "OZ034731.1" = "CHR1",
  "OZ034732.1" = "CHR2",
  "OZ034733.1" = "CHR3",
  "OZ034734.1" = "CHR4",
  "OZ034735.1" = "CHR5",
  "OZ034736.1" = "CHR6",
  "OZ034737.1" = "CHR7",
  "OZ034738.1" = "CHR8",
  "OZ034739.1" = "CHR9",
  "OZ034740.1" = "CHR10",
  "OZ034741.1" = "CHR11",
  "OZ034742.1" = "CHR12",
  "OZ034743.1" = "CHR13",
  "OZ034744.1" = "CHR14"
)

### Import Fst and Dxy files -------------
fst_files <- list.files("./outputs/angsd/fst/50kb", full.names = T)
dxy_files <- list.files("./outputs/angsd/dxy", full.names = T)

# Fst file is non-overlapping 50000 bp window size, 50000 bp steps 
angsd_fst_colnames <- c("REGION", "CHR", "wincenter", "Nsites", "Fst")
angsd_fst_cols <- cols("REGION" = col_character(),
                       "CHR" = col_character(),
                       "wincenter" = col_double(), # midPos
                       "Nsites" = col_double(),
                       "Fst" = col_double())
# Dxy file is all sites
dxy_colnames <- c("CHR", "position", "Dxy")

for (i in 1:length(fst_files)){
  
fst_file <- fst_files[i]
fst_df <- read_delim(fst_file, delim = "\t", skip = 1, col_names = angsd_fst_colnames, col_types = angsd_fst_cols)
#
dxy_file <- dxy_files[i]
dxy_df <- read_delim(dxy_file, delim = "\t", skip = 1, col_names = dxy_colnames)

### -------------
dxy.out <- c(rep(-999, nrow(fst_df)))
wsize <- 50000

### -------------
#pos = window center
#+1 is to make sure windows aren't touching
#this code takes a long time so I would suggest making a tester file but more importantly what does it do
#so the fst values are calc on variant sites, the dxy and pi values are calc on all sites including invariant so you have to average over the same windows
#which is an absolute mess
#so this code averages your dxy values over the same windows as the fst values 
#line37 adds the number of dxy values then divides by that value within that window (dxy values/number of sites)

for (j in 1:nrow(fst_df)){
  
# setup windows from Fst data
  scaf <- fst_df[j,2]$CHR
  pos <- fst_df[j,3]$wincenter
  start_pos <- (pos - (wsize/2))
  end_pos <- pos + (wsize/2) - 1
  
# calculate average Dxy across windows
  tdf <- dxy_df[which(dxy_df$CHR == scaf & dxy_df$position >= start_pos & dxy_df$position <= end_pos), ] # & dxy$Dxy > 0 or 2e-6
  tdf.len <- nrow(tdf)
  dxy.out[i] <- (sum(tdf$Dxy))/tdf.len
  
}

#now add the dxy.out to the fst file:
fst_dxy_df <- fst_df
fst_dxy_df$Dxy <- dxy.out

# write dxy output table: REGION, CHR, wincenter, NSites / NSitesDxy, dxy
comp <- str_extract(fst_files[i], "/(p.*_p.*)") %>% str_remove(., "_50kb.*$") %>% str_remove(., "^/")
filename <- paste0("./outputs/angsd/", comp, "_fst_dxy_50kb.txt") 

write_delim(fst_dxy_df, file = filename, delim = "\t")

# make demo window plot

dxy_plot <- ggplot(fst_dxy_df, aes(wincenter, Dxy)) + geom_point() + facet_grid(CHR ~ .)

fst_plot <- ggplot(fst_dxy_df, aes(wincenter, Fst)) + geom_point() + facet_grid(CHR ~ .)

double_window_plot <- fst_plot / dxy_plot 

plotfilename <- filename <- paste0("./outputs/angsd/", comp, "_fst_dxy_50kb_window.pdf")

pdf(file = plotfilename)
print(double_window_plot)
dev.off()

}