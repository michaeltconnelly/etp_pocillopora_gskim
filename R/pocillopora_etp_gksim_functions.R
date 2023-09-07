#pocillopora_etp_gksim_functions.R
#author: "Mike Connelly"
#date: "05/11/2023"

# project-wide colors and theme
region_colors <- c("gold1", "steelblue", "darkorange", "darkturquoise")
# 
site_colors <- c("gold1", "blue", "darkblue", "steelblue", "lightblue", "darkorange3","darkorange", "orange", "orangered", "cyan", "cyan3", "cyan4", "turquoise", "lightgreen")
names(site_colors)
#
spp_colors <- c("purple", "darkorange1", "gold", "turquoise2", "turquoise4")

theme_set(theme_bw())

# notin operator
`%notin%` <- Negate(`%in%`)

# PCA visualization function
###------------------------------------------------------------------------------------
gg_pcangsd <- function(e, samples, pc = c(1,2)) {
  #
  pcs <- as.data.frame(e$vectors) %>% cbind(samples)
  #
  pcvar <- round(e$values, 2)
  #
  pca_plot <- pcs %>% ggplot(aes(x = pcs[,pc[1]], y = pcs[,pc[2]])) +
    geom_point(aes(fill = `Region:Site`), shape = 21, size = 2, alpha = 0.75) +
    scale_fill_manual(values = site_colors) +
    coord_fixed(pcvar[pc[2]]/pcvar[pc[1]]) + 
    xlab(paste0( "PC", pc[1], " (", pcvar[pc[1]], "%)")) + 
    ylab(paste0( "PC", pc[2], " (", pcvar[pc[2]], "%)")) 
  #
  pca_plot
}

# Basic phylogenetic tree visualization
###------------------------------------------------------------------------------------
gg_tree <- function(treefile) {
  ggtree(treefile) + 
    geom_tiplab() +
    geom_treescale(width = 0.01, x = 0.05) +
    geom_rootedge(rootedge = 0.01)
}
