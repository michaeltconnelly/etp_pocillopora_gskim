#pocillopora_etp_gksim_functions.R
#author: "Mike Connelly"
#date: "05/11/2023"

# project-wide theme, colors and factors
###----------------------------------------------------------------------------------------------------------------
theme_set(theme_bw())

region_order_longitude <- c("Clipperton", "Galapagos Northern", "Galapagos Central", "Coco", "Bahia Culebra", "Golfo Dulce", "Chiriqui", "Panama", "Gorgona")

region_colors<- c("Bahia Culebra"="aquamarine3",
                  "Chiriqui"="coral",
                  "Clipperton"="goldenrod1",
                  "Coco"="dodgerblue",
                  "Galapagos Northern"="lightslateblue",
                  "Galapagos Central"="cyan3",
                  "Golfo Dulce"="lightcoral",
                  "Gorgona"="mediumaquamarine",
                  "Panama"="cornflowerblue")

region_colors<- c("Bahia Culebra"="darkgoldenrod3",
                  "Chiriqui"="coral2",
                  "Clipperton"="deepskyblue",
                  "Coco"="lightseagreen",
                  "Galapagos Northern"="dodgerblue",
                  "Galapagos Central"="dodgerblue3",
                  "Golfo Dulce"="goldenrod1",
                  "Gorgona"="orchid",
                  "Panama"="lightcoral")

# 
site_colors <- c("gold1", "blue", "darkblue", "steelblue", "lightblue", "darkorange3","darkorange", "orange", "orangered", "cyan", "cyan3", "cyan4", "turquoise", "lightgreen")
names(site_colors)
#
mtorf_colors <- c("darkorange", "purple3", "turquoise2", "turquoise4", "darkslategrey")
#
spp_colors <- c("purple3", "darkorange", "gold2", "turquoise2", "turquoise4")
#
ngsadmix_clusters <- c("Cluster 1",
                       "Cluster 2",
                       "Cluster 3A",
                       "Cluster 3B", 
                       "Cluster 4")
#
ngsadmix_spp <- c("P. effusa",
                  "P. meandrina",
                  "P. grandis - Offshore",
                  "P. grandis - Continent",
                  "P. verrucosa 3a",
                  "Unassigned")
#
ngsadmix_spp_clusters <- c("Cluster 1" = "P. effusa",
                       "Cluster 2" = "P. meandrina",
                       "Cluster 3A" = "P. grandis - Offshore",
                       "Cluster 3B" = "P. grandis - Continent", 
                       "Cluster 4" = "P. verrucosa 3a")
#
ngsadmix_pop_colors <- c("P. effusa"="purple3",
                         "P. meandrina"="gold2",
                         "P. grandis - Offshore"="darkorange",
                         "P. grandis - Continent"="darkorange3",
                         "P. verrucosa 3a"="turquoise2",
                         "Unassigned"="grey90")

ngsadmix_pop_colors_k6 <- c("P. effusa"="purple3",
                         "P. meandrina"="gold2",
                         "P. grandis - Offshore"="darkorange",
                         "P. grandis - Continent"="darkorange3",
                         "P. verrucosa 3a"="turquoise2",
                         "P. verrucosa 3b"="turquoise4",
                         "Unassigned"="grey90")

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
