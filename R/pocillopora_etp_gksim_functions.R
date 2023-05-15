#pocillopora_etp_gksim_functions.R
#author: "Mike Connelly"
#date: "05/11/2023"

gg_pcangsd <- function(e, pc = c(1,2)) {
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
