sites <- readr::read_tsv("~/computing/projects/etp_pocillopora_gskim/outputs/angsd/final_noclones_pgra_ibs05.sites.txt", col_names = c("chr", "site"))

N_sites <- sites %>% group_by(chr) %>% count()

N_sites_classes <- N_sites %>% 
  mutate("class" = ifelse(n < 10000, "scaffold", "chromosome"))

N_sites_classes %>% group_by(class) %>% count()
# 14 chromosome-scale contigs, 128 unplaced scaffolds

class_sums <- N_sites_classes %>% group_by(class) %>% summarise(sumn = sum(n)) %>% arrange(class)

sum(class_sums$sumn)

classes <- sort(unique(N_sites_classes$class))

# write out files with chr, pos info for all unplaced scaffolds

for (i in seq(2:(length(classes)))){
  print(classes[i])
  class_sites <- sites[sites$chr %in% N_sites_classes[N_sites_classes$class == classes[i], ]$chr, ]
  print(nrow(class_sites))
  class <- classes[i]
  filename <- paste0("~/Desktop/LD_pruning/class_", i, "_", class, "_sites.txt")
  print(filename)
  write_tsv(class_sites, filename, col_names = F)
}

# write out files with chr, pos info for scaffolds with >10K SNPs
tail(classes, 1)
chrs <- N_sites_classes[N_sites_classes$class == tail(classes, 1), ]$chr

for (i in chrs){
  print(i)
  chr_sites <- sites[sites$chr == i, ]
  filename <- paste0("~/Desktop/LD_pruning/", i, "_sites.txt")
  print(filename)
  write_tsv(chr_sites, filename, col_names = F)
}
