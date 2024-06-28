sites <- readr::read_tsv("~/Desktop/final_noclones_ibs05.sites.txt")

N_sites <- sites %>% group_by(chr) %>% count()

N_sites_classes <- N_sites %>% 
  mutate("class" = ifelse(n < 200, "A-<200",
                          ifelse(n >=200 & n < 400, "A-200-400",
                          ifelse(n >= 400 & n < 600, "A-400-600",
                          ifelse(n >= 600 & n < 800, "A-600-800",
                          ifelse(n >= 800 & n < 1000, "A-800-1000",
                          ifelse(n >= 1000 & n < 1200, "B-1000-1200",
                          ifelse(n >= 1200 & n < 1400, "B-1200-1400",
                          ifelse(n >= 1400 & n < 1600, "B-1400-1600",
                          ifelse(n >= 1600 & n < 1800, "B-1600-1800",
                          ifelse(n >= 1800 & n < 2000, "B-1800-2000",   
                          ifelse(n >= 2000 & n < 2200, "C-2000-2200",  
                          ifelse(n >= 2200 & n < 2400, "C-2200-2400",  
                          ifelse(n >= 2400 & n < 2600, "C-2400-2600",  
                          ifelse(n >= 2600 & n < 2800, "C-2600-2800",   
                          ifelse(n >= 2800 & n < 3000, "C-2800-3000",        
                          ifelse(n >= 3000 & n < 3200, "D-3000-3200",
                          ifelse(n >= 3200 & n < 3400, "D-3200-3400",
                          ifelse(n >= 3400 & n < 3600, "D-3400-4000",
                          ifelse(n >= 3600 & n < 3800, "D-3600-3800",
                          ifelse(n >= 3800 & n < 4000, "D-3800-4000",
                          ifelse(n >= 4000 & n < 4200, "E-4000-4200",
                          ifelse(n >= 4200 & n < 4400, "E-4200-4400",
                          ifelse(n >= 4400 & n < 4600, "E-4400-4600",       
                          ifelse(n >= 4600 & n < 4800, "E-4600-4800",
                          ifelse(n >= 4800 & n < 5000, "E-4800-5000",
                          ifelse(n >= 5000 & n < 5200, "F-5000-5200",
                          ifelse(n >= 5200 & n < 5400, "F-5200-5400",
                          ifelse(n >= 5400 & n < 5600, "F-5400-5600",       
                          ifelse(n >= 5600 & n < 5800, "F-5600-5800",
                          ifelse(n >= 5800 & n < 6000, "F-5800-6000",
                          ifelse(n >= 6000 & n < 6200, "G-6000-6200",
                          ifelse(n >= 6200 & n < 6400, "G-6200-6400",
                          ifelse(n >= 6400 & n < 6600, "G-6400-6600",       
                          ifelse(n >= 6600 & n < 6800, "G-6600-6800",
                          ifelse(n >= 6800 & n < 7000, "G-6800-7000",    
                          ifelse(n >= 7000 & n < 7200, "H-7000-7200",
                          ifelse(n >= 7200 & n < 7400, "H-7200-7400",
                          ifelse(n >= 7400 & n < 7600, "H-7400-7600",       
                          ifelse(n >= 7600 & n < 7800, "H-7600-7800",
                          ifelse(n >= 7800 & n < 8000, "H-7800-8000",    
                          ifelse(n >= 8000 & n < 8200, "I-8000-8200",
                          ifelse(n >= 8200 & n < 8400, "I-8200-8400",
                          ifelse(n >= 8400 & n < 8600, "I-8400-8600",       
                          ifelse(n >= 8600 & n < 8800, "I-8600-8800",
                          ifelse(n >= 8800 & n < 8000, "I-8800-9000",    
                          ifelse(n >= 9000 & n < 9200, "J-9000-9200",
                          ifelse(n >= 9200 & n < 9600, "J-9200-9600",
                          ifelse(n >= 9600 & n < 10000, "J-9600-10000",       
                          "K->10K")))))))))))))))))))))))))))))))))))))))))))))))))

N_sites_classes %>% group_by(class) %>% count()

class_sums <- N_sites_classes %>% group_by(class) %>% summarise(sumn = sum(n)) %>% arrange(class)

sum(class_sums$sumn)

classes <- sort(unique(N_sites_classes$class))

# write out files with chr, pos info for different classes

for (i in seq(1:(length(classes)-1))){
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
