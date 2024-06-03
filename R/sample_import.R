# Sample import script

metadata_col_types <- cols(
  `Sample ID` = col_character(),
  `Genotype ID` = col_character(),
  `Cattle Tag Number` = col_character(),
  `Biorepository Number` = col_character(),
  Country = col_character(),
  Region = col_character(),
  Site = col_character(),
  `Sampling Date` = col_character(),
  Latitude = col_double(),
  Longitude = col_double(),
  Set = col_character(),
  STRI = col_logical(),
  RSMAS = col_logical(),
  TECHNICAL_REPLICATE = col_logical(),
  SINGLE_GENOTYPE = col_logical(),
  CLONAL_GENOTYPE = col_logical(),
  NOCLONES_SAMPLE = col_logical(),
  `Region:Site` = col_character(),
  haplotype = col_character(),
  mtorf_type = col_character(),
  Species = col_character()
)

# Read in all sample metadata for the 3 sequencing rounds
samples_r1 <- readr::read_csv("./data/pocillopora_samples_r1_metadata.csv", col_types = metadata_col_types) # Panama, Isla del Coco, Clipperton
samples_r2 <- readr::read_csv("./data/pocillopora_samples_r2_metadata.csv", col_types = metadata_col_types) # Colombia, Galapagos, technical replicates
samples_r3 <- readr::read_csv("./data/pocillopora_samples_r3_metadata.csv", col_types = metadata_col_types) # Golfo Dulce, Bahia Culebra

# Bind all samples together into master data frame
all_samples <- dplyr::bind_rows(samples_r1, samples_r2, samples_r3, .id = "Sequence Round") %>%
  relocate(`Sequence Round`, .after = `Sample ID`) %>%
  arrange(`Sample ID`) 

# Filter to analysis samples
samples_analysis <- all_samples %>% filter(ANALYSIS) %>%
  arrange(`Sample ID`) 

# Filter to "no-clones" samples only
samples_noclones <- all_samples %>% filter(NOCLONES_SAMPLE) %>%
  arrange(`Sample ID`) 

# Join NGSAdmix species assignments for K=5 and K=6 to data frame
samples_k5_assignments <- read_csv("./data/samples_noclones_ngsadmix_assignments_k5.csv")
samples_k6_assignments <- read_csv("./data/samples_noclones_ngsadmix_assignments_k6.csv")
samples_noclones_ngsadmix <- left_join(samples_noclones, samples_k5_assignments, by = "Sample ID") %>% left_join(., samples_k6_assignments, by = "Sample ID") %>%
  arrange(`Sample ID`) 

