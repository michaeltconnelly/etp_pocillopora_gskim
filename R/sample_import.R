# Sample import script

metadata_col_types <- cols(
  `Sample ID` = col_character(),
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
  PANAMA_MANUSCRIPT = col_logical(),
  REPLICATE = col_logical(),
  `Region:Site` = col_character(),
  haplotype = col_character(),
  mtorf_type = col_character(),
  Species = col_character()
)

samples_r1 <- readr::read_csv("./data/pocillopora_samples_r1_metadata.csv", col_types = metadata_col_types)
samples_r2 <- readr::read_csv("./data/pocillopora_samples_r2_metadata.csv", col_types = metadata_col_types)
samples_r3 <- readr::read_csv("./data/pocillopora_samples_r3_metadata.csv", col_types = metadata_col_types)
all_samples <- dplyr::bind_rows(samples_r1, samples_r2, samples_r3, .id = "Sequence Round") %>% relocate(`Sequence Round`, .after = `Sample ID`) %>% arrange(`Sample ID`) 
