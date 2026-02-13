library(tidyverse)

# Load environment variables from .env file
if (file.exists(".env")) {
  env_vars <- readLines(".env")
  env_vars <- env_vars[!grepl("^#", env_vars) & nchar(env_vars) > 0]
  for (var in env_vars) {
    parts <- strsplit(var, "=", fixed = TRUE)[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      value <- trimws(gsub("^\"|\"$", "", parts[2]))
      do.call(Sys.setenv, setNames(list(value), key))
    }
  }
}

# Read paths from environment variables
adnimerge_path <- Sys.getenv("ADNIMERGE_DATA")
fam_file <- paste0(Sys.getenv("INPUT_PLINK_SUFFIX"), ".fam")
output_dir <- "temp"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Install and load adnimerge once
# install.packages(adnimerge_path, repos = NULL, type = "source")
library(ADNIMERGE)

# Load ADNI data
data(adnimerge)
data(gdscale)
data(npi)

cat("=== Data loaded ===\n")
cat("adnimerge dimensions:", dim(adnimerge), "\n")
cat("gdscale dimensions:", dim(gdscale), "\n")

# Load genetic data (fam file) to identify subjects with genetic data
fam_data <- read.table(fam_file, header = FALSE, 
                       col.names = c("FID", "IID", "father", "mother", "sex", "phenotype"))

cat("\nGenetic data (fam file) loaded\n")
cat("Number of individuals with genetic data:", nrow(fam_data), "\n")
cat("Sample IIDs:\n")
print(head(fam_data$IID))

# Extract PTID from IID (remove prefix before first '_')
fam_data <- fam_data %>%
  mutate(PTID = sub("^[^_]+_", "", IID)) %>%
  select(IID, PTID)  # Keep only IID and PTID for joining

cat("\nSample PTIDs after prefix removal:\n")
print(head(fam_data$PTID))

# Select relevant variables from adnimerge
adni_subset <- adnimerge %>%
  select(RID, COLPROT, ORIGPROT, PTID, VISCODE, DX, DX.bl, AGE, 
         PTGENDER, PTEDUCAT, SITE, APOE4) %>%
  mutate(
    RID = as.character(RID),
    VISCODE = as.character(VISCODE)
  ) %>% distinct(RID, VISCODE, .keep_all = TRUE)  # Keep only unique RID-VISCODE combinations

# Select relevant variables from gdscale (GDS5 = GDS >= 5, GDS10 = GDS >= 10)
gds_subset <- gdscale %>%
  select(RID, VISCODE, GDTOTAL) %>%
  filter(!is.na(GDTOTAL)) %>%  # Keep only rows with non-missing GDTOTAL
  mutate(
    RID = as.character(RID),
    VISCODE = as.character(VISCODE),
    GDS5 = if_else(GDTOTAL >= 5, 1, 0),
    GDS10 = if_else(GDTOTAL >= 10, 1, 0)
  ) %>% 
  # Convert 'sc' to 'bl' to join baseline visits
  mutate(VISCODE = if_else(VISCODE == "sc", "bl", VISCODE)) %>%  
  distinct(RID, VISCODE, .keep_all = TRUE)  # Keep only unique RID-VISCODE combinations

# Select NPI variables from  npi (NPI_DPR = NPI Depression, NPI_ANX = NPI Anxiety)
npi_subset <- npi %>%
  select(RID, VISCODE, NPID, NPIE) %>%
  filter(!is.na(NPID) & !is.na(NPIE)) %>%  # Keep only rows with non-missing NPID and NPIE
  mutate(
    RID = as.character(RID),
    VISCODE = as.character(VISCODE),
    NPI_DPR = case_when(NPID=='Yes' ~ 1, NPID=='No' ~ 0, TRUE ~ NA_real_),
    NPI_ANX = case_when(NPIE=='Yes' ~ 1, NPIE=='No' ~ 0, TRUE ~ NA_real_)
  ) %>% 
  distinct(RID, VISCODE, .keep_all = TRUE)  # Keep only unique RID-VISCODE combinations


cat("\n=== Joining datasets ===\n")

# Join adnimerge and gdscale
pheno <- adni_subset %>%
  left_join(gds_subset, by = c("RID", "VISCODE")) %>%
  left_join(npi_subset, by = c("RID", "VISCODE")) %>%
  mutate(PTID = as.character(PTID))

cat("After joining adnimerge, gdscale, and npi:", nrow(pheno), "rows\n")

# Keep only subjects with genetic data and add IID
pheno <- pheno %>%
  filter(PTID %in% fam_data$PTID) %>%
  left_join(fam_data, by = "PTID")

cat("After filtering for genetic data availability:", nrow(pheno), "rows\n")
cat("Unique subjects:", length(unique(pheno$PTID)), "\n")

# Create timepoint variable
pheno <- pheno %>%
  mutate(timepoint = case_when(
    VISCODE == "bl" ~ 0,
    grepl("^m", VISCODE) ~ as.numeric(gsub("m", "", VISCODE)),
    TRUE ~ as.numeric(VISCODE)
  )) %>% filter(!is.na(timepoint))  # Keep only rows with valid timepoints

cat("\n=== Timepoint distribution ===\n")
print(table(pheno$timepoint, useNA = "ifany"))

# Outcome distribution checks
cat("\n=== Depression variable distribution ===\n")
cat("\nGDS >= 5:\n")
print(table(pheno$GDS5, useNA = "ifany"))
cat("\nGDS >= 10:\n")
print(table(pheno$GDS10, useNA = "ifany"))
cat("\n=== NPI Depression ===\n")
print(table(pheno$NPI_DPR, useNA = "ifany"))
cat("\n=== NPI Anxiety ===\n")
print(table(pheno$NPI_ANX, useNA = "ifany"))

# Distribution over timepoints
cat("\n=== GDS >= 5 by timepoint ===\n")
print(table(pheno$timepoint, pheno$GDS5, useNA = "ifany"))

cat("\n=== GDS >= 10 by timepoint ===\n")
print(table(pheno$timepoint, pheno$GDS10, useNA = "ifany"))

cat("\nGDS availability by DX.bl:\n")
pheno %>% filter(!is.na(GDS5)) %>%
  with(table(timepoint, DX.bl)) %>% print()

cat("\n=== NPI Depression by timepoint ===\n")
print(table(pheno$timepoint, pheno$NPI_DPR, useNA = "ifany"))

cat("\n=== NPI Anxiety by timepoint ===\n")
print(table(pheno$timepoint, pheno$NPI_ANX, useNA = "ifany"))

cat("\nNPI_DPR availability by DX.bl:\n")
pheno %>% filter(!is.na(NPI_DPR)) %>%
  with(table(timepoint, DX.bl)) %>% print()

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("\nGDTOTAL summary:\n")
print(summary(pheno$GDTOTAL))

cat("\nMissing data summary:\n")
missing_summary <- pheno %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(pheno) * 100, 2)) %>%
  arrange(desc(n_missing))
print(missing_summary)

# Check protocol distribution
cat("\n=== Protocol distribution ===\n")
cat("ORIGPROT:\n")
print(table(pheno$ORIGPROT, useNA = "ifany"))
cat("\nCOLPROT:\n")
print(table(pheno$COLPROT, useNA = "ifany"))

# Diagnosis distribution
cat("\n=== Diagnosis distribution ===\n")
print(table(pheno$DX, useNA = "ifany"))

# Demographics
cat("\n=== Demographics ===\n")
cat("Age: Mean =", round(mean(pheno$AGE, na.rm = TRUE), 1), 
    "SD =", round(sd(pheno$AGE, na.rm = TRUE), 1), "\n")
cat("Gender:\n")
print(table(pheno$PTGENDER, useNA = "ifany"))
cat("Education: Mean =", round(mean(pheno$PTEDUCAT, na.rm = TRUE), 1),
    "SD =", round(sd(pheno$PTEDUCAT, na.rm = TRUE), 1), "\n")
cat("APOE4:\n")
print(table(pheno$APOE4, useNA = "ifany"))

cat("\n=== Data preview ===\n")
print(head(pheno, 10))

# Based on the distribution of timepoints and depression variables, 
# we will focus on the 12-month timepoint and 24-month cross-sectional data
pheno %>% filter(timepoint==12) %>%
    select(IID, GDS5, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
    filter(complete.cases(.)) %>% 
    distinct(IID, .keep_all = TRUE) %>%
    write.csv(file.path(output_dir, "gds5_y1.csv"), row.names = FALSE)

pheno %>% filter(timepoint==24) %>%
    select(IID, GDS5, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
    filter(complete.cases(.)) %>% 
    distinct(IID, .keep_all = TRUE) %>%
    write.csv(file.path(output_dir, "gds5_y2.csv"), row.names = FALSE)
# for DPI_DPR and NPI_ANX, we will use timepoint 72
pheno %>% filter(timepoint==72) %>%
    select(IID, NPI_DPR, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
    filter(complete.cases(.)) %>% 
    distinct(IID, .keep_all = TRUE) %>%
    write.csv(file.path(output_dir, "npi_dpr_y6.csv"), row.names = FALSE)

pheno %>% filter(timepoint==72) %>%
    select(IID, NPI_ANX, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
    filter(complete.cases(.)) %>% 
    distinct(IID, .keep_all = TRUE) %>%
    write.csv(file.path(output_dir, "npi_anx_y6.csv"), row.names = FALSE)


# NPI_DPR_EVER: Create a variable indicating if the subject ever had NPI_DPR=1
# keep the timepoint of first NPI_DPR or last follow-up if never
npi_depr_ever <- pheno %>%
  filter(!is.na(NPI_DPR)) %>%
  group_by(IID) %>%
  arrange(timepoint) %>%
  slice(if(any(NPI_DPR == 1)) which(NPI_DPR == 1)[1] else n()) %>%
  ungroup() %>%
  rename(NPI_DPR_EVER = NPI_DPR, timepoint_censored = timepoint) %>%
  select(IID, timepoint_censored, NPI_DPR_EVER, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
  filter(complete.cases(.))

cat("\nNPI_DPR_EVER distribution:\n")
print(table(npi_depr_ever$NPI_DPR_EVER, useNA = "ifany"))
npi_depr_ever %>%
  write.csv(file.path(output_dir, "npi_dpr_ever.csv"), row.names = FALSE)

# NPI_ANX_EVER: Create a variable indicating if the subject ever had NPI_ANX=1
# keep the timepoint of first NPI_ANX or last follow-up if never
npi_anx_ever <- pheno %>%
  filter(!is.na(NPI_ANX)) %>%
  group_by(IID) %>%
  arrange(timepoint) %>%
  slice(if(any(NPI_ANX == 1)) which(NPI_ANX == 1)[1] else n()) %>%
  ungroup() %>%
  rename(NPI_ANX_EVER = NPI_ANX, timepoint_censored = timepoint) %>%
  select(IID, timepoint_censored, NPI_ANX_EVER, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
  filter(complete.cases(.))

cat("\nNPI_ANX_EVER distribution:\n")
print(table(npi_anx_ever$NPI_ANX_EVER, useNA = "ifany"))
npi_anx_ever %>%
  write.csv(file.path(output_dir, "npi_anx_ever.csv"), row.names = FALSE)


cat('Output files saved in:', output_dir, "\n")
cat("gds5_y1.csv: Year 1 (12-month) GDS >= 5 phenotype\n")
cat("gds5_y2.csv: Year 2 (24-month) GDS >= 5 phenotype\n")
cat("npi_dpr_y6.csv: Year 6 (72-month) NPI Depression phenotype\n")
cat("npi_anx_y6.csv: Year 6 (72-month) NPI Anxiety phenotype\n")
cat("npi_dpr_ever.csv: NPI Depression EVER phenotype with timepoint of censoring (first occurrence or last follow-up)\n")
cat("npi_anx_ever.csv: NPI Anxiety EVER phenotype with timepoint of censoring (first occurrence or last follow-up)\n")  