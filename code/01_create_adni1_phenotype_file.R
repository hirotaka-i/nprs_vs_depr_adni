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
  )

# Select relevant variables from gdscale
gds_subset <- gdscale %>%
  select(RID, VISCODE, GDTOTAL) %>%
  mutate(
    RID = as.character(RID),
    VISCODE = as.character(VISCODE)
  )

cat("\n=== Joining datasets ===\n")

# Join adnimerge and gdscale
pheno <- adni_subset %>%
  left_join(gds_subset, by = c("RID", "VISCODE")) %>%
  mutate(PTID = as.character(PTID))

cat("After joining adnimerge and gdscale:", nrow(pheno), "rows\n")

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
  ))

cat("\n=== Timepoint distribution ===\n")
print(table(pheno$timepoint, useNA = "ifany"))

# Create depression variables
pheno <- pheno %>%
  mutate(
    DEPR1 = case_when(
      is.na(GDTOTAL) ~ NA_real_,
      GDTOTAL >= 5 ~ 1,
      GDTOTAL < 5 ~ 0
    ),
    DEPR2 = case_when(
      is.na(GDTOTAL) ~ NA_real_,
      GDTOTAL >= 10 ~ 1,
      GDTOTAL < 10 ~ 0
    )
  )

cat("\n=== Depression variable distribution ===\n")
cat("\nDEPR1 (GDS >= 5):\n")
print(table(pheno$DEPR1, useNA = "ifany"))

cat("\nDEPR2 (GDS >= 10):\n")
print(table(pheno$DEPR2, useNA = "ifany"))

# Distribution over timepoints
cat("\n=== DEPR1 by timepoint ===\n")
print(table(pheno$timepoint, pheno$DEPR1, useNA = "ifany"))

cat("\n=== DEPR2 by timepoint ===\n")
print(table(pheno$timepoint, pheno$DEPR2, useNA = "ifany"))

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
print(head(pheno, 20))

# Based on the distribution of timepoints and depression variables, 
# we will focus on the 12-month timepoint and 24-month cross-sectional data
pheno %>% filter(timepoint==12) %>%
    select(IID, DEPR1,
    AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
    filter(complete.cases(.)) %>% 
    distinct(IID, .keep_all = TRUE) %>%
    write.csv(file.path(output_dir, "depr1_y1.csv"), row.names = FALSE)

pheno %>% filter(timepoint==24) %>%
    select(IID, DEPR1,
    AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl) %>%
    filter(complete.cases(.)) %>% 
    distinct(IID, .keep_all = TRUE) %>%
    write.csv(file.path(output_dir, "depr1_y2.csv"), row.names = FALSE)