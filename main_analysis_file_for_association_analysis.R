
setwd("gotodirectory")

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(splitstackshape)
#library(PheWAS)
library(broom)
library(forestplot)

# ---- Read genotype data ----
#using plink to create dosage file and make sure you also have pt_id column
load("gotodirectory/genotype_data.raw.RData")

print(colnames(dat))
# [1] "FID"                  "IID"                  "PAT"                  "MAT"                  "SEX"                 
# [6] "PHENOTYPE"            "chr1.98096704.T.C_T"  "chr1.177916103.C.T_C" "chr1.188520686.A.C_A" "pt_id" 

# ---- Read covariate data ----
load("gotodirectory/covariate_data.RData")

head(covariate_file)
#       pt_id    sex  age  bmi
# 1 PT0000001   Male 56.6 21.3
# 2 PT0000002 Female 63.0 10.7
# 3 PT0000003 Female 48.8 24.6
# 4 PT0000004 Female 46.6 27.7
# 5 PT0000005 Female 41.8 29.1
# 6 PT0000006   Male 48.9 25.0

# ---- Read phenotype data ----
raw_file <- "gotodirectory/soma1100.csv"

# ---- Identify interested biomarker columns from raw_file ----
colnames_raw <- names(fread(raw_file, nrows = 0))
phenotype_columns <- colnames_raw[4:length(colnames_raw)]  #4 is location of the column for the proteomics data to begin with
print(phenotype_columns)
# [1] "Biomarker_1"    "Biomarker_2"    "Biomarker_3"    "Biomarker_4"    "Biomarker_5"    "Biomarker_6"    "Biomarker_7"   
# [8] "Biomarker_8"    "Biomarker_9"    "Biomarker_10"   "Biomarker_11"   "Biomarker_12"   "Biomarker_13"   "Biomarker_14"  

phenotype_data <- fread(raw_file, select = c("pt_id", "Biomarker_3"))  #using Biomarker_3 as example 
names(phenotype_data)[2] <- "Phenotype"  # standardize column name for model

# Replace values ≤ 0 with NA
phenotype_data$Phenotype[phenotype_data$Phenotype <= 0] <- NA

# # Merge with covariate data
phenotype_data <- merge(covariate_file[, c("pt_id", "sex", "age", "bmi")], phenotype_data, by = "pt_id")
female_only <- phenotype_data[phenotype_data$sex == "Female", ]$pt_id

# Drop rows with missing biomarker data
phenotype_data <- phenotype_data[!is.na(phenotype_data$Phenotype), ]

#merge with genotype data
dat_select <- dat[dat$pt_id %in% phenotype_data$pt_id, ]
#optional step to limit sample to female only
dat_select_female <- dat_select[dat_select$pt_id %in% female_only, ]

head(dat_select)
# FID       IID PAT MAT SEX PHENOTYPE chr1.98096704.T.C_T chr1.177916103.C.T_C chr1.188520686.A.C_A     pt_id
# 1 561 ID0000001   1   0   2         2                   2                    1                    0 PT0000001
# 2 997 ID0000002   0   0   1         1                   0                    0                    2 PT0000002
# 3 321 ID0000003   0   1   2         2                   2                    0                    2 PT0000003
# 4 153 ID0000004   1   1   1         1                   0                    0                    0 PT0000004
# 5  74 ID0000005   1   1   2         2                   1                    0                    0 PT0000005
# 6 228 ID0000006   1   0   2         1                   0                    0                    2 PT0000006

head(dat_select_female)
# FID       IID PAT MAT SEX PHENOTYPE chr1.98096704.T.C_T chr1.177916103.C.T_C chr1.188520686.A.C_A     pt_id
# 2 997 ID0000002   0   0   1         1                   0                    0                    2 PT0000002
# 3 321 ID0000003   0   1   2         2                   2                    0                    2 PT0000003
# 4 153 ID0000004   1   1   1         1                   0                    0                    0 PT0000004
# 5  74 ID0000005   1   1   2         2                   1                    0                    0 PT0000005
# 7 146 ID0000007   1   0   2         1                   1                    0                    2 PT0000007
# 8 634 ID0000008   1   1   2         1                   2                    1                    1 PT0000008

#Before you are going to the following step, you need to see the dat_select and dat_select_female something like above. This is sanity check.

###########################################################################################################
#create a list of names for the candidate variant with allele dosage for the allele of interested
variants <- c("chr1.98096704.T.C_T", "chr1.177916103.C.T_C", "chr1.188520686.A.C_A") #this name have to match the genotyping dosage file you created

result_list <- list()
row_index <- 1

for (i in seq_along(phenotype_columns)) {
  phenotype_name <- phenotype_columns[i]
  
  # Load phenotype data
  phenotype_data <- fread(raw_file, select = c("pt_id", phenotype_name))  
  
  # Clean and transform phenotype values
  phenotype_data[[phenotype_name]] <- as.numeric(phenotype_data[[phenotype_name]])
  phenotype_data[[phenotype_name]][phenotype_data[[phenotype_name]] <= 0] <- NA
  
  # log transformation
  # phenotype_data[[phenotype_name]] <- log(phenotype_data[[phenotype_name]])
  
  # INT transformation
  phenotype_data[[phenotype_name]] <- rank_inverse_normal(phenotype_data[[phenotype_name]])
  
  # Merge with covariates
  phenotype_data_select <- merge(
    covariate_file[, c("pt_id", "sex", "age", "bmi")],
    phenotype_data,
    by = "pt_id"
  )
  
  # Select female-only individuals
  female_only <- phenotype_data_select[phenotype_data_select$sex == "Female", ]$pt_id
  
  #Merge with genotype data
  dat_select <- dat[dat$pt_id %in% phenotype_data_select$pt_id, ]
  # dat_select_female <- dat_select[dat_select$pt_id %in% female_only, ]
  
  # Prepare unrelated individuals if needed
  # dat_select <- dat_select[!(dat_select$pt_id %in% c(c("PT1234567", "12345678", "123456789"))), ] #remove a list of related individuals
  # dat_select_female <- dat_select_female[!(dat_select_female$pt_id %in% c("PT1234567", "12345678", "123456789")), ]
  
  for (j in seq_along(variants)) {
    variant_name <- variants[j]
    
    tryCatch({
      message(paste0("Running: phenotype = ", phenotype_name, " | variant = ", variant_name))
      
      # Merge phenotype + covariates + genotype
      # phenotype_data_select_dat <- merge(phenotype_data_select, dat_select_female, by = "pt_id") # female only
      phenotype_data_select_dat <- merge(phenotype_data_select, dat_select, by = "pt_id") #both sex
      
      if (!(variant_name %in% colnames(phenotype_data_select_dat))) {
        stop(paste("Column", variant_name, "not found in merged dataset"))
      }
      
      # Build formula with covariates
      # formula <- as.formula(paste(phenotype_name, "~", variant_name, "+ age + bmi")) # for female only
      formula <- as.formula(paste(phenotype_name, "~", variant_name, "+ age + sex + bmi"))
      
      # Fit linear regression model
      model <- lm(formula = formula, data = phenotype_data_select_dat)
      summary_model <- summary(model)
      
      # Extract results for variant
      pvalue <- coef(summary_model)[variant_name, "Pr(>|t|)"]
      coeff <- coef(summary_model)[variant_name, "Estimate"]
      conf_int <- suppressMessages(confint(model))[variant_name, ]
      coeff_lower <- conf_int[1]
      coeff_upper <- conf_int[2]
      
      result_list[[row_index]] <- c(phenotype_name, variant_name, pvalue, coeff, coeff_lower, coeff_upper)
      
    }, error = function(e) {
      message(paste("❌ Error in:", phenotype_name, "vs", variant_name, "—", e$message))
      result_list[[row_index]] <- c(phenotype_name, variant_name, NA, NA, NA, NA)
    })
    
    row_index <- row_index + 1
    gc()
  }
}

# Combine results into a data frame
result_df <- as.data.frame(do.call(rbind, result_list), stringsAsFactors = FALSE)
colnames(result_df) <- c("phenotype", "genotype", "pvalue", "coeff", "coeff_lower", "coeff_upper")

head(result_df)
# phenotype             genotype             pvalue                coeff          coeff_lower          coeff_upper
# 1 Biomarker_1  chr1.98096704.T.C_T 0.0827789127837476   0.0298494054783445 -0.00387592644161281   0.0635747373983017
# 2 Biomarker_1 chr1.177916103.C.T_C  0.774125225058055 -0.00492177421824657  -0.0385415727964441   0.0286980243599509
# 3 Biomarker_1 chr1.188520686.A.C_A 0.0413096622492923   -0.035566293243477  -0.0697298114635412 -0.00140277502341288
# 4 Biomarker_2  chr1.98096704.T.C_T  0.279059971476284   0.0186184037449892  -0.0150984755924819   0.0523352830824604
# 5 Biomarker_2 chr1.177916103.C.T_C  0.817507885490272  0.00395568443647339  -0.0296496025652501   0.0375609714381969
# 6 Biomarker_2 chr1.188520686.A.C_A  0.872232048909384  0.00280267948242533  -0.0313601401117873    0.036965499076638

#create a forest plot for the candidate variants associated with biomarker named Biomarker_3 as an example
result_df_biomarker <- result_df[result_df$phenotype == "Biomarker_3", ]
result_df_biomarker$logp <- -log10(as.numeric(result_df_biomarker$pvalue))
forestplot <- ggplot(result_df_biomarker, aes(x = as.numeric(coeff), y = reorder(genotype, coeff))) +
  geom_point(aes(size = logp), color = "lightblue") +
  geom_errorbarh(aes(xmin = as.numeric(coeff_lower), xmax = as.numeric(coeff_upper)), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_size_continuous(name = expression(-log[10](pvalue))) +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  labs(
    title = "Coefficients of top hits for Biomarker_3",
    x = "Coefficient",
    y = "Variants"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0.5)
  ) +
  guides(color = guide_legend(title = NULL))
print(forestplot)

