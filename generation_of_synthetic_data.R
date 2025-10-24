#simulation section for synthetic data to run the code

# Simulate genotype data
set.seed(42)
n <- 5000
pt_id <- sprintf("PT%07d", 1:n)

dat <- data.frame(
  FID = sample(1:1000, n, replace = TRUE),
  IID = sprintf("ID%07d", 1:n),
  PAT = sample(0:1, n, replace = TRUE),
  MAT = sample(0:1, n, replace = TRUE),
  SEX = sample(1:2, n, replace = TRUE),
  PHENOTYPE = sample(1:2, n, replace = TRUE),
  chr1.98096704.T.C_T = sample(0:2, n, replace = TRUE),
  chr1.177916103.C.T_C = sample(0:2, n, replace = TRUE),
  chr1.188520686.A.C_A = sample(0:2, n, replace = TRUE),
  pt_id = pt_id
)

save(dat, file = "gotodirectory/simulated_data/genotype_data.raw.RData")

# Simulate covariate data
covariate_file <- data.frame(
  pt_id = pt_id,
  sex = sample(c("Female", "Male"), n, replace = TRUE),
  age = round(rnorm(n, mean = 50, sd = 10), 1),
  bmi = round(rnorm(n, mean = 25, sd = 4), 1)
)

save(covariate_file, file = "gotodirectory/simulated_data/covariate_data.RData")

#Simulate phenotype data
# set.seed(42)
# n <- 5000
# pt_id <- sprintf("PT%07d", 1:n)

# Create metadata columns
phenotype_data <- data.frame(
  pt_id = pt_id,
  meta1 = sample(c("A", "B", "C"), n, replace = TRUE),
  meta2 = sample(100:200, n, replace = TRUE)
)

# Add 1100 biomarker columns with skewed values between ~100 and ~10,000
for (i in 1:1100) {
  biomarker_name <- paste0("Biomarker_", i)
  raw_values <- rlnorm(n, meanlog = 7, sdlog = 0.5)  # log-normal distribution
  phenotype_data[[biomarker_name]] <- pmin(pmax(raw_values, 100), 10000)
}

# Save to CSV
write.csv(phenotype_data, file = "gotodirectory/simulated_data/soma1100.csv", row.names = FALSE)
