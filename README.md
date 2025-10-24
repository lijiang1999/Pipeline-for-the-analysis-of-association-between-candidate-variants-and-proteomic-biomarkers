# Pipeline-for-the-analysis-of-association-between-candidate-variants-and-proteomic-biomarkers
This is a part of pipeline codes to determine the association between candidate variants and proteomic biomarkers using a linear regression model after adjustment for covariates

1. run 'generation_of_synthetic_data.R' to get genotype, phenotype, and covariate files
2. run 'main_analysis_file_association_analysis.R' to do the analysis and create forest plot to show the association.
3. There is an example of output image file named 'Rplot_forestplot_example1.png'

During the process of the association analysis, you would see the progress like the following: 
Running: phenotype = Biomarker_1 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_1 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_1 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_2 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_2 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_2 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_3 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_3 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_3 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_4 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_4 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_4 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_5 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_5 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_5 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_6 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_6 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_6 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_7 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_7 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_7 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_8 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_8 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_8 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_9 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_9 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_9 | variant = chr1.188520686.A.C_A
Running: phenotype = Biomarker_10 | variant = chr1.98096704.T.C_T
Running: phenotype = Biomarker_10 | variant = chr1.177916103.C.T_C
Running: phenotype = Biomarker_10 | variant = chr1.188520686.A.C_A
