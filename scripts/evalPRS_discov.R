library(ggplot2)

sessionInfo()

options(warn = -1)

option_list = list(
  make_option(c('-p', '--phenotype'), action='store', type='character', help='.tsv file with phenotype (first 2 columns could be FID and IID, 3rd column should be the phenotype. If more than 1 phenotype is present, specify the column to be considered using pheno_col)', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/validated_input/processed.pheno.validated.txt'),
  make_option(c('--pheno_col'), action='store', type='character', help='If phenotype is provided, but contains more than 1 phenotype, the phenotype that is to be used can be specified here.', default='S1701_Z49'),
  make_option(c('-c', '--covariates'), action='store', type='character', help='.tsv file with covariates (first 2 columns could be FID and IID, 3rd column onwards will be considered as covariates.)', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/validated_input/cov_newset_chp_F20208_Long_axis_heart_images_DICOM_H5v3_NOnoise.cov.validated.txt'),
  make_option(c('--cov_cols'), action='store', type='character', help="Coma-seperated list of covaraite columns. If left blank, all columns will be used.", default='Sex,MRI_Visit,Age,MRI_Date,MRI_Centre,BSA'),
  make_option(c('--cov_cat'), action='store', type='character', help="Coma-seperated list of categorical covaraite columns. If left blank, all columns will be considered continous.", default='Sex,MRI_Visit,MRI_Date,MRI_Centre'),
  make_option(c('--cov_logical'), action='store', type='character', help="Coma-seperated list of logical covaraite columns. If left blank, only cov_cols and cov_cat will be considered."),
  make_option(c('--cov_PC'), action='store', type='character', help=".rds file containing the principal components of the genotype. If supplied, this will be used during the training of the sparse linear regression model.", default = "/project/ukbblatent/clinicaldata/v1.1.0_seventh_basket/genPC_82779_MD_01_03_2024_00_05_30.tsv"),
  make_option(c('--cov_nPC'), action='store', type='numeric', help='Number of PCs to include from the supplied cov_PC (Default: 20, same as the PLR code]', default=20),
  make_option(c('--subs2include'), action='store', type='character', help='txt file (or, output of plink) with FID and IID columns (tab-separated), containing subjects to include (typically used for relatedness filtering)', default='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/king_cutoff_0p0625_cond_plus_plink_maf1p_geno10p_caucasian_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.king.cutoff.in.id'),
  
  
  
  make_option(c('--disease_csv'), action='store', type='character', help='csv file (or, output of plink) with IID column (coma-separated), and BinCAT_Disease containing 1 or 0 for the particular disease', default='/project/ukbblatent/clinicaldata/binary_disease_cohorts/F20208v3_DiffAE/discov/atrial-fibrillation_flutter.csv'),
  
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]', default="/group/glastonbury/soumick/PRS/CnT/initial_external_FreqFilt_extGWAS_S1701_Z49")
)

opt = parse_args(OptionParser(option_list=option_list))


now<-Sys.time()
message('[',now,'][Message] performing pre-flight checks') 
check_output(opt)
opt$output<-paste0(opt$output,'.seed', opt$seed)
message('[',now,'][Message] done') 

#Read the phenotype and covariates
pheno<-data.frame(fread(opt$phenotype))

pheno_col_index <- which(names(pheno) == opt$pheno_col)
pheno <- pheno[, c(1:2, pheno_col_index)]

cov<-data.frame(fread(opt$covariates))
if (!is.null(opt$cov_cols)){
  message('[',now,'][Message] cov_cols supplied, selecting only the provided columns from the covariates table')
  cov <- cov[c(names(cov)[1:2], strsplit(opt$cov_cols, ",")[[1]])]
}
cov_cols <- colnames(cov)[-c(1, 2)]

merged_data <- merge(pheno, cov, by = names(pheno)[1:2])
rownames(merged_data) <- merged_data$IID

#prepare categorical covariates
if (!is.null(opt$cov_cat)){
  message('[',now,'][Message] cov_cat supplied, converting the provided columns as factors')
  categorical_covars <- strsplit(opt$cov_cat, ",")[[1]]
  merged_data[categorical_covars] <- lapply(merged_data[categorical_covars], factor)
}

#prepare logical covariates
if (!is.null(opt$cov_logical)){
  message('[',now,'][Message] cov_logical supplied, converting the provided columns as logical variables')
  logical_covars <- strsplit(opt$cov_logical, ",")[[1]]
  merged_data[logical_covars] <- lapply(merged_data[logical_covars], as.logical)
}

#drop covs where there is only one unique value
to_drop <- sapply(merged_data[cov_cols], function(x) length(unique(x)) < 2)
merged_data <- merged_data[, !(names(merged_data) %in% names(to_drop)[to_drop])]
cov_cols <- cov_cols[!to_drop[cov_cols]]

now<-Sys.time()
message('[',now,'][Message] phenotype and covariates tables read and processed')

if (!is.null(opt$subs2include)){
  message('[',now,'][Message] subs2include supplied, subsetting the pheno+cov table to only include the provided subject IDs')
  filt_subIDs <- read.delim(opt$subs2include, header = TRUE, sep = "\t")
  merged_data <- merged_data[merged_data$FID %in% filt_subIDs[, 1] & merged_data$IID %in% filt_subIDs[, 2], ]
}

# #if cov_PC is supplied, read and process
# if (!is.null(opt$cov_PC)){
#   
#   message('[',now,'][Message] cov_PC supplied, adding the principal components as a covariate (will not be used for baseline model)')
#   PC <- load_covPC(opt$cov_PC, opt$cov_nPC)
#   cov_PC_cols <- colnames(PC)
#   PC$IID <- rownames(PC)
#   merged_data <- merge(merged_data, PC, by="IID")
#   
# }


disease <- read.table(opt$disease_csv, sep = ",", header=1)
disease <- disease[,c("IID","BinCAT_Disease")]

###LDPred2 & PLR
data.full <- readRDS("/group/glastonbury/soumick/PRS/LDPred2/initial_external_basic_FiltAF_S1701_Z49.seed1701.fullDS.resNdata.basic.LDPred2.rds")
# data.full <- readRDS("/group/glastonbury/soumick/PRS/filteredPLR/initial_external_exGWAS_MfiltMAF_S1701_Z49.seed.pred.external.rds")
# data.full <- readRDS("/group/glastonbury/soumick/PRS/CnT/initial_external_FreqFilt_extGWAS_S1701_Z49.seed1701.prs.fullDS.rds")
# data.full <- readRDS("/group/glastonbury/soumick/PRS/simple/initial_external_FreqFilt_simple_S1701_Z49.seed1701.resNdata.simplePRS.rds")
prs <- data.full[,c("IID","prs_select")]
colnames(prs) <- c("IID", "PRS")

prs_disease <- merge(prs, disease, by="IID")
###############

#box plot
ggplot(prs_disease, aes(x = as.factor(BinCAT_Disease), y = PRS, fill = as.factor(BinCAT_Disease))) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("blue", "red"), name = "BinCAT_Disease") +
  labs(title = "Comparison of PRS Distributions",
       x = "BinCAT_Disease",
       y = "PRS") +
  theme_minimal()

#violin plot
ggplot(prs_disease, aes(x = as.factor(BinCAT_Disease), y = PRS, fill = as.factor(BinCAT_Disease))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +  # Adding a narrow box plot inside for median/IQR
  scale_fill_manual(values = c("blue", "red"), name = "BinCAT_Disease") +
  labs(title = "Comparison of PRS Distributions",
       x = "BinCAT_Disease",
       y = "PRS") +
  theme_minimal()

#density plot
ggplot(prs_disease, aes(x = PRS, fill = as.factor(BinCAT_Disease))) +
  geom_density(alpha = 0.5) + # Adjust alpha for transparency
  scale_fill_manual(values = c("blue", "red"), name = "BinCAT_Disease") +
  labs(title = "Density Plot of PRS by BinCAT_Disease Status",
       x = "PRS",
       y = "Density") +
  theme_minimal() +
  guides(fill=guide_legend(title="BinCAT_Disease"))

phenocov <- data.full[,c("IID",opt$pheno_col,cov_cols)]

prs_disease_phenocov <- merge(prs_disease, phenocov, by = "IID")

#model_formula <- as.formula(paste("BinCAT_Disease ~", opt$pheno_col, "+", paste(cov_cols, collapse = "+")))
model_formula <- paste("BinCAT_Disease ~", opt$pheno_col)
model_without_PRS <- glm(model_formula, data = prs_disease_phenocov, family = binomial())

#model_formula <- as.formula(paste("BinCAT_Disease ~", opt$pheno_col, "+", paste(cov_cols, collapse = "+"), "+ PRS"))
model_formula <- paste("BinCAT_Disease ~", opt$pheno_col, "+ PRS")
model_with_PRS <- glm(model_formula, data = prs_disease_phenocov, family = binomial())

#model_formula <- as.formula(paste("BinCAT_Disease ~", "PRS", "+", paste(cov_cols, collapse = "+")))
model_formula <- paste("BinCAT_Disease ~", "PRS")
model_PRS_only <- glm(model_formula, data = prs_disease_phenocov, family = binomial())

prediction_without_PRS <- predict(model_without_PRS, prs_disease_phenocov, type = "response")
roc_without_PRS <- roc(prs_disease_phenocov$BinCAT_Disease, prediction_without_PRS)
auc_without_PRS <- auc(roc_without_PRS)

prediction_with_PRS <- predict(model_with_PRS, prs_disease_phenocov, type = "response")
roc_with_PRS <- roc(prs_disease_phenocov$BinCAT_Disease, prediction_with_PRS)
auc_with_PRS <- auc(roc_with_PRS)

prediction_PRS_only <- predict(model_PRS_only, prs_disease_phenocov, type = "response")
roc_PRS_only <- roc(prs_disease_phenocov$BinCAT_Disease, prediction_PRS_only)
auc_PRS_only <- auc(roc_PRS_only)

# Compare AUC values
print(paste("AUC without PRS:", auc_without_PRS))
print(paste("AUC with PRS:", auc_with_PRS))
print(paste("AUC PRS only:", auc_PRS_only))

test_result <- roc.test(roc_without_PRS, roc_with_PRS)
