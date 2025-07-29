source("scripts/utils.R")

sessionInfo()

options(warn = -1)

option_list = list(
  make_option(c('-i', '--input'), action='store', type='character', help='.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]', default='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/cond_plus_plink_maf1p_geno10p_caucasian_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.bgen'),
  make_option(c('-p', '--phenotype'), action='store', type='character', help='.tsv file with phenotype (first 2 columns could be FID and IID, 3rd column should be the phenotype. If more than 1 phenotype is present, specify the column to be considered using pheno_col)', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/validated_input/processed.pheno.validated.txt'),
  make_option(c('--pheno_col'), action='store', type='character', help='If phenotype is provided, but contains more than 1 phenotype, the phenotype that is to be used can be specified here.', default='S1701_Z49'),
  make_option(c('-c', '--covariates'), action='store', type='character', help='.tsv file with covariates (first 2 columns could be FID and IID, 3rd column onwards will be considered as covariates.)', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/validated_input/cov_newset_chp_F20208_Long_axis_heart_images_DICOM_H5v3_NOnoise.cov.validated.txt'),
  make_option(c('--cov_cols'), action='store', type='character', help="Coma-seperated list of covaraite columns. If left blank, all columns will be used.", default='Sex,MRI_Visit,Age,MRI_Date,MRI_Centre,BSA'),
  make_option(c('--cov_cat'), action='store', type='character', help="Coma-seperated list of categorical covaraite columns. If left blank, all columns will be considered continous.", default='Sex,MRI_Visit,MRI_Date,MRI_Centre'),
  make_option(c('--cov_logical'), action='store', type='character', help="Coma-seperated list of logical covaraite columns. If left blank, only cov_cols and cov_cat will be considered."),
  make_option(c('--cov_PC'), action='store', type='character', help=".rds file containing the principal components of the genotype. If supplied, this will be used during the training of the sparse linear regression model.", default = "/project/ukbblatent/clinicaldata/v1.1.0_seventh_basket/genPC_82779_MD_01_03_2024_00_05_30.tsv"),
  make_option(c('--cov_nPC'), action='store', type='numeric', help='Number of PCs to include from the supplied cov_PC (Default: 20, same as the PLR code]', default=20),
  make_option(c('--subs2include'), action='store', type='character', help='txt file (or, output of plink) with FID and IID columns (tab-separated), containing subjects to include (typically used for relatedness filtering)', default='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/king_cutoff_0p0625_cond_plus_plink_maf1p_geno10p_caucasian_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.king.cutoff.in.id'),
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]', default="/group/glastonbury/soumick/PRS/PLR/initial_test_multiAlp_S1701_Z49"),
  make_option(c('--threads'), action='store', type='numeric', help='computing threads [1]', default=5),
  make_option(c('--seed'), action='store', type='numeric', help='set seed (to be used for train-test split', default=1701),
  make_option(c('-t', '--train'), action='store', type='double', help='train percentage for training-testing - internal validation', default = 0.8)
)

opt = parse_args(OptionParser(option_list=option_list))

now<-Sys.time()
message('[',now,'][Message] performing pre-flight checks') 
is_bgen<-check_input(opt$input)
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

now<-Sys.time()
message('[',now,'][Message] phenotype and covariates tables read and processed')

if (!is.null(opt$subs2include)){
  message('[',now,'][Message] subs2include supplied, subsetting the pheno+cov table to only include the provided subject IDs')
  filt_subIDs <- read.delim(opt$subs2include, header = TRUE, sep = "\t")
  merged_data <- merged_data[merged_data$FID %in% filt_subIDs[, 1] & merged_data$IID %in% filt_subIDs[, 2], ]
}

#if cov_PC is supplied, read and process
if (!is.null(opt$cov_PC)){
  
  message('[',now,'][Message] cov_PC supplied, adding the principal components as a covariate (will not be used for baseline model)')
  PC <- load_covPC(opt$cov_PC, opt$cov_nPC)
  cov_PC_cols <- colnames(PC)
  PC$IID <- rownames(PC)
  merged_data <- merge(merged_data, PC, by="IID")
  
}

#Read and prepare genotype data
now<-Sys.time()
message('[',now,'][Message] loading data') 
message('[',now,'][Message] reading  .bed/.bgen')

if (!is_bgen) {
  
  obj.bigSNP<-load_bed(opt$input, threads=opt$threads)
  
} else {
  
  obj_bgen<-load_bgen(opt$input,threads=opt$threads, subset_subIDs=merged_data$IID) #subsetting the genotypes to keep only the subjects present in the phenotypes table
  obj.bigSNP<-obj_bgen[[1]]
  obj.sample<-obj_bgen[[2]]
  
}

G <- tryCatch({
  snp_fastImputeSimple(obj.bigSNP$genotypes)
}, error = function(e) {
  obj.bigSNP$genotypes
})


if (is_bgen) {
  
  obj.bigSNP$fam <- snp_fake(n = nrow(G), m = 1)$fam
  
  obj.bigSNP$fam$family.ID <- obj.sample$ID_1
  obj.bigSNP$fam$sample.ID <- obj.sample$ID_2
}

CHR <- as.integer(obj.bigSNP$map$chromosome) #this is somehow necessary for .bgen files, not for bed
POS <- obj.bigSNP$map$physical.pos

now<-Sys.time()
message('[',now,'][Message] genotypes loaded and processed')

#read the sample and match 
merged_data <- subset(merged_data, IID %in% obj.sample$ID_2) #discarding phenotypes and covariates for the subjects not present in the genotype data
data.full <- merged_data[match(obj.sample$ID_2, merged_data$IID), ] #ensure that both the genotype and phenotype tables are in the same order
ind.subIDs <- data.full$IID

#Create train-test split
set.seed(opt$seed)
ind.train <- sort(sample(length(ind.subIDs), round(length(ind.subIDs) * opt$train)))
ind.test <- setdiff(seq_along(ind.subIDs), ind.train)
saveRDS(list(ind.subIDs=ind.subIDs, ind.train=ind.train, ind.test=ind.test), file=file.path(paste0(opt$output, ".ind.rds")))

data.train <- subset(data.full, IID %in% ind.subIDs[ind.train])
data.test <- subset(data.full, IID %in% ind.subIDs[ind.test])
saveRDS(list(data.full=data.full,data.train=data.train,data.test=data.test), file=file.path(paste0(opt$output, ".data.rds")))

#baseline model with linear regression
baseline_formula <- as.formula(paste(opt$pheno_col, "~", paste(cov_cols, collapse = "+")))
sink(file.path(paste0(opt$output,'.summary.mod.base.txt')))
summary(mod.base <- lm(baseline_formula, data=data.train))
sink()
pred.base <- predict(mod.base, data.full)

#Train the Sparse linear regression model
if (is.null(opt$cov_PC)){
  mod <- big_spLinReg(X = G, 
                      y.train = data.train[[opt$pheno_col]], 
                      ind.train = ind.train,
                      base.train = pred.base[ind.train],
                      dfmax = Inf,
                      ncores = opt$threads)
  pred.mod.train <- predict(mod, G,  ind.row = ind.train)
  pred.mod.test <- predict(mod, G,  ind.row = ind.test)
} else {
  cov.PC.train <- covar_from_df(data.train[, cov_PC_cols])
  cov.PC.test <- covar_from_df(data.test[, cov_PC_cols])
  mod <- big_spLinReg(X = G, 
                      y.train = data.train[[opt$pheno_col]], 
                      ind.train = ind.train,
                      covar.train = cov.PC.train,
                      base.train = pred.base[ind.train],
                      alphas = seq(0.1, 1, by=0.1),
                      dfmax = Inf,
                      ncores = opt$threads)
  pred.mod.train <- predict(mod, G, ind.row = ind.train, covar.row = cov.PC.train, base.row = pred.base[ind.train])
  pred.mod.test <- predict(mod, G, ind.row = ind.test, covar.row = cov.PC.test, base.row = pred.base[ind.test])
}

saveRDS(list(mod=mod,mode.base=mod.base), file=file.path(paste0(opt$output, ".mod.rds")))

pdf(file.path(paste0(opt$output,'.mod.pdf')))
plot(mod)
dev.off()

sink(file.path(paste0(opt$output,'.summary.mod.txt')))
summary(mod)
sink()

pred.combo.train <- pred.base[ind.train] + pred.mod.train
pred.combo.test <- pred.base[ind.test] + pred.mod.test
saveRDS(list(pred.base=pred.base,pred.mod.train=pred.mod.train,pred.mod.test=pred.mod.test,pred.combo.train=pred.combo.train,pred.combo.test=pred.combo.test), file=file.path(paste0(opt$output, ".pred.rds")))

writeLines(sprintf(
  "R^2 scores:::\nBaseline: ---\nTrain: %f\nTest: %f\nPLR: ---\nTrain: %f\nTest: %f\nCombo: ---\nTrain: %f\nTest: %f",
  cor(pred.base[ind.train], data.train[[opt$pheno_col]])^2, # Train Baseline
  cor(pred.base[ind.test], data.test[[opt$pheno_col]])^2, # Test Baseline
  cor(pred.mod.train, data.train[[opt$pheno_col]])^2, # Train PLR
  cor(pred.mod.test, data.test[[opt$pheno_col]])^2, # Test PLR
  cor(pred.combo.train, data.train[[opt$pheno_col]])^2, # Train Combo
  cor(pred.combo.test, data.test[[opt$pheno_col]])^2 # Test Combo
), file.path(paste0(opt$output,'.scores.txt')))
