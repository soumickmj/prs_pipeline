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
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]', default="/group/glastonbury/soumick/PRS/simple/initial_test_FreqFilt_simple_S1701_Z49"),
  make_option(c('--selectSNPs'), action='store', type='character', help='[Only if filtered is 1] specify path to a .csv or .tsv file containing the CHR and POS of the selected SNPs. If not supplied, GWAS-driven filtering will be performed.', default='/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/results/gwas/independent/genome-wide_significant_hits_post_cojo.csv'),
  make_option(c('--col_selectSNPs'), action='store', type='character', help='[Only if filtered is 1 and selectSNPs is supplied] Coma sperated list of column names in the selectSNPs file containing CHR and POS', default = 'Chr,bp'),
  make_option(c('--ext_sumstats'), action='store', type='character', help='Path to external sumstats. In this case, GWAS will not be performed and the top SNPs from the sumstats will be considered', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/results/gwas/S1701_Z49.gwas.regenie.gz'),
  make_option(c('--ext_col_sumstats'), action='store', type='character', help='.json file defining columns to use from the extarnal sumstats', default = '/home/soumick.chatterjee/Codes/GitLab/tricorder/PRS/davide_fede/sumcols_UKBB_regenie.json'),
  make_option(c('--filtAF'), action='store', type='numeric', help='[Only if filtered is 1 and ext_sumstats is supplied] Filter sumstats to remove SNPs with MAF < 0.01', default=1),
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

rsIDs <- obj.bigSNP$map$rsid
CHR <- as.integer(obj.bigSNP$map$chromosome) #this is somehow necessary for .bgen files, not for bed
POS <- obj.bigSNP$map$physical.pos

now<-Sys.time()
message('[',now,'][Message] genotypes loaded and processed')

#read the sample and match 
merged_data <- subset(merged_data, IID %in% obj.sample$ID_2) #discarding phenotypes and covariates for the subjects not present in the genotype data
data.full <- merged_data[match(obj.sample$ID_2, merged_data$IID), ] #ensure that both the genotype and phenotype tables are in the same order

if (tolower(file_ext(opt$selectSNPs)) == "csv") {
  selectSNPs <- read.csv(opt$selectSNPs)
}  else {
  selectSNPs <- read.table(opt$selectSNPs)
}
cols_selectSNPs <- strsplit(opt$col_selectSNPs, ",")[[1]]
selectSNPs_keys <- paste(selectSNPs[[cols_selectSNPs[1]]], selectSNPs[[cols_selectSNPs[2]]], sep = "-")

message('[',now,'][Message] external sumstats supplied, loading the sumstats. GWAS will not be performed.')
stats<-load_summary(opt$ext_sumstats, opt$ext_col_sumstats, opt$threads)
sumstats<-stats[[1]]

matched_indices <- match(rsIDs, sumstats$rsid)
sumstats_ordered <- sumstats[na.omit(matched_indices), ]
      
if (opt$filtAF == 1){
  
  sumstats_ordered <- subset(sumstats_ordered, a1freq > 0.01)
  include_SNPs <- na.omit(match(sumstats_ordered$rsid, rsIDs))
  
} else {
  
  include_SNPs <- 1:nrow(sumstats_ordered)

}

sumstats_ordered$key <- paste(sumstats_ordered$chr, sumstats_ordered$pos, sep = "-")


#PRS
prs_all <- rep(0, nrow(data.full))
prs_sig <- rep(0, nrow(data.full))
prs_select <- rep(0, nrow(data.full))

pb <- txtProgressBar(min = 0, max = nrow(sumstats_ordered), style = 3)
for(i in 1:nrow(sumstats_ordered)) {
  snp_info <- sumstats_ordered[i, ]
  ind_G <- include_SNPs[i]
  dosage <- G[, ind_G]  
  snp_contribution <- dosage * snp_info$beta
  prs_all <- prs_all + snp_contribution
  
  if(snp_info$p < 5e-8) {
    prs_sig <- prs_sig + snp_contribution
  }
  
  if(snp_info$key %in% selectSNPs_keys) {
    prs_select <- prs_select + snp_contribution
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

data.full$PRS_all <- prs_all
data.full$prs_sig <- prs_sig
data.full$prs_select <- prs_select
  
saveRDS(data.full, file=file.path(paste0(opt$output, ".resNdata.simplePRS.rds")))

get_R2 <- function(pred) {
  mylm <- lm(y ~ pred + COVAR, data.frame(pred = pred, y = data.full[[opt$pheno_col]], COVAR = I(covar_from_df(data.full[, c(cov_cols, cov_PC_cols)]))))
  summary(mylm)$r.squared
}

cor_PRS_all <- get_R2(data.full$PRS_all)
cor_PRS_sig <- get_R2(data.full$prs_sig)
cor_PRS_select <- get_R2(data.full$prs_select)

output <- capture.output(paste("R^2:--- PRS_all: ", cor_PRS_all, "PRS_sig: ", cor_PRS_sig, "PRS_select: ", cor_PRS_select))
writeLines(output, file.path(paste0(opt$output, ".simplePRS.fullDS.results.txt")))

now<-Sys.time()
message('[',now,'][Message] done')
