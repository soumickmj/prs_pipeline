source("scripts/utils.R")

sessionInfo()

options(warn = -1)

option_list = list(
  make_option(c('-i', '--input'), action='store', type='character', help='.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]', default='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/cond_plus_plink_maf1p_geno10p_caucasian_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.bgen'),
  make_option(c('-c', '--covariates'), action='store', type='character', help='.tsv file with covariates (first 2 columns could be FID and IID, 3rd column onwards will be considered as covariates.)', default = '/group/glastonbury/soumick/PRS/inputs/cov_nonMRI_cohort.tsv'),
  make_option(c('--cov_cols'), action='store', type='character', help="Coma-seperated list of covaraite columns. If left blank, all columns will be used.", default='Genotype_batch,Sex,Age,Ethnic_gen,BSA'),
  make_option(c('--cov_cat'), action='store', type='character', help="Coma-seperated list of categorical covaraite columns. If left blank, all columns will be considered continous.", default='Genotype_batch,Sex,Ethnic_gen'),
  make_option(c('--cov_logical'), action='store', type='character', help="Coma-seperated list of logical covaraite columns. If left blank, only cov_cols and cov_cat will be considered."),
  make_option(c('--cov_PC'), action='store', type='character', help=".rds file containing the principal components of the genotype. If supplied, this will be used during the training of the sparse linear regression model.", default = "/project/ukbblatent/clinicaldata/v1.1.0_seventh_basket/genPC_82779_MD_01_03_2024_00_05_30.tsv"),
  make_option(c('--cov_nPC'), action='store', type='numeric', help='Number of PCs to include from the supplied cov_PC (Default: 20, same as the PLR code]', default=20),
  make_option(c('--subs2include'), action='store', type='character', help='txt file (or, output of plink) with FID and IID columns (tab-separated), containing subjects to include (typically used for relatedness filtering)'),
  make_option(c('--rsids2include'), action='store', type='character', help='txt file (or, output of plink) containing a list of rsIDs'),
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]', default="/group/glastonbury/soumick/PRS/LDPred2/initial_nonMRI_basic_FiltAF_S1701_Z49"),
  make_option(c('--ext_sumstats'), action='store', type='character', help='Path to external sumstats. In this case, GWAS will not be performed and the top SNPs from the sumstats will be considered', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/results/gwas/S1701_Z49.gwas.regenie.gz'),
  make_option(c('--ext_col_sumstats'), action='store', type='character', help='.json file defining columns to use from the extarnal sumstats', default = '/home/soumick.chatterjee/Codes/GitLab/tricorder/PRS/davide_fede/sumcols_UKBB_regenie.json'),
  make_option(c('--filtMAF'), action='store', type='numeric', help='[Only if filtered is 1 and ext_sumstats is supplied] Filter sumstats to remove SNPs with MAF < 0.01', default=1),
  make_option(c('--threads'), action='store', type='numeric', help='computing threads [1]', default=5)
  )

opt = parse_args(OptionParser(option_list=option_list))

now<-Sys.time()
message('[',now,'][Message] performing pre-flight checks') 
is_bgen<-check_input(opt$input)
check_output(opt)
message('[',now,'][Message] done') 

cov<-data.frame(fread(opt$covariates))
if (!is.null(opt$cov_cols)){
  message('[',now,'][Message] cov_cols supplied, selecting only the provided columns from the covariates table')
  cov <- cov[c(names(cov)[1:2], strsplit(opt$cov_cols, ",")[[1]])]
}
cov_cols <- colnames(cov)[-c(1, 2)]

merged_data <- cov
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

if (!is.null(opt$rsids2include)){
  message('[',now,'][Message] rsids2include supplied, analysis will be limited to the provided rsIDs')
  filt_rsIDs <- readLines(opt$rsids2include)
} else {
  filt_rsIDs <- NULL
}

#Read and prepare genotype data
now<-Sys.time()
message('[',now,'][Message] loading data') 
message('[',now,'][Message] reading  .bed/.bgen')

if (!is_bgen) {
  
  obj.bigSNP<-load_bed(opt$input, threads=opt$threads)
  
} else {
  
  obj_bgen<-load_bgen(opt$input,threads=opt$threads, subset_subIDs=merged_data$IID, rsIDs=filt_rsIDs) #subsetting the genotypes to keep only the subjects present in the phenotypes table
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

now<-Sys.time()
message('[',now,'][Message] loading sumstats.')
stats<-load_summary(opt$ext_sumstats, opt$ext_col_sumstats, opt$threads)
sumstats<-stats[[1]]

matched_indices <- match(rsIDs, sumstats$rsid)
sumstats_ordered <- sumstats[na.omit(matched_indices), ]
      
if (opt$filtMAF == 1){
  
  sumstats_ordered$maf <- pmin(sumstats_ordered$a1freq, 1 - sumstats_ordered$a1freq)
  sumstats_ordered <- subset(sumstats_ordered, maf > 0.01)
  
} 

sumstats_ordered$key <- paste(sumstats_ordered$chr, sumstats_ordered$pos, sep = "-")

now<-Sys.time()
message('[',now,'][Message] sumstats loaded.')

#Calculate LDMatrix
message('[',now,'][Message] calculating the LDMatrix.')
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
fam.order <- NULL
map <- obj.bigSNP$map[c("chromosome", "rsid", "physical.pos", "allele1", "allele2")]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
map$chr <- as.integer(map$chr)
info_snp <- snp_match(sumstats_ordered, map, join_by_pos = FALSE)
genotype <- obj.bigSNP$genotypes
CHR <- map$chr
POS <- map$pos
for (chr in 1:22) { 
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = opt$threads,
    infos.pos = POS[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))
# saveRDS(ld, file=file.path(paste0(opt$output, ".LD.LDPred2.rds"))) 
# saveRDS(corr, file=file.path(paste0(opt$output, ".corr.LDPred2.rds")))
now<-Sys.time()
message('[',now,'][Message] LDMatrix ready.')

#Perform LD score regression
message('[',now,'][Message] performing LD score regression.')
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld),
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff,
                    blocks = NULL)
# saveRDS(ldsc, file=file.path(paste0(opt$output, ".LDSC.LDPred2.rds")))
h2_est <- ldsc[["h2"]]
now<-Sys.time()
message('[',now,'][Message] done.')

if (h2_est <= 0) {
  message('[',now,'][Message] resulted in negative heritability. LDPred2 cannot work with it. So, replacing it with 0.01.')
  h2_est <- 0.01
}


now<-Sys.time()
message('[',now,'][Message] done.')


#infinitesimal model
message('[',now,'][Message] infinitesimal model starting...')
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
genotype <- obj.bigSNP$genotypes
ind.test <- 1:nrow(genotype)
pred_inf <- big_prodVec(    genotype,
                            beta_inf,
                            ind.row = ind.test,
                            ind.col = info_snp$`_NUM_ID_`)
saveRDS(list(beta_inf=beta_inf, pred_inf=pred_inf), file=file.path(paste0(opt$output, ".fullDS.inf.mod.LDPred2.rds")))
now<-Sys.time()
message('[',now,'][Message] done.')


#grid model
message('[',now,'][Message] grid model starting....')
p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
grid.param <- expand.grid(p = p_seq,
              h2 = h2_seq,
              sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = opt$threads)
genotype <- obj.bigSNP$genotypes
ind.test <- 1:nrow(genotype)
pred_grid <- big_prodMat(   genotype,
                            beta_grid,
                            ind.col = info_snp$`_NUM_ID_`)
saveRDS(list(beta_grid=beta_grid, pred_grid=pred_grid), file=file.path(paste0(opt$output, ".fullDS.grid.mod.LDPred2.rds")))
now<-Sys.time()
message('[',now,'][Message] done.')

#auto model
message('[',now,'][Message] auto model starting....')
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = opt$threads),
  ncores = opt$threads
)
beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)
genotype <- obj.bigSNP$genotypes
ind.test <- 1:nrow(genotype)
pred_auto <-
  big_prodMat(genotype,
              beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
  rowMeans(beta_auto[,
                     abs(pred_scaled -
                           median(pred_scaled)) <
                       3 * mad(pred_scaled)])
pred_auto <-
  big_prodVec(genotype,
              final_beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)
saveRDS(list(multi_auto=multi_auto, beta_auto=beta_auto, pred_auto=pred_auto, pred_scaled=pred_scaled, final_beta_auto=final_beta_auto, pred_auto=pred_auto), file=file.path(paste0(opt$output, ".fullDS.auto.mod.LDPred2.rds")))
now<-Sys.time()
message('[',now,'][Message] done.')

#final save
saveRDS(data.full, file=file.path(paste0(opt$output, ".fullDS.resNdata.basic.LDPred2.rds")))

now<-Sys.time()
message('[',now,'][Message] done')
