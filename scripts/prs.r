#!/usr/bin/env Rscript

source("scripts/utils.R")

sessionInfo()

# run functions before running the pipeline

check_model<-function(x) {

  if (!x %in% c('automatic', 'grid')) {

    now<-Sys.time()
    stop(red('[',now,'][Error] specified model can be either automatic or grid'))

  }

}

check_phenotype<-function(x) {

  if (is.null(x$phenotype)) {

    now<-Sys.time()
    message(yellow('[',now,'][Warning] missing table of phenotypes/covariates. Model coherced to automatic even when grid is specified'))
    x$model<-"automatic"

    #if (! is.null(x$train)) {

      #now<-Sys.time()
      #stop('[',now,'][Error] cannot peform train/test evaluation with automatic model') 
    
    #}

  } else {

    if (! file.exists(file.path(x$phenotype))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] if provided, table of phenotypes/covariates must exist'))

    }
  }

  x$model
}

check_predictions<-function(x,pheno) {

  sub_pheno<-pheno[,c(3:ncol(pheno))]

  if (class(sub_pheno) == "numeric") { #this is a vector

    sub_pheno<-cbind(sub_pheno)
    colnames(sub_pheno)<-"V3"

  }

  sub_val<-data.frame(cbind(sub_pheno,x))
  params<-list()

  if (length(unique(sub_val$V3)) == 2) { #binary trait  #maybe change to length(unique(subpheno$V3))

    res<-summary(glm(V3 ~ ., data=sub_val, family="binomial"))
    params$score<-res$coef["x",3]
    params$prt<-res$coef["x",4]
    params$r2<-1-(res$deviance/res$null.deviance) #calculate (adjusted) r2

  } else {
    
    res<-summary(lm(V3 ~ ., data=sub_val))
    params$score<-res$coef["x",3]
    params$prt<-res$coef["x",4]
    params$r2<-res$adj.r.squared
  
  }

  list(res,params)

}

get_complement<-function(nuc) {

  x<-"A"
  if (nuc=="A") x<-"T"
  if (nuc=="T") x<-"A"
  if (nuc=="C") x<-"G"
  if (nuc=="G") x<-"C"
  x

}

check_correlation<-function(x) {

  has_corr<-TRUE

  if (is.null(x)) {

    has_corr<-FALSE

  } else {

    if (! file.exists(file.path(x))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] if provided, correlation matrix must exist'))

    } else {

        if (file_ext(file.path(x)) != "rds") {

            now<-Sys.time()
            stop(red('[',now,'][Error] if provided, correlation matrix must be provided in .rds file'))
        }

        if (! file.exists(str_replace(file.path(x), ".rds", ".sbk"))) {

            now<-Sys.time()
            stop(red('[',now,'][Error] if provided, correlation matrix must have a supporting .sbk file - same name, different extension - in the same directory'))

        }

        if (! file.exists(str_replace(file.path(x), ".corr.rds", ".ld.rds"))) {

            now<-Sys.time()
            stop(red('[',now,'][Error] if provided, correlation matrix must have a matchind .ld.rds file - different name, same extension - in the same directory'))

        }

    }

  }

  has_corr
}

# start the prs script 
# create the list with all necessary files
options(warn = -1)

# option_list = list(
#   make_option(c('-m', '--model'), action='store', type='character', help='model to use, either automatic or grid. Default to automatic', default='automatic'),
#   make_option(c('-i', '--input'), action='store', type='character', help='.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]'),
#   make_option(c('-s', '--summary'), action='store', type='character', help='GWAS summary stats or pre-calculated .tsv (with header) containing beta scores [required]'),
#   make_option(c('--summarycols'), action='store', type='character', help='.json file defining columns to use'),
#   make_option(c('--filter_SNPs'), action='store', type='character', help='.txt file containing a list of SNPs to be used to filter the summary stats and used in the analysis'),
#   make_option(c('-p', '--phenotype'), action='store', type='character', help='.tsv file with phenotype (first 2 columns could be FID and IID, 3rd column should be the phenotype. If more than 1 phenotype is present, specify the column to be considered using pheno_col)'),
#   make_option(c('--pheno_col'), action='store', type='character', help='If phenotype is provided, but contains more than 1 phenotype, the phenotype that is to be used can be specified here.'),
#   make_option(c('-c', '--covariates'), action='store', type='character', help='.tsv file with covariates (first 2 columns could be FID and IID, 3rd column onwards will be considered as covariates.) If both phenotype and covariates are supplied, first 2 columns must match.'),
#   make_option(c('--heritability'), action='store', type='numeric', help='heritability, if known'),
#   make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]'),
#   make_option(c('--threads'), action='store', type='numeric', help='computing threads [1]', default=1),
#   make_option(c('--correlation'), action='store', type='character', help='the correlation matrix provided as a pre-computed .rds object'),
#   make_option(c('--reference'), action='store', type='character', help='reference panel in .rds format'),
#   make_option(c('--index'), action='store', type='character', help='indexes for the individuals in the reference panel in .tsv format')
#   #make_option(c('-t', '--train'), action='store', type='character', help='train percentage for training-testing - internal validation'),
#   #make_option(c('-x', '--external'), action='store', type='character', help='external validation set. Comma-separated .bed (or .bgen) and .pheno tsv. .bed should have accompanying .bim and .fam')
# )

option_list = list(
  make_option(c('-m', '--model'), action='store', type='character', help='model to use, either automatic or grid. Default to automatic', default='automatic'),
  make_option(c('-i', '--input'), action='store', type='character', help='.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]', default='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/conv_cond_plus_plink_maf1p_geno10p_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.bgen'),
  make_option(c('-s', '--summary'), action='store', type='character', help='GWAS summary stats or pre-calculated .tsv (with header) containing beta scores [required]', default='/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/results/gwas/S1701_Z49.gwas.regenie.gz'),
  make_option(c('--summarycols'), action='store', type='character', help='.json file defining columns to use', default='/home/soumick.chatterjee/Codes/GitLab/tricorder/PRS/davide_fede/sumcols_UKBB_regenie.json'),
  make_option(c('--filter_SNPs'), action='store', type='character', help='.txt file containing a list of SNPs to be used to filter the summary stats and used in the analysis', default ='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/cond_plus_plink_maf1p_geno10p_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.txt'),
  make_option(c('-p', '--phenotype'), action='store', type='character', help='.tsv file with phenotype (first 2 columns could be FID and IID, 3rd column should be the phenotype. If more than 1 phenotype is present, specify the column to be considered using pheno_col)', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/validated_input/processed.pheno.validated.txt'),
  make_option(c('--pheno_col'), action='store', type='character', help='If phenotype is provided, but contains more than 1 phenotype, the phenotype that is to be used can be specified here.', default='S1701_Z49'),
  make_option(c('-c', '--covariates'), action='store', type='character', help='.tsv file with covariates (first 2 columns could be FID and IID, 3rd column onwards will be considered as covariates.)', default = '/group/glastonbury/GWAS/F20208v3_DiffAE/select_latents_r80/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/nNs_Qntl_INF30_DiffAE128_5Sd_r80_discov_fullDSV3/validated_input/cov_newset_chp_F20208_Long_axis_heart_images_DICOM_H5v3_NOnoise.cov.validated.txt'),
  make_option(c('--heritability'), action='store', type='numeric', help='heritability, if known'),
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]', default='/group/glastonbury/soumick/PRS/davide_fede/initial_test_S1701_Z49_ldprune_prova'),
  make_option(c('--threads'), action='store', type='numeric', help='computing threads [1]', default=10),
  make_option(c('--correlation'), action='store', type='character', help='the correlation matrix provided as a pre-computed .rds object'),
  make_option(c('--reference'), action='store', type='character', help='reference panel in .rds format'),
  make_option(c('--index'), action='store', type='character', help='indexes for the individuals in the reference panel in .tsv format')
  #make_option(c('-t', '--train'), action='store', type='character', help='train percentage for training-testing - internal validation'),
  #make_option(c('-x', '--external'), action='store', type='character', help='external validation set. Comma-separated .bed (or .bgen) and .pheno tsv. .bed should have accompanying .bim and .fam')
)

opt = parse_args(OptionParser(option_list=option_list))

#quickly check all input. This functions do not load anything in r

now<-Sys.time()
message('[',now,'][Message] performing pre-flight checks') 

check_model(opt$model)
is_bgen<-check_input(opt$input)
check_summary(opt$summary)
opt$model<-check_phenotype(opt)
check_output(opt)
has_corr<-check_correlation(opt$correlation)
#train_test<-check_train(opt$train)
#external_validation<-check_external(opt$external)
has_index<-check_index(opt$reference)
has_index2<-check_index(opt$index)

now<-Sys.time()
message('[',now,'][Message] done') 
message('[',now,'][Message] loading data') 
message('[',now,'][Message] reading  .bed/.bgen')

if (!is_bgen) {

  obj.bigSNP<-load_bed(opt$input, threads=opt$threads)

} else {

  obj_bgen<-load_bgen(opt$input,threads=opt$threads)
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
message('[',now,'][Message] done')

if (!is.null(opt$phenotype)) {

  message('[',now,'][Message] reading table of phenotypes')

  pheno<-data.frame(fread(opt$phenotype))
  
  if (!is.null(opt$pheno_col)) {
    message('[',now,'][Message] pheno_col has been supplied. Discarding all the columns except for the first two and pheno_col')
    pheno_col_index <- which(names(pheno) == opt$pheno_col)
    pheno <- pheno[, c(1:2, pheno_col_index)]
  }

  now<-Sys.time()
  message('[',now,'][Message] done')
    
}

if (!is.null(opt$covariates)) {
  message('[',now,'][Message] reading table of covariates')
  
  cov<-data.frame(fread(opt$covariates))
  
  if (exists('pheno')) {
    message('[',now,'][Message] phenotype table is also supplied and has been already read.')
    message('[',now,'][Message] phenotype table has ', nrow(pheno), ' elements, while covariates table has ',nrow(cov),' elements.')
    pheno <- merge(pheno, cov, by = names(pheno)[1:2])
    message('[',now,'][Message] After merging phenotype and covariates tables, there are ', nrow(pheno), ' common elements.')
  } else {
    message('[',now,'][Message] phenotype table is not supplied. Considering only the covariates table.')
    pheno <- cov
  }
  
  now<-Sys.time()
  message('[',now,'][Message] done')
}

if (exists('pheno')) {
  message('[',now,'][Message] read phenotype and/or covariates table(s), and now changing the column names to make them consistent.')
  colnames(pheno)<-paste0("V", c(1:ncol(pheno))) #be sure to have proper name for pheno columns
  
  now<-Sys.time()
  message('[',now,'][Message] done')
}

#not tested yet
#if (train_test) {

  #now<-Sys.time()
  #message('[',now,'][Warning] provided a percentage for training. Model coherced to grid even when automatic is specified')
  #opt$model<-"grid" 
  #n_val<-round(nrow(G)/100*as.numeric(opt$train))
  #ind.val<-sample(nrow(G), n_val)
  #ind.test<-setdiff(rows_along(G), ind.val)
  
#} else {

ind.test<-ind.val<-c(1:nrow(G))

#}

message('[',now,'][Message] reading summary statistics')

stats<-load_summary(opt$summary, opt$summarycols, opt$threads)
sumstats<-stats[[1]]
beta_is_precomp<-stats[[2]]

if (!is.null(opt$filter_SNPs)) {
  message('[',now,'][Message] filter_SNPs has been supplied. Hence, the rsIDs from this file will be read and the sumstats will be filtered to only keep the provided SNPs.')
  rsIDs <- readLines(opt$filter_SNPs)
  rsIDs <- trimws(rsIDs) 
  rsIDs <- rsIDs[rsIDs != ""] 
  sumstats <- sumstats[sumstats$rsid %in% rsIDs, ]
  message('[',now,'][Message] sumstats filtering based on the provided list of SNPs: done')
}

sumstats$rsid <- sub("_.*", "", sumstats$rsid) #remove the alleles have they are present in the sumstats

now<-Sys.time()
message('[',now,'][Message] done')
message('[',now,'][Message] matching variants between genotype data and summary statistics - or previously computed beta scores')

if (is_bgen) {

  map<-data.frame(obj.bigSNP$map)
  map$chromosome<-as.numeric(map$chromosome) #this is character otherwise
  map<-map[c("chromosome", "rsid", "physical.pos", "allele2", "allele1")]

} else { #is .bed

  map<-obj.bigSNP$map
  map<-map[c("chromosome", "marker.ID", "physical.pos", "allele1", "allele2")]

}

colnames(map)<-c("chr", "rsid", "pos", "a1", "a0")

#take care of incorrect rsid (?) how frequently does this happen?

rsid_map<-map$rsid
idsl<-strsplit(rsid_map, "_")
ids_<-sapply(idsl,"[[",1)
map$rsid<-ids_

##filter based on panel, first

sumstats<-recover_missing_cols(sumstats,map)
sumstats$chr<-as.integer(sumstats$chr)

df_beta<- snp_match(sumstats, map,join_by_pos = FALSE) #match by rsid for the time being

now<-Sys.time()
message('[',now,'][Message] done')

if (has_index) {

  now<-Sys.time()
  message('[',now,'][Message] subsetting to variants in the reference panel provided')

  reference_ext <- file_ext(opt$reference)

  if (reference_ext == "rds") {
    map_ldref <- readRDS(opt$reference)
    message("[", Sys.time(), "][Message] Loaded .rds reference panel")
  } 
  else if (reference_ext == "bed") {
    temp_rds_file <- tempfile(fileext = ".rds") # Temporary file to hold the .rds conversion
    convert_bed_to_rds(opt$reference, temp_rds_file, threads=opt$threads)
    message("[", Sys.time(), "][Message] Converted .bed reference panel to temporary .rds file")
    map_ldref <- readRDS(temp_rds_file)
    message("[", Sys.time(), "][Message] Loaded temporary .rds reference panel")
  }

  map_panel<-map_ldref$map
  map_panel<-map_panel[c("chromosome", "marker.ID", "physical.pos", "allele1", "allele2")]
  colnames(map_panel)<-c("chr", "rsid", "pos", "a1", "a0")
  keep<-which(df_beta$rsid %in%map_panel$rsid)
  df_beta<-df_beta[keep,]
  now<-Sys.time()
  message('[',now,'][Message] ', nrow(df_beta), ' variants will be used')
  message('[',now,'][Message] done')

}

if (!beta_is_precomp) {

  if (!has_corr) { 

    POS2 <- snp_asGeneticPos(CHR, POS, dir = dirname(opt$output), ncores = opt$threads)
    
    if (has_index2) {

      now<-Sys.time()
      message('[',now,'][Message] loading indexes to subset correlation matrix')
      individual_idx<-as.integer(fread(opt$index)$V1)
      now<-Sys.time()
      message('[',now,'][Message] done')

    } else {
    
      individual_idx<-c(1:nrow(G))

    }

    now<-Sys.time()
    message('[',now,'][Message] computing correlation between variants')

    for (chr in 1:22) {

        # print(chr)
        ind.chr <- which(df_beta$chr == chr) #subsample G

        if (length(ind.chr) == 0) next

        ## indices in 'G'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, ind.row=individual_idx, infos.pos = POS2[ind.chr2], ncores = opt$threads)

        if ((chr == 1) | (length(unique(df_beta$chr)) == 1)) { #also if we have a single chromosome
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, backingfile=file.path(paste0(opt$output, ".corr")), compact = TRUE)
        } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
        }
    }

    now<-Sys.time()
    message('[',now,'][Message] done')
    message('[',now,'][Message] storing correlation matrix and ld values to files')

    saveRDS(corr, file=file.path(paste0(opt$output, ".corr.rds")))
    saveRDS(ld, file=file.path(paste0(opt$output, ".ld.rds")))

  } else {

    now<-Sys.time()
    message('[',now,'][Message] loading precomputed correlation values')
    corr<-readRDS(file.path(opt$correlation))
    ld<-readRDS(file.path(str_replace(opt$correlation, ".corr.rds", ".ld.rds")))
  
  }

  now<-Sys.time()
  message('[',now,'][Message] done')

}


if (!beta_is_precomp) {

	h2_tab_file<-file.path(paste0(opt$output, ".h2.tsv"))

	if (is.null(opt$heritability)) {

            now<-Sys.time()
            message('[',now,'][Message] estimating h2')

            ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL))
            h2_est <- ldsc[["h2"]]

            now<-Sys.time()
            message('[',now,'][Message] done')

        } else {

            h2_est<-opt$heritability

        }

        h2_tab<-cbind(h2_est)
        colnames(h2_tab)<-"h2"
        fwrite(data.frame(h2_tab), file=h2_tab_file, sep="\t", col.names=TRUE, row.names=F)

}


if (!beta_is_precomp) {

  if (opt$model == "automatic") {

    beta_tab_file<-file.path(paste0(opt$output,'.auto.beta_scores.tsv'))
    prs_tab_file<-file.path(paste0(opt$output,'.auto.prs.tsv'))
    
    summary_file<-file.path(paste0(opt$output,'.auto.summary.rds')) #use only if pheno provided
    summarytab_file<-file.path(paste0(opt$output,'.auto.summary.tsv')) #use only if pheno provided

    auto_stats_file<-file.path(paste0(opt$output, ".auto.params.rds"))

    #summary_file2<-file.path(paste0(opt$output,'.auto.summary.external.rds')) #use only if ext
    #summarytab_file2<-file.path(paste0(opt$output,'.auto.summary.external.tsv')) #use only if ext

    now<-Sys.time()
    message('[',now,'][Message] running automatic model')

    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                  vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                  allow_jump_sign = FALSE, shrink_corr = 0.95,
                  ncores = opt$threads)
    saveRDS(multi_auto, file=auto_stats_file)
    range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
    keep <- (range > (0.9 * quantile(range, 0.9)))
    beta<-rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    pred <- big_prodVec(G, beta, ind.col = df_beta[["_NUM_ID_"]])

  } else { #is grid

    beta_tab_file<-file.path(paste0(opt$output,'.grid.beta_scores.tsv'))
    prs_tab_file<-file.path(paste0(opt$output,'.grid.prs.tsv'))
    
    #summary_file<-file.path(paste0(opt$output,'.grid.summary.rds'))
    #summarytab_file<-file.path(paste0(opt$output,'.grid.summary.tsv'))
    
    rank_file<-file.path(paste0(opt$output,'.grid.allsummary.tsv'))
    best_model_file<-file.path(paste0(opt$output,'.grid.bestsummary.rds'))

    #summary_file2<-file.path(paste0(opt$output,'.grid.summary.external.rds')) #use only if ext
    #summarytab_file2<-file.path(paste0(opt$output,'.grid.summary.external.tsv')) #use only if ext

    now<-Sys.time()
    message('[',now,'][Message] running grid model')

    h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
    p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
    params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = opt$threads)
    pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
    sub_pheno<-pheno[ind.val,c(3:ncol(pheno))] #it may be equal to pheno if no split

    if (length(unique(pheno$V3)) == 2) { #binary trait 

      res <- apply(pred_grid[ind.val, ], 2, function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          sub_val<-cbind(sub_pheno,x)
          colnames(sub_val)[1]<-"V3"
          summary(glm(V3 ~ ., data=data.frame(sub_val), family="binomial"))
        }
      })

      params$score<-do.call(c,lapply(res, function(x) x$coef["x",3]))
      params$prt<-do.call(c,lapply(res, function(x) x$coef["x",4]))
      params$r2<-1-(do.call(c,lapply(res, function(x) x$deviance)/do.call(c,lapply(res, function(x) x$null.deviance)))) #calculate (adjusted) r2

    } else { #quantitative/continuous trait

      res <- apply(pred_grid[ind.val, ], 2, function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          sub_val<-cbind(sub_pheno,x)
          colnames(sub_val)[1]<-"V3"
          summary(lm(V3 ~ ., data=data.frame(sub_val)))
        }
      })

      params$score<-do.call(c,lapply(res, function(x) x$coef["x",3]))
      params$prt<-do.call(c,lapply(res, function(x) x$coef["x",4]))
      params$r2<-do.call(c,lapply(res, function(x) x$adj.r.squared))

    }

    rank<- params %>%
      mutate(id = row_number()) %>%
      # filter(sparse) %>% not used. But leave here for future reference  
      arrange(desc(score))

    rank$snps<-nrow(df_beta)

    best_beta_grid_id<- rank %>%
      slice(1) %>%
      pull(id)    

    summary_best_beta_grid<-res[[best_beta_grid_id]]

    now<-Sys.time()
    message('[',now,'][Message] storing tested and best model summary to .tsv and .rds object')

    fwrite(rank, file=rank_file, col.names=T, row.names=F, sep="\t")
    saveRDS(summary_best_beta_grid, file = best_model_file)
    
    beta<-beta_grid[,best_beta_grid_id]

    #if (train_test) {

      #pred_test<-big_prodVec(G, beta, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
      #pred_train<-big_prodVec(G, beta, ind.row = ind.val, ind.col = df_beta[["_NUM_ID_"]])
      #pred<-c(pred_train,pred_test)
      #type_pred<-c(rep("TRAIN", length(pred_train)), rep("TEST", length(pred_test)))     
    
    #} else { #no validation requested or file with external validation given

      pred<-big_prodVec(G, beta, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]]) #here ind.test==ind.val

    #}

  }

} else { #betas as pre-computed values

  now<-Sys.time()
  message('[',now,'][Message] using pre-computed beta scores')

  beta_tab_file<-file.path(paste0(opt$output,'.precomp.beta_scores.tsv'))
  prs_tab_file<-file.path(paste0(opt$output,'.precomp.prs.tsv'))
  
  summary_file<-file.path(paste0(opt$output,'.precomp.summary.rds')) #use only if pheno provided
  summarytab_file<-file.path(paste0(opt$output,'.precomp.summary.tsv')) #use only if pheno provided
  
  #summary_file2<-file.path(paste0(opt$output,'.precomp.summary.external.rds')) #use only if ext
  #summarytab_file2<-file.path(paste0(opt$output,'.precomp.summary.external.tsv')) #use only if ext

  beta<-df_beta$beta
  pred <- big_prodVec(G, beta, ind.col = df_beta[["_NUM_ID_"]])

}

beta_tab<-data.frame(chr=df_beta$chr, rsid=df_beta$rsid,  pos=df_beta$pos, a0=df_beta$a0, a1=df_beta$a1, beta=beta)

#if (train_test) {

  #prs_tab<-data.frame(FID=obj.bigSNP$fam$family.ID,IID=obj.bigSNP$fam$sample.ID,PRED=pred, TYPE=type_pred)

#} else {

prs_tab<-data.frame(FID=obj.bigSNP$fam$family.ID,IID=obj.bigSNP$fam$sample.ID,PRED=pred, TYPE="ALL")

#}

now<-Sys.time()
message('[',now,'][Message] done')

now<-Sys.time()
message('[',now,'][Message] storing beta scores and predictions to .tsv files')

fwrite(beta_tab, file=beta_tab_file,col.names=T, row.names=F, sep="\t")
fwrite(prs_tab, file=prs_tab_file, col.names=T, row.names=F, sep="\t")

now<-Sys.time()
message('[',now,'][Message] done')

if (exists('pheno') & (opt$model != "grid" | beta_is_precomp)) {
  pheno_prs_tab <- merge(pheno, prs_tab[, c("FID", "IID", "PRED")], by.x = c("V1", "V2"), by.y = c("FID", "IID"))
  pred <- pheno_prs_tab$PRED
  pheno <- pheno_prs_tab[, -ncol(pheno_prs_tab)]

  now<-Sys.time()
  message('[',now,'][Message] evaluating model on provided phenotype')

  res_params<-check_predictions(pred,pheno)
  res<-res_params[[1]]
  params<-res_params[[2]]
  pheno_tab<-data.frame(score=params$score, prt=params$prt, r2=params$r2, snps=nrow(df_beta))

  now<-Sys.time()
  message('[',now,'][Message] done')

  now<-Sys.time()
  message('[',now,'][Message] storing performances to .rds and .tsv files')

  saveRDS(res, file = summary_file)
  fwrite(pheno_tab, file=summarytab_file, sep="\t", col.names=TRUE, row.names=F)

  now<-Sys.time()
  message('[',now,'][Message] done')

}
