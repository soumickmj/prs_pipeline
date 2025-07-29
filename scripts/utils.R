#!/usr/bin/env Rscript

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(data.table)
library(bigreadr)
library(magrittr)
library(stringr)
library(fmsb)
library(dplyr)
library(optparse)
library(tools)
library(rjson)
library(crayon)
library(digest)

check_input<-function(x) {
  
  if (is.null(x)) {
    
    now<-Sys.time()
    stop(red('[',now,'][Error] a .bed/.bgen file must be provided as input'))
    
  } else {
    
    if (! file.exists(file.path(x))) {
      
      now<-Sys.time()
      stop(red('[',now,'][Error] the .bed/.bgen file must exist'))
      
    } else {
      
      if (file_ext(x) == "bed") {
        
        now<-Sys.time()
        message('[',now,'][Message] detected .bed input') 
        
        is_bgen<-FALSE
        bim_f<-str_replace(x, ".bed", ".bim")
        fam_f<-str_replace(x, ".bed", ".fam")
        
        if ((! file.exists(file.path(bim_f))) | (! file.exists(file.path(fam_f)))) {
          
          now<-Sys.time()
          stop(red('[',now,'][Error] .bim and .fam files must exist in .bed file location'))
          
        }
        
      } else if (file_ext(x) == "bgen") {
        
        now<-Sys.time()
        message('[',now,'][Message] detected .bgen input')
        is_bgen<-TRUE
        
        if (!file.exists(paste0(x, ".bgi"))) {
          
          now<-Sys.time()
          stop(red('[',now,'][Error] .bgen file must have a companion .bgi index. Use bgenix with .bgen file to create a proper index'))
        }
        
      } else {
        
        now<-Sys.time()
        stop(red('[',now,'][Error] unknown input file extension'))
        
      }
      
    }
    
  }
  
  is_bgen
  
}


check_output<-function(x) {
  
  if (is.null(x$output)) {
    
    now<-Sys.time()
    stop('[',now,'][Error] output prefix must be specified') 
    
  } else {
    
    #x$output<-normalizePath(file.path(x$output))
    out_dir<-dirname(x$output)
    dir.create(out_dir, recursive=TRUE,showWarning=FALSE)
  }
  
}

check_summary<-function(x) {
  
  if (is.null(x)) {
    
    now<-Sys.time()
    stop(red('[',now,'][Error] summary statistics or tsv with pre-computed beta scores must be provided'))
    
  } else {
    
    if (! file.exists(file.path(x))) {
      
      now<-Sys.time()
      stop(red('[',now,'][Error] summary statistics or tsv with pre-computed beta scores file must exist'))
      
    } 
  } 
}

check_train<-function(x) {
  train_test<-FALSE

  if (!is.null(x)) { #requested train-test

    val_type<-as.numeric(x)

    if (is.na(val_type) | val_type > 100) { #incorrect number

      now<-Sys.time()
      stop('[',now,'][Error] invalid percentage for training') 

    } else { #is valid number

      train_test<-TRUE

    }

  }

  train_test

}

load_summary<-function(x, cols, threads) {
  
  stat_file<-file.path(x)
  cols_file<-file.path(cols)
  jfile<-fromJSON(file = cols_file)
  names_cols<-names(jfile)[c(1:length(jfile)-1)]
  select_values<-as.character(jfile)[c(1:length(jfile)-1)]
  precomp_betas<-type.convert(as.character(jfile)[length(jfile)],as.is=TRUE)
  
  if (class(jfile$n_eff) == "numeric") {
    
    n_val<-jfile$n_eff
    jfile<-jfile[which(names(jfile) != "n_eff")]
    select_values<-as.character(jfile)[c(1:length(jfile)-1)]
    names_cols<-names(jfile)[c(1:length(jfile)-1)]
    
  }
  
  sumstats<-data.frame(fread2(stat_file, nThread=threads, select=select_values))
  colnames(sumstats)<-names_cols
  
  if ("mlog10p" %in% colnames(sumstats)) {
    sumstats$p <- 10^(-sumstats$mlog10p)
  }
  
  if (! "n_eff" %in% colnames(sumstats)) {
    
    sumstats$n_eff<-n_val
    
  }
  
  if (!isTRUE(precomp_betas)) { 
    
    now<-Sys.time()
    message('[',now,'][Message] assuming summary stats file')
    
  } else {
    
    now<-Sys.time()
    message('[',now,'][Message] assuming pre-computed beta scores')
    
  }
  
  list(sumstats,precomp_betas)
  
}

load_bed<-function(x, threads) {
  
  now<-Sys.time()
  message('[',now,'][Message] reading .bed/.bim/.fam files')
  
  bed_file<-file.path(x)
  backing_file<-str_replace(bed_file, ".bed", "") #cut extension out
  bk_file<-file.path(str_replace(bed_file, ".bed", ".bk"))
  rds_file<-file.path(str_replace(bed_file, ".bed", ".rds"))
  
  if (!file.exists(bk_file)) {
    
    snp_readBed2(bed_file, backingfile=backing_file, ncores=threads) #this generate a .rds object that can be loaded into the environment
    
  }
  
  obj.bigSNP <- snp_attach(rds_file)
  obj.bigSNP
  
}


load_bgen<-function(x,threads,subset_subIDs=NULL, rsIDs=NULL) {
  
  now<-Sys.time()
  message('[',now,'][Message] reading .bgen file')
  
  bgen_file<-file.path(x)
  backing_file<-str_replace(bgen_file, ".bgen", "") #cut extension out
  bk_file<-file.path(str_replace(bgen_file, ".bgen", ".bk"))
  rds_file<-file.path(str_replace(bgen_file, ".bgen", ".rds"))
  bgi_file<-file.path(str_replace(bgen_file, ".bgen", ".bgen.bgi"))
  sample_file<-file.path(str_replace(bgen_file, ".bgen", ".sample"))

  if (!is.null(subset_subIDs)) {

    subset_digest <- digest(subset_subIDs)
    if (!is.null(rsIDs)) {
      subset_digest <- paste(subset_digest, digest(rsIDs), sep="_")
    }

    subset_bk_file <- file.path(str_replace(bk_file, ".bk", paste0(".subset.", subset_digest, ".bk")))
    subset_rds_file <- file.path(str_replace(rds_file, ".rds", paste0(".subset.", subset_digest, ".rds")))
    subset_backing_file <- str_replace(subset_bk_file, ".bk", "")

    if (!file.exists(subset_bk_file)) { 
        
      bgi<-snp_readBGI(bgi_file,snp_id=NULL)
      if (!is.null(rsIDs)) {
        message('[',now,'][Message] subsetting BGEN file to only include the supplied rsIDs')
        bgi <- bgi[bgi$rsid %in% rsIDs, ]
      }
      snps_ids<-list(paste(bgi$chromosome, bgi$position, bgi$allele1, bgi$allele2, sep="_")) #do we want this or an external table?
      sampleIDs<-fread2(sample_file)[-1, ]$ID_2
      ind_row<-sort(match(subset_subIDs, sampleIDs))
      snp_readBGEN(bgen_file, backingfile=subset_backing_file, list_snp_id=snps_ids, ind_row = ind_row, ncores=threads, read_as ="dosage")
      
    } else {
      
      now<-Sys.time()
      message('[',now,'][Message] subsetted backing files of the BGEN file already exists')
      
    }

    rds_file<-subset_rds_file

  } else {

    if (!file.exists(bk_file)) {

      bgi<-snp_readBGI(bgi_file,snp_id=NULL)
      if (!is.null(rsIDs)) {
        message('[',now,'][Message] subsetting BGEN file to only include the supplied rsIDs')
        bgi <- bgi[bgi$rsid %in% rsIDs, ]
      }
      snps_ids<-list(paste(bgi$chromosome, bgi$position, bgi$allele1, bgi$allele2, sep="_")) #do we want this or an external table?
      snp_readBGEN(bgen_file, backingfile=backing_file, list_snp_id=snps_ids, ncores=threads, read_as ="dosage")

    } else {

      message("[", Sys.time(), "] [Message] .bk file already exists")

    }
  }
  
  obj.bigSNP <- snp_attach(rds_file)
  obj.sample <- fread2(sample_file)
  obj.sample <- obj.sample[-1, ]

  if (!is.null(subset_subIDs)) {
    
    ind_row<-sort(match(subset_subIDs, obj.sample$ID_2))
    obj.sample <- obj.sample[ind_row, ]
      
  }

  list(obj.bigSNP, obj.sample)
  
}

recover_missing_cols<-function(sumstats,map) {
  
  av_cols<-colnames(sumstats)
  #error if a1 not there ?   
  
  if ("rsid" %in% colnames(sumstats)) {
    
    sub_sumstats<-sumstats[which(sumstats$rsid  %in% map$rsid),]
    sub_map<-map[which(map$rsid  %in% sumstats$rsid),]
    sub_sumstats <- sub_sumstats[match(sub_map$rsid,sub_sumstats$rsid),]
    
    if (! "chr" %in% colnames(sumstats)) { #derive from rsid
      
      sub_sumstats$chr <- sub_map$chr
      
    }
    
    if (! "a0" %in% colnames(sumstats)) { #derive from rsid and a1
      
      a0<-rep(".", nrow(sub_map))
      
      for (i in 1:nrow(sub_map)) {
        
        a0m<-sub_map$a0[i]
        a1m<-sub_map$a1[i]
        
        if (sub_sumstats$a1[i] %in% c(a0m,a1m)) { #they appear as is, no need to reverse
          
          a0[i]<-ifelse(sub_sumstats$a1[i] == a0m, a1m, a0m)
          
        } else {
          
          a0[i]<-ifelse(get_complement(sub_sumstats$a1[i]) == a0m, get_complement(a1m), get_complement(a0m))
          
        }
        
      }
      
      sub_sumstats$a0<-a0
      
    }
    
    if (! "pos" %in% colnames(sumstats)) { #derive from rsid
      
      sub_sumstats$pos<-sub_map$pos
      
    }
    
  } else { #we do not have rsid but we have chr pos
    
    #create fake id
    map$fake_id<-paste(map$chr,map$pos,sep="_")
    sumstats$fake_id<-paste(sumstats$chr,sumstats$pos,sep="_")
    
    #exclude those that are duplicated
    sumstats<-sumstats[!(duplicated(sumstats$fake_id) | duplicated(sumstats$fake_id, fromLast=T)),]
    map<-map[!(duplicated(map$fake_id) | duplicated(map$fake_id, fromLast=T)),]
    
    #use fake_id instead of rsid
    sub_sumstats<-sumstats[which(sumstats$fake_id  %in% map$fake_id),]
    sub_map<-map[which(map$fake_id  %in% sumstats$fake_id),]
    sub_sumstats <- sub_sumstats[match(sub_map$fake_id,sub_sumstats$fake_id),]
    
    sub_sumstats$rsid<-sub_map$rsid
    
    if (! "a0" %in% colnames(sumstats)) { #derive from rsid and a1
      
      a0<-rep(".", nrow(sub_map))
      
      for (i in 1:nrow(sub_map)) {
        
        a0m<-sub_map$a0[i]
        a1m<-sub_map$a1[i]
        
        if (sub_sumstats$a1[i] %in% c(a0m,a1m)) { #they appear as is, no need to reverse
          
          a0[i]<-ifelse(sub_sumstats$a1[i] == a0m, a1m, a0m)
          
        } else {
          
          a0[i]<-ifelse(get_complement(sub_sumstats$a1[i]) == a0m, get_complement(a1m), get_complement(a0m))
          
        }
        
      }
      
      sub_sumstats$a0<-a0
      
    }
    
  }
  
  sub_sumstats
  
}

check_index<-function(x) {
  
  has_index<-FALSE
  
  if (is.null(x)) {
    
    now<-Sys.time()
    message(yellow('[',now,'][Warning] no reference .bed provided'))
    
  } else {
    
    if (! file.exists(file.path(x))) {
      
      now<-Sys.time()
      stop(red('[',now,'][Error] if provided, reference .bed must exists'))
      
    } else {
      
      has_index<-TRUE
      
    }
  }
  
  has_index
  
}

load_covPC<-function(x, nPC=NULL) {
  
  now<-Sys.time()
  
  if (endsWith(x, ".rds")) {
    
    message('[',now,'][Message] .rds file found as PC covariate (cov_PC), implying the file is an object of big_SVD (partial SVD function from bigstatsr). Left singular vectors (i.e. U) will be used as covariates.')
    obj.svd <- readRDS(x)
    covPC <- obj.svd$u
    
  } else {
    
    if (endsWith(x, ".tsv") | endsWith(x, ".csv")) {
      
      message('[',now,'][Message] .tsv/.csv file found as PC covariate (cov_PC), implying the file is a table containing N principal components, where the IDs must be the first column, followed by the principal components.')
      if (is.null(nPC)) {
        
        covPC <- fread2(x)
        
      } else {
        
        covPC <- fread2(x, select=c(1, 2:(nPC+1)))
        
      }
      colnames(covPC) <- c("IID", paste0("PC", 1:(ncol(covPC)-1)))
      rownames(covPC) = covPC$IID
      covPC$IID <- NULL
      
    } else {
      
      stop(red('[',now,'][Error] if provided, cov_PC must be .rds/.tsv/.csv'))
      
    }
  }
  
  covPC <- na.omit(covPC)
  covPC
  
}

custom_snp_match<-function(sumstats_ordered, map, join_by_pos = FALSE) {
  
  now<-Sys.time()
  
  message('[',now,'][Message] ',nrow(sumstats_ordered),' SNPs present in the sumstats')
  message('[',now,'][Message] ',nrow(map),' SNPs present in the genotypes file')
  
  sumstats_ordered$chr <- as.character(sumstats_ordered$chr)
  map$chr <- as.character(map$chr)
  sumstats_ordered$pos <- as.integer(sumstats_ordered$pos)
  map$pos <- as.integer(map$pos)
  
  if (join_by_pos) {
    matched_snps <- merge(sumstats_ordered, map, by = c("chr", "pos"), suffixes = c(".ss", ""))
  } else {
    matched_snps <- merge(sumstats_ordered, map, by = "rsid", suffixes = c(".ss", ""))
  }
  message('[',now,'][Message] ',nrow(matched_snps),' SNPs matched')
  
  matched_snps$`_NUM_ID_.ss` <- seq_len(nrow(sumstats_ordered))[match(matched_snps$rsid, sumstats_ordered$rsid)]
  matched_snps$`_NUM_ID_` <- seq_len(nrow(map))[match(matched_snps$rsid, map$rsid)]
  
  needs_flipping <- with(matched_snps, a1.ss == a0 & a0.ss == a1)
  message('[',now,'][Message] ',sum(needs_flipping),' SNPs need to be flipped')
  matched_snps$a1[needs_flipping] <- matched_snps$a0[needs_flipping]
  matched_snps$a0[needs_flipping] <- matched_snps$a1.ss[needs_flipping]
  
  matched_snps_final <- matched_snps[c("chr.ss", "rsid", "a0.ss", "a1.ss", "pos.ss", "a1freq", "beta", 
                                       "mlog10p", "n_eff", "beta_se", "p", "maf", "key", "_NUM_ID_.ss", 
                                       "pos", "_NUM_ID_")]
  colnames(matched_snps_final) <- c("chr", "rsid", "a0", "a1", "pos.ss", "a1freq", "beta", 
                                    "mlog10p", "n_eff", "beta_se", "p", "maf", "key", "_NUM_ID_.ss", 
                                    "pos", "_NUM_ID_")
  
  matched_snps_final <- matched_snps_final[order(as.integer(matched_snps_final$chr), matched_snps_final$pos), ]
  
  matched_snps_final
}