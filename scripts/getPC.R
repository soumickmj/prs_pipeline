#!/usr/bin/env Rscript

#script to perform PCA, using Truncated SVD while limiting LD (cf. https://privefl.github.io/bigsnpr/reference/snp_autoSVD.html)

source("/home/soumick.chatterjee/Codes/GitHub/prs_pipeline/scripts/utils.R")

sessionInfo()

options(warn = -1)

option_list = list(
  make_option(c('-i', '--input'), action='store', type='character', help='.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]', default='/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/cond_plus_plink_maf1p_geno10p_caucasian_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4.bgen'),
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]', default="/group/glastonbury/soumick/PRS/inputs/F20208v3_DiffAE_select_latents_r80_discov_INF30/PC/prova_cond_plus_plink_maf1p_geno10p_caucasian_prune_250_5_r0p5_ukbb_autosomes_mac100_info0p4"),
  make_option(c('--threads'), action='store', type='numeric', help='computing threads [1]', default=32)
)

opt = parse_args(OptionParser(option_list=option_list))

now<-Sys.time()
message('[',now,'][Message] performing pre-flight checks') 
is_bgen<-check_input(opt$input)
check_output(opt)

#Read and prepare genotype data
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
message('[',now,'][Message] genotypes loaded and processed')

# obj.svd <- snp_autoSVD(G, 
#                        CHR,
#                        POS, 
#                        ncores = opt$threads)

obj.svd <- snp_autoSVD(G, 
                       CHR,
                       POS, 
                       ncores = opt$threads,
                       k=1,
                       size=10,
                       max.iter=1,
                       thr.r2=NA,
                       min.mac=100)

now<-Sys.time()
message('[',now,'][Message] PCA finished!')

pdf(file.path(paste0(opt$output,'.truncSVD.PC.pdf')))
plot(obj.svd)
plot(obj.svd, type = "scores") 
dev.off()

now<-Sys.time()
message('[',now,'][Message] Saving the PCs..')


PC_file<-file.path(paste0(opt$output,'.geno.truncSVD.PC.rds'))
saveRDS(obj.svd, PC_file)