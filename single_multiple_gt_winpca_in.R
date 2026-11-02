#!/usr/bin/env Rscript
# Use a subsampled, non-redundant pon loci file!

library(rtracklayer)
library(data.table)

btn1_ids <- c("EICE18_889H", "EUCE18_1024H", "FRCE17_701H", "PACE17_421H1",
              "PACE17_433H", "PVCE17_1247H", "PVCE17_1402H")
btn2_ids <- "UGCE17_2401H"

args <- commandArgs(trailingOnly = TRUE)
#loci_file <- args[1] # non-redudant loci file containing chrom, pos, eg. from A0_A1
vaf_file  <- args[1]
out_file <- args[2]
#group     <- args[3]   # btn1 or btn2

#if (group == "btn1") {
#  ids <- btn1_ids
#} else if (group == "btn2") {
#  ids <- btn2_ids
#}
ids <-c(btn1_ids,btn2_ids)

##################
# LOAD BED FILES #
##################

bed_dir <- "/mnt/beegfs/scratch/genetics/baezortega/xx822/scubacancers_data/bed_files"

files <- list.files(bed_dir, pattern = "\\.segmentation\\.bed$",recursive = TRUE, full.names = TRUE)

names(files) <- sub("\\.segmentation\\.bed$", "", basename(files))

beds <- lapply(files[ids], import, format = "BED")

####################################
# CONVERT LOCI INTO GRANGES OBJECT #
####################################

vaf <- read.csv(vaf_file, sep = '\t', header = TRUE)
loci_gr <- GRanges(vaf$CHROM, IRanges(vaf$POS, width = 1))

snp_panel_key <- paste(vaf$CHROM, vaf$POS, sep = ":")
stopifnot(!anyDuplicated(snp_panel_key))  

##################
# WRITE CN TABLE #
##################


cn <- data.frame(chrom = vaf$CHROM, pos = vaf$POS)
for (s in ids) {
  seg <- beds[[s]]
  hit <- findOverlaps(loci_gr, seg, select = "first")
  cn[[s]] <- seg$score[hit]         
}


############################
# FILTER FOR CONSISTENT CN #
############################

btn1_cn_mat <- as.matrix(cn[btn1_ids])
btn1_consistent <- rowSums(btn1_cn_mat == btn1_cn_mat[, 1]) == ncol(btn1_cn_mat)
btn1_consistent[is.na(btn1_consistent)] <- FALSE

cn_btn1 <- btn1_cn_mat[, 1]                  
keep <- btn1_consistent & !is.na(cn_btn1) & cn_btn1 %in% 3:8
keep <- keep & !is.na(cn[[btn2_ids]])        

################################
# DETERMINE >1 COPY OR =1 COPY #
################################


# initialise the matrix with -1 for na
wpca_in <- 
  matrix(-1, nrow = nrow(vaf), ncol = length(ids) * 2,dimnames = list(snp_panel_key,
  c(paste0(ids, "_single"), paste0(ids, "_multiple"))))

for (sample in ids) {
  v <- vaf[[sample]]
  thres <- 1 /cn[[sample]]                  # per-sample CN, row-aligned with vaf
  ok  <- keep & !is.na(thres) & !is.na(v)
  
  single   <- ok & v > 0& v <= thres    # 1 copy
  mult_het <- ok & v > thres & v <= 0.9  # >1 copy, heterozygous
  mult_hom <- ok & v > thres & v >  0.9      # >1 copy, homozygous
  
  wpca_in[single,   paste0(sample, "_single")]   <- 1
  wpca_in[mult_het, paste0(sample, "_multiple")] <- 1
  wpca_in[mult_hom, paste0(sample, "_multiple")] <- 2
}

out <- data.table(chrom = vaf$CHROM, pos = vaf$POS, as.data.table(wpca_in))
fwrite(out, out_file, sep = "\t")
