rm(list=ls(all=TRUE))

source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")

## create bed and bed12 file for features
if(0){

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

  outdir <- "C:/GREENBLATT/genomic_feature"
  setwd(outdir)
  feature_source <- "gencode.v19.annotation.gtf"
  #
  for(featureName in c("transcript", "utr5", "cds", "utr3", "exon", "intron")){
    gtf_to_bed_longest_tx(txdb, featureName, feature_source, export=TRUE)
  }


  tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
  tl <- tl[tl$cds_len > 0, ]
  longest_tx <- extract_longest_tx(txdb, plot=TRUE)
  longest_tx <- longest_tx[order(longest_tx$cds_len),]
  tl <- tl[tl$tx_name %in% longest_tx$tx_name, ]
  tl <- tl[order(tl$cds_len),]

  means <- c(promoter=1000, apply(tl[,c(7,6,8)], 2, mean), TTS=1000)
  scaled_bins <- round(means*1000/sum(means))

}


if(0){
gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
txdb <-  makeTxDbFromGFF(gtffile)
ext <- 400
wd <- "C:/GREENBLATT/Nabeel/kristie/clip_reads_around_junction"
setwd(wd)
#queryfiles <- c("combined_crosslink_YTHDF2.uniq.bed") #"MR_Ac4c_merged.thUni.bed",
queryfiles <- c("TR_Ac4c_merged.thUni.bam", "MR_Ac4c_merged.thUni.bam", "NAT10_merged.thUni.bam")
querylabels <- c("TR", "MR", "NAT10")
pdf("TR_MR_NAT10_merged.bam_metagene.pdf")
data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, sn=0, norm=TRUE, CLIP_reads=F)
dev.off()

featureName <- "cds"
hl <- c(-100, 0, 0, 100) ## boundaries for highlight, 1,2 for start, 3,4 for end
sn <- 100000
ran <- FALSE

pdf("TR_Ac4c_merged.thUni_intron.pdf")
plot_start_end_gtf(queryfiles, querylabels, txdb, featureName, ext, hl, sn, ran)
dev.off()
}

## down-sampling noChange ZNF281 ####

if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")
  qapa_dir <- "C:/GREENBLATT/Nabeel/Nujhat/qapa/differential_analysis_ZNF281_KD10_vs_siNT_0_100"
  downGene <- readLines(file.path(qapa_dir, "PAU_analysis_siZNF281_down_list.txt"))
  upGene <- readLines(file.path(qapa_dir, "PAU_analysis_siZNF281_up_list.txt"))
  noChangeGene <- readLines(file.path(qapa_dir, "PAU_analysis_siZNF281_noChange_list.txt"))

  rnaseq_dir <- "C:/GREENBLATT/Nabeel/Nujhat/RNAseq"
  rna_s1 <- read.delim(file.path(rnaseq_dir, "Greenblatt_003_Strip1_C01_siZNF281_S1_S3.genes.results"), sep="\t", header=T)
  rna_s2 <- read.delim(file.path(rnaseq_dir, "Greenblatt_004_Strip1_D01_siZNF281_S2_S4.genes.results"), sep="\t", header=T)

  tpm1 <- rna_s1[, c("gene_id", "TPM")]
  tpm2 <- rna_s2[, c("gene_id", "TPM")]
  tpm <- merge(tpm1, tpm2, by=("gene_id"))
  mean_tpm <- apply(tpm[, 2:3], 1, mean)
  names(mean_tpm) <- unlist(lapply(tpm$gene_id, function(x) unlist(strsplit(x, split="_"))[2]))

  up_down_rna <- mean_tpm[names(mean_tpm) %in% c(downGene, upGene)]
  noChange_rna <- mean_tpm[names(mean_tpm) %in% noChangeGene]
  up_down_rna[is.na(up_down_rna)] <- 0
  noChange_rna[is.na(noChange_rna)] <- 0

  sub_rna <- sub_sample(noChange_rna, up_down_rna, 0.95)
  for(utr in c("putr", "dutr")){
    #utr <- "putr"
    noChange <- read.delim(file.path(qapa_dir, paste("PAU_analysis_siZNF281_noChange_list_",utr,".bed", sep="")), header=F)
    sub_noChange <- noChange[noChange[,4] %in% names(sub_rna),]
    write.table(sub_noChange, file.path(qapa_dir, paste("PAU_analysis_siZNF281_noChange_list_",utr,"_sub.bed", sep="")), sep="\t", col.names=F, row.names=F, quote=F)
  }


}


## ZNF281 clip reads around cleavage sites ####
if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

  chainFile <- "C:/GREENBLATT/genomic_feature/hg38ToHg19.over.chain"

  ext <- 200
  wd <- "C:/GREENBLATT/Nabeel/Nujhat/profile_at_cleavageD"
  setwd(wd)

  qapa_dir <- "C:/GREENBLATT/Nabeel/Nujhat/qapa/differential_analysis_ZNF281_KD20_vs_siNT_0_100"


### all 3'_UTR
  centerfiles <- file.path(qapa_dir, "UTR3_coordinates_liftOver.bed")
  centerlabels <- c("all_3utr")

  queryfiles <- c("ZNF281.merged.bam", "NUDT21.merged.bam", "CPSF7.merged.bam", "CSTF2_GSM917676_275_+_density.bw", "CSTF2T_merged_+_density.bw", "Im68_merged_+_density.bw" )
  querylabels <- c("ZNF281", "NUDT21", "CPSF7", "CSTF2", "CSTF2T", "CPSF6")

  if(1){
    pdf(paste("ZNF281_CPSF_CSTF_clip_read_starts_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=T, stats=F, scale=T)
    dev.off()
  }
  if(1){
    pdf(paste("ZNF281_CPSF_CSTF_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=F, stats=F, scale=T)
    dev.off()
  }

  queryfiles <- c("ZNF281.merged.bam", "NUDT21.merged.bam", "CPSF7.merged.bam", "CPSF160.merged.bam", "CPSF2.merged.bam", "FIP1.merged.bam")
  querylabels <- c("ZNF281", "NUDT21", "CPSF7", "CPSF160", "CPSF2", "FIP1")

  if(1){
    pdf(paste("ZNF281_CPSF_clip_read_starts_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=T, stats=F, scale=T)
    dev.off()
  }
  if(1){
    pdf(paste("ZNF281_CPSF_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=F, stats=F, scale=T)
    dev.off()
  }

  queryfiles <- c("ZNF281.ABC_crosslinkregions_wInputwCL.bed", "NUDT21_CDE_crosslinkregions_wInputwCL.bed", "CPSF7_AB_crosslinkregions_wInputwCL.bed", "CPSF160_AB_crosslinkregions_wInputwCL.bed", "CPSF2_AB_crosslinkregions_wInputwCL.bed", "FIP1_AB_crosslinkregions_wInputwCL.bed")
  querylabels <- c("ZNF281", "NUDT21", "CPSF7", "CPSF160", "CPSF2", "FIP1")

  if(1){
    pdf(paste("ZNF281_CPSF_Peaks_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", smo=T, norm=F, stats=F, scale=T)
    dev.off()
  }

### PAU regulated cleavage sites
  for(utr in c("putr", "dutr")){
    centerfiles_hg38 <- c(file.path(qapa_dir, paste("PAU_analysis_siZNF281_up_list_",utr,".bed", sep="")),
                     file.path(qapa_dir, paste("PAU_analysis_siZNF281_down_list_",utr,".bed", sep="")),
                      file.path(qapa_dir, paste("PAU_analysis_siZNF281_noChange_list_",utr,"_sub.bed", sep="")))
    centerfiles <- NULL
    for(queryfile in centerfiles_hg38){
      newFile <- liftOverBed(chainFile, queryfile)
      centerfiles <- c(centerfiles, newFile)
    }
    centerlabels <- c("Up", "Down", "noChange")

    hl <- c(-50, 50) ## c(-50, 50), c(-100, 0)
    hln <- paste(hl, collapse="-")
    ext <- 200
    # reads


    if(1){
      queryfiles <- c("ZNF281.merged.bam")
      querylabels <- c("ZNF281")
      Inputfile <- "Input.SP1.A.thUni.bam"

      op <- paste("ZNF281_merged_clip_ratioOverInput_around_cleavage_sites", hln, utr, sep="_")
      plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=Inputfile, ratioOverInput=T, outPrefix=op)

    }


    if(1){
      queryfiles <- c("NUDT21.merged.bam")
      querylabels <- c("NUDT21")
      Inputfile <- "Input.NUDT21.A.thUni.bam"

      op <- paste("NUDT21_merged_clip_ratioOverInput_around_cleavage_sites", hln, utr, sep="_")
      plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=Inputfile, ratioOverInput=T, outPrefix=op)

    }


    if(1){
      queryfiles <- c("CPSF7.merged.bam")
      querylabels <- c("CPSF7")
      Inputfile <- "Input.CPA.A.thUni.bam"

      op <- paste("CPSF7_merged_clip_ratioOverInput_around_cleavage_sites", hln, utr, sep="_")
      plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=Inputfile, ratioOverInput=T, outPrefix=op)

    }

    if(1){
      queryfiles <- c("CSTF2_GSM917676_275_+_density.bw", "CSTF2T_merged_+_density.bw", "Im68_merged_+_density.bw")
      names(queryfiles) <- c("CSTF2" ,"CSTF2T", "CPSF6")
      for(querylabels in names(queryfiles)){
        queryfiles <- queryfiles[querylabels]

        op <- paste(querylabels,"_clip_reads_around_cleavage_sites", hln, utr, sep="_")
        plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=NULL, ratioOverInput=F, outPrefix=op)

      }
    }
  }

  ### reads over peaks
  centerfiles <- c("ZNF281.ABC_crosslinkregions_wInputwCL.bed")
  centerlabels <- c("ZNF281_pureclip_peaks")


  queryfiles <- c("ZNF281.merged.bam")
  querylabels <- c("ZNF281ABC")

  if(1){
    pdf(paste("ZNF281_merged_clip_reads_around_pureclip_peaks.pdf", sep=""), height=8, width=12)
    plot_reference_locus_with_random(txdb=txdb, queryfiles=queryfiles, centerfiles=centerfiles, ext=ext, hl=hl, querylabels=querylabels, centerlabels=centerlabels, CLIP_reads = F)
    dev.off()
  }


  ### non-regulated putr dutr sutr

  centerfiles_hg38 <- c(file.path(qapa_dir, "PAU_analysis_putr.bed"),
                        file.path(qapa_dir, "PAU_analysis_dutr.bed"),
                        file.path(qapa_dir, "PAU_analysis_sutr.bed"))
                        #file.path(qapa_dir, "UTR3_coordinates.bed"))

  centerfiles <- NULL
  for(queryfile in centerfiles_hg38){
    newFile <- liftOverBed(chainFile, queryfile)
    centerfiles <- c(centerfiles, newFile)
  }
  centerlabels <- c("putr", "dutr", "sutr")

  hl <- c(-50, 50) ## c(-50, 50), c(-100, 0)

  queryfiles <- c("CPSF7.merged.bam")
  querylabels <- c("CPSF7AB")
  if(1){
    pdf(paste("CPSF7_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  queryfiles <- c("Im68_GSM917665_273_+_density.wig", "Im68B_GSM917664_274_+_density.wig")
  querylabels <- c("Im68A", "Im68B")
  if(1){
    pdf(paste("Im68_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }


  # reads
  queryfiles <- c("ZNF281.A.thUni.bam", "ZNF281.B.thUni.bam", "ZNF281.C.thUni.bam")
  querylabels <- c("ZNF281A", "ZNF281B", "ZNF281C")
  Inputfile <- "Input.SP1.A.thUni.bam"

  if(0){
    pdf(paste("ZNF281_individual_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }
  if(0){
    pdf(paste("ZNF281_individual_clip_ratioOverInput_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", inputfile=Inputfile, ratioOverInput = T)
    dev.off()
  }

  queryfiles <- c("ZNF281.AC.merged.bam", "ZNF281.merged.bam")
  querylabels <- c("ZNF281AC", "ZNF281ABC")

  if(1){
    pdf(paste("ZNF281_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  if(1){
    pdf(paste("ZNF281_merged_clip_ratioOverInput_around_cleavage_site.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site", inputfile=Inputfile, ratioOverInput = T)
    dev.off()
  }

  queryfiles <- c("NUDT21.merged.bam")
  querylabels <- c("NUDT21CDE")

  if(1){
    pdf(paste("NTDT21_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }


  queryfiles <- c("CSTF2_GSM917676_275_+_density.bw", "CSTF2T_A_GSM917677_283_+_density.bw")
  querylabels <- c("CSTF2", "CSTF2T")

  if(1){
    pdf(paste("CSTF2_CSTF2T_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  queryfiles <- c("CSTF2T_merged_+_density.bw")
  querylabels <- c("CSTF2TAB")

  if(1){
    pdf(paste("CSTF2T_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  queryfiles <- c("Im68_merged_+_density.bw")
  querylabels <- c("Im68AB")

  if(1){
    pdf(paste("Im68_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  if(1){
    pdf(paste("metagene_profile_of_sutrs.pdf", sep=""), height=8, width=12)
    data_df <- plot_5parts_metagene(centerfiles[3], centerlabels[3], txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=F)
    dev.off()
  }

}


## ZNF281 clip reads metagene plot ####
if(1){
    gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
    txdb <-  makeTxDbFromGFF(gtffile)
    wd <- "C:/GREENBLATT/Nabeel/Nujhat/profile_at_cleavage"
    setwd(wd)

    queryfiles <- c("ZNF281.merged.bam")#, "Input.SP1.A.thUni.bam")
    querylabels <- c("ZNF281ABC")#, "Input")

    pdf("metagene_profile_of_ZNF281_merged_bam_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, CLIP_reads=F, smo=T)
    dev.off()

    queryfiles <- c("NUDT21.merged.bam")#, "Input.NUDT21.A.thUni.bam")
    querylabels <- c("NUDT21CDE")#, "Input")

    pdf("metagene_profile_of_NUDT21_merged_bam_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, CLIP_reads=F, smo=T)
    dev.off()

    queryfiles <- c("CPSF7.merged.bam")#, "Input.NUDT21.A.thUni.bam")
    querylabels <- c("CPSF7AB")#, "Input")

    pdf("metagene_profile_of_CPSF7_merged_bam_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, CLIP_reads=F, smo=T)
    dev.off()


    pdf("metagene_profile_of_ZNF281_merged_bam_500bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=500, CLIP_reads=F)
    dev.off()

    queryfiles <- c("ZNF281_over_Input_DESeq2_insignificant_crossLink_site.bed", "ZNF281_over_Input_DESeq2_significant_crossLink_site.bed")
    querylabels <- c("insignificant_peaks", "significant_peaks")

    pdf("metagene_profile_of_ZNF281_significant_peaks_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F)
    dev.off()

    pdf("metagene_profile_of_ZNF281_significant_peaks_500bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=500, CLIP_reads=F)
    dev.off()

    queryfiles <- c("ZNF281.ABC_crosslinkregions_wInputwCL.bed")
    querylabels <- c("pureclip_peaks")

    pdf("metagene_profile_of_ZNF281_pureclip_peaks_100bin_scaled.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T, scale=T)
    dev.off()

    pdf("metagene_profile_of_ZNF281_pureclip_peaks_500bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=500, CLIP_reads=F, smo=T)
    dev.off()

    queryfiles <- c("CSTF2_GSM917676_275_+_density.bw")
    querylabels <- c("CSTF2")

    if(1){
      pdf(paste("metagene_profile_of_CSTF2_downloaded_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    queryfiles <- c("Im68_GSM917665_273_+_density.bw", "Im68B_GSM917664_274_+_density.bw")
    querylabels <- c("Im68A", "Im68B")

    if(1){
      pdf(paste("metagene_profile_of_CFIm68_downloaded_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    if(1){
      pdf(paste("metagene_profile_of_Input_bam_100bin.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(Inputfile, "Input", txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    queryfiles <- c("CSTF2T_merged_+_density.bw")
    querylabels <- c("CSTF2TAB")

    if(1){
      pdf(paste("metagene_profile_of_CSTF2T_merged_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    queryfiles <- c("Im68_merged_+_density.bw")
    querylabels <- c("Im68AB")

    if(1){
      pdf(paste("metagene_profile_of_CFIm68_merged_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }
}

## down-sampling noChange ZBTB7A ####

if(1){
  qapa_dir <- "C:/GREENBLATT/Nabeel/RNAseq/PAU_results/differential_analysis_C2H2_KD_vs_siNT_0_100"
  downGene <- readLines(file.path(qapa_dir, "PAU_analysis_ZBTB7A_down_list.txt"))
  upGene <- readLines(file.path(qapa_dir, "PAU_analysis_ZBTB7A_up_list.txt"))
  noChangeGene <- readLines(file.path(qapa_dir, "PAU_analysis_ZBTB7A_noChange_list.txt"))

  rnaseq_dir <- "C:/GREENBLATT/Nabeel/RNAseq"
  rna_s1 <- read.delim(file.path(rnaseq_dir, "Greenblatt_009_Strip_2_A02_ZBTB7A_S1_R_S9.genes.results"), sep="\t", header=T)
  rna_s2 <- read.delim(file.path(rnaseq_dir, "Greenblatt_010_Strip_2_B02_ZBTB7A_S2_R_S10.genes.results"), sep="\t", header=T)

  tpm1 <- rna_s1[, c("gene_id", "TPM")]
  tpm2 <- rna_s2[, c("gene_id", "TPM")]
  tpm <- merge(tpm1, tpm2, by=("gene_id"))
  mean_tpm <- apply(tpm[, 2:3], 1, mean)
  names(mean_tpm) <- unlist(lapply(tpm$gene_id, function(x) unlist(strsplit(x, split="_"))[2]))

  up_down_rna <- mean_tpm[names(mean_tpm) %in% c(downGene, upGene)]
  noChange_rna <- mean_tpm[names(mean_tpm) %in% noChangeGene]
  up_down_rna[is.na(up_down_rna)] <- 0
  noChange_rna[is.na(noChange_rna)] <- 0

  sub_rna <- sub_sample(noChange_rna, up_down_rna, 0.90)
  for(utr in c("putr", "dutr")){
    #utr <- "putr"
    noChange <- read.delim(file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_noChange_list_",utr,".bed", sep="")), header=F)
    sub_noChange <- noChange[noChange[,4] %in% names(sub_rna),]
    write.table(sub_noChange, file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_noChange_list_",utr,"_sub.bed", sep="")), sep="\t", col.names=F, row.names=F, quote=F)
  }


}



## ZBTB7A clip reads around cleavage sites ####
if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

  chainFile <- "C:/GREENBLATT/genomic_feature/hg38ToHg19.over.chain"

  ext <- 200
  wd <- "C:/GREENBLATT/Nabeel/ZBTB7A"
  setwd(wd)

  qapa_dir <- "C:/GREENBLATT/Nabeel/ZBTB7A"

  queryfile <- file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_up_list_","putr",".bed", sep=""))
  for(utr in c("putr", "dutr")){
    centerfiles_hg38 <- c(file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_up_list_",utr,".bed", sep="")),
                          file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_down_list_",utr,".bed", sep="")),
                          file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_noChange_list_",utr,"_sub.bed", sep="")))
    centerfiles <- NULL
    for(queryfile in centerfiles_hg38){
      newFile <- liftOverBed(chainFile, queryfile)
      centerfiles <- c(centerfiles, newFile)
    }
    centerlabels <- c("Up", "Down", "noChange")

    hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

    # reads
    queryfiles <- c("ZBTB7A.A.thUni.bam", "ZBTB7A.B.thUni.bam", "ZBTB7A.ABCD.merged.bam")
    querylabels <- c("ZBTB7A.A", "ZBTB7A.B", "ZBTB7A.ABCD")
    pdf(paste("ZBTB7A_clip_reads_around_cleavage_sites3_",utr,".pdf", sep=""), height=8, width=12)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()

  }

  ### plot intron retention
  ext <- 500
  hls <- list(start=c(0, 500), end=c(-500,0))
  wd <- "C:/GREENBLATT/Nabeel/ZBTB7A"
  setwd(wd)

  vast_dir <- "C:/GREENBLATT/Nabeel/VAST-TOOLS/ZBTB7A"
  centerfiles <- extract_intron_VAST("C:/GREENBLATT/Nabeel/VAST-TOOLS/ZBTB7A/Intron_retention_regulated_events.tab", length_filter = 500L)
  centerlabels <- c("IR_up", "IR_down", "IR_noChange")

  queryfiles <- c("ZBTB7A.A.thUni.bam", "ZBTB7A.B.thUni.bam", "ZBTB7A.ABCD.merged.bam")
  querylabels <- c("ZBTB7A.A", "ZBTB7A.B", "ZBTB7A.ABCD")
  for(refp in c("start", "end")){
    hl <- hls[[refp]]
    op <- paste0("ZBTB7A_clip_reads_around_regulated_intron_", refp)
    plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint=refp, Xlab=refp, stat=T, scale=F, smo=T, norm=F, outPrefix=op)
  }
}


## NAT10 reads around AC4C peaks ####
if(0){

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 400
  wd <- "C:/GREENBLATT/Nabeel/kristie/NAT10/reads_metagene_profile"
  setwd(wd)

  centerfiles <- c("TR_merged_N_at_CIMS.bed", "MR_merged_N_at_CIMS.bed")
  centerlabels <- c("TR_AC4C", "MR_AC4C")

  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  # reads
  queryfiles <- c("NAT10_merged.thUni.bam")
  querylabels <- c("NAT10")

  pdf("reads_distance_to_center_AC4C_with_random_narrow.pdf", height=8, width=12)
  plot_reference_locus_with_random(txdb, queryfiles, centerfiles, ext, hl, querylabels, centerlabels, binsize=10, refPoint="center", CLIP_reads=FALSE, stranded=TRUE)
  dev.off()

  ## peaks
  queryfiles <- c("NAT10_AB_crosslinkregions_wInputwCL.bed")
  querylabels <- c("NAT10")
  pdf("peaks_distance_to_center_AC4C.pdf", height=8, width=12)
  plot_reference_locus(queryfiles, centerfiles, ext, hl, sn, querylabels, centerlabels)
  dev.off()

  pdf("peaks_distance_to_center_AC4C_with_random.pdf", height=8, width=12)
  plot_reference_locus_with_random(txdb, queryfiles, centerfiles, ext, hl, querylabels, centerlabels)
  dev.off()
}

if(1){
  refdir <- "C:/GREENBLATT/genomic_feature"
  #gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  #txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 100
  wd <- "C:/GREENBLATT/Nabeel/kristie/clip_reads_around_junction"
  setwd(wd)
  centerfiles <- file.path(refdir, c("gencode.v19.protein_coding_transcript_EXON.bed", "gencode.v19.protein_coding_transcript_INTRON.bed"))
  centerlabels <- c("Exon", "Intron")
  queryfiles <- c("TR_Ac4c_merged.thUni.bed", "MR_Ac4c_merged.thUni.bed")
  querylabels <- c("TR_AC4C", "MR_AC4C")

  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  pdf("distance_to_ExonCenter_100_AC4C.pdf", height=8, width=12)
  plot_reference_locus(queryfiles, centerfiles, ext, hl, sn, querylabels, centerlabels)
  dev.off()
}

## Nabeel m6A ####
if(1){

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 400
  wd <- "C:/GREENBLATT/Nabeel/m6A/cits_peak"
  setwd(wd)
  queryfiles <- c("m6A_FTO_merged.thUni.crossLink_site.0.05.bed", "m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed")
  centerlabels <- c("Hek")
  centerfiles <- c("m6A_Hek_merged.thUni.crossLink_site.0.05.bed")
  querylabels <- c("FTO", "ZBTB48")
  sn <- 100000
  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  pdf("distance_to_Hek_peaks.pdf", height=8, width=12)
  plot_reference_locus(queryfiles, centerfiles, ext, hl, sn, querylabels, centerlabels)
  dev.off()

  ## plot for gene feature boundaries
  queryfiles <- c("m6A_FTO_merged.thUni.crossLink_site.0.05.bed",
                "m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed",
                "m6A_Hek_merged.thUni.crossLink_site.0.05.bed")
  querylabels <- c("FTO", "ZBTB48", "Hek")

  queryfiles <- c("m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed")
  querylabels <- c("ZBTB48")

  pdf("metagene_profile_of_peaks_longest_transcript.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, sn=0, norm=FALSE, longest=TRUE)
  dev.off()

  queryfiles <- c("m6A_FTO_merged.thUni.bed")
  querylabels <- c("FTO_reads")
  pdf("metagene_profile_of_merged_FTO_bed.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, sn=0, norm=F)
  dev.off()

  for(featureName in c("utr5", "utr3", "cds")){

    ran <- ifelse(featureName %in% c("exon", "intron"), TRUE, FALSE)
    ## hl is boundaries for highlight, 1,2 for start, 3,4 for end
    hl <- c(0, 100, -100, 0)
    if(featureName %in% c("intron")) hl <- c(-100, 0, 0, 100)

    pdf(paste("m6A_FTO_reads_at_", gsub("'", "",featureName), ".pdf", sep=""), height=8, width=12)
    plot_start_end_gtf(queryfiles, querylabels, txdb, featureName, ext, hl, sn, ran)
    dev.off()
  }
}


## plot metagene profile for NAT10 reads #############
if(1){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/kristie/NAT10/reads_metagene_profile"
  setwd(wd)
  queryfiles <- "NAT10_merged.thUni.bam"
  querylabels <- "NAT10_merged"

  pdf("metagene_profile_of_NAT10_merged_bam_100bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, sn=0, CLIP_reads=F)
  dev.off()

  pdf("metagene_profile_of_NAT10_merged_bam_500bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=500, sn=0, CLIP_reads=F)
  dev.off()
}

## plot metagene profile for NAT10 peaks #############
if(1){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/kristie/NAT10/peak_metagene_profile"
  setwd(wd)
  queryfiles <- "NAT10_AB_crosslinkregions_wInputwCL.bed"
  querylabels <- "NAT10_merged_peaks"

  pdf("metagene_profile_of_NAT10_merged_peaks_100bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, sn=0, CLIP_reads=F)
  dev.off()

  pdf("metagene_profile_of_NAT10_merged_peaks_500bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=500, sn=0, CLIP_reads=F)
  dev.off()
}



## Kristie AC4C ####

if(0){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 400
  wd <- "C:/GREENBLATT/Nabeel/kristie/clip_reads_around_junction"
  setwd(wd)
  #queryfiles <- c("combined_crosslink_YTHDF2.uniq.bed") #"MR_Ac4c_merged.thUni.bed",
  queryfiles <- c("TR_Ac4c_merged.thUni.bed", "MR_Ac4c_merged.thUni.bed")
  querylabels <- c("TR_AC4C", "MR_AC4C")

  featureName <- "5'utr"
  hl <- c(-100, 0, 0, 100) ## boundaries for highlight, 1,2 for start, 3,4 for end
  sn <- 100000
  ran <- TRUE

  pdf("TR_Ac4c_merged.thUni_intron.pdf")
  plot_start_end_gtf(queryfiles, txdb, featureName, ext, hl, sn, ran)
  dev.off()
}

if(0){

  #gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  #txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 400
  wd <- "C:/GREENBLATT/Nabeel/kristie/clip_reads_around_junction"
  setwd(wd)
  centerfiles <- c("NAT10_merged.thUni.bed", "combined_crosslink_YTHDF2.uniq.bed")
  centerlabels <- c("NAT10", "YTHDF2_peaks")
  queryfiles <- c("TR_Ac4c_merged.thUni.bed", "MR_Ac4c_merged.thUni.bed")
  querylabels <- c("TR_AC4C", "MR_AC4C")
  sn <- 10000
  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  pdf("distance_to_center_AC4C.pdf")
  plot_reference_locus(queryfiles, centerfiles, ext, hl, sn, querylabels, centerlabels)
  dev.off()
}


## Jingwen ChIPseq reads #############
if(1){
  options(error=recover)
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Jingwen/ChIPseq/bam_nabeel"
  qapa <- "C:/GREENBLATT/Jingwen/GO_results_25_75"
  setwd(wd)

  queryfiles <- c("SP1_ChIP-12_merged.bam") ## "SP1_SRR30831245_merged.bam",
  querylabels <- c("SP1_ChIP_Nabeel") ## "SP1_ChIP_SRR3083124",
  inputfiles <- c("SP1_Input-12_merged.bam") ## "SP1_Input-12_merged1.bam",
  inputlabels <- c("SP1_Input_Nabeel") ## "SP1_Input_Nabeel1",

  qTF <- "siSP1"
  xlabs <- c("Distal Cleavage Site", "Proximal Cleavage Site")
  names(xlabs) <- c("dutr", "putr")
  for(utr in c("dutr", "putr")){
    ext <- c(-150, 150)
    hl <- c(-150,0)
    xlab <- xlabs[utr]
    centerfiles <- c(file.path(qapa, paste("PAU_analysis_",qTF, "_down_list_",utr,".bed", sep="")),
                     file.path(qapa, paste("PAU_analysis_",qTF, "_up_list_",utr,".bed", sep=""))
                     #file.path(qapa, paste("PAU_analysis_",qTF, "_noChange_list_",utr,".bed", sep=""))
    )
    #file.path(qapa_dir, "UTR3_coordinates.bed"))

    centerlabels <- c("Lenthening genes", "Shortening genes")

    outPrefix <- paste("SP1_ChIPseq_reads_around_",utr,"_cleavage_sites", sep="")
    plot_reference_locus(queryfiles=queryfiles, centerfiles=centerfiles, ext=ext, hl=hl, querylabels=querylabels, centerlabels=centerlabels,
                         inputfiles=inputfiles, inputlabels=inputlabels, refPoint="end", stats=T, norm=T, smo=T, Xlab=xlab, outPrefix=outPrefix)
  }



  TFs <- c("SP1", "CTCF", "MAZ", "OSR2", "ZBTB48", "ZNF281")

  cl <- start_parallel(length(TFs))

  parLapply(cl, TFs, function(TF){

    source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")
    gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
    txdb <-  makeTxDbFromGFF(gtffile)
    wd <- "C:/GREENBLATT/Jingwen/ChIPseq"
    qapa <- "C:/GREENBLATT/Jingwen/GO_results_25_75"
    setwd(wd)

    hl_list <- list(c(-500, -250), c(-250,0), c(0,250), c(250,500))## c(-50, 50), c(-100, 0)

    if(0){
      queryfiles <- paste(TF, ".sorted.bam", sep="")
      pdf(paste("metagene_profile_of_",TF, "_ChIPseq_100bin.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(queryfiles, TF, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F)
      dev.off()

      if(TF %in% c("MAZ")){
        queryfiles <- paste(TF, ".merged_clip.bam", sep="")
        pdf(paste("metagene_profile_of_",TF, "_CLIP_100bin.pdf", sep=""), height=8, width=12)
        data_df <- plot_5parts_metagene(queryfiles, TF, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F)
        dev.off()
      }
    }

    for(i in 1:length(hl_list)){
      hl <- hl_list[[i]]


      qTF <- ifelse(TF %in% c("SP1", "MAZ"), paste("si",TF,sep=""), TF)
      for(utr in c("dutr", "putr")){
        ext <- 500
        centerfiles <- c(file.path(qapa, paste("PAU_analysis_",qTF, "_down_list_",utr,".bed", sep="")),
                         file.path(qapa, paste("PAU_analysis_",qTF, "_up_list_",utr,".bed", sep=""))
                         #file.path(qapa, paste("PAU_analysis_",qTF, "_noChange_list_",utr,".bed", sep=""))
                         )
        #file.path(qapa_dir, "UTR3_coordinates.bed"))

        centerlabels <- paste(c("PAU-down", "PAU-up"), utr, sep="_")


        if(1){
          queryfiles <- paste(TF, ".sorted.bam", sep="")
          querylabels <- paste(TF, "ChIPseq")
          outPrefix <- paste(TF,i, "_ChIPseq_reads_around_",utr,"_cleavage_sites", sep="")
          plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", stats=T, norm=F, smo=F, Xlab="Cleavage site",outPrefix=outPrefix)



          if(TF %in% c("MAZ")){
            queryfiles <- paste(TF, ".merged_clip.bam", sep="")
            querylabels <- paste(TF, "CLIP")
            outPrefix <- paste(TF,i, "_CLIP_reads_around_",utr,"_cleavage_sites", sep="")
            plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="end", stats=T, norm=F, smo=F, Xlab="Cleavage site",outPrefix=outPrefix)


          }
        }


  ## find the transcripts that correspond to the utr3
        if(utr == "putr"){
          centerfiles <- unlist(lapply(centerfiles, function(x)utr3_to_transcript(x, txdb)))
          ext <- 1000

          if(1){
            queryfiles <- paste(TF, ".sorted.bam", sep="")
            querylabels <- paste(TF, "ChIPseq")
            outPrefix <- paste(TF,i, "_ChIPseq_reads_around_TSS", sep="")
            plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="start", stats=T, norm=F, smo=F, Xlab="TSS",outPrefix=outPrefix)



            if(TF %in% c("MAZ")){
              queryfiles <- paste(TF, ".merged_clip.bam", sep="")
              querylabels <- paste(TF, "CLIP")
              outPrefix <- paste(TF,i, "_CLIP_reads_around_TSS", sep="")
              plot_reference_locus(queryfiles, centerfiles, ext, hl, querylabels, centerlabels, refPoint="start", stats=T, norm=F, smo=F, Xlab="TSS",outPrefix=outPrefix)

            }
          }
        }
      }
    }
  })
  stop_parallel(cl)
}

### Nabeel TDP43 #######################
if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

### R_loop
  wd_r <- "C:/GREENBLATT/Nabeel/TDP43/R_loop"
  wd_chip <- "C:/GREENBLATT/Nabeel/TDP43/TDP43_ChIP"
  wd_clip <- "C:/GREENBLATT/Nabeel/TDP43/CLIP"
  wd_cor <- "C:/GREENBLATT/Nabeel/TDP43/Correlations"
  setwd(wd_r)

  queryfiles <- list.files(wd, pattern=".p.bw")
  querylabels <- paste(rep(c("V5", "Input"), 3), rep(c("rep1", "rep2", "rep3"),each=2), sep="_")

  op <- "metagene_profile_of_R_loop"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op)


  queryfiles <- c("GSM2550993_chipseq.HEK293.D210N_V5.p.bw", "GSM2551003_chipseq.HEK293.D210N_input.p.bw")
  querylabels <- c("V5", "Input")

  op <- "metagene_profile_of_R_loop_merged"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op)

  centerfiles <- c("level9_R-loop_zones.bed", "level5_R-loop_zones.bed")
  centerlabels <- c("R_loop_level9", "R_loop_level5")

  op  <- "R_loop_signal_around_level9level5_peaks"
  system.time(
    plot_reference_locus(queryfiles, centerfiles, 500, c(-50,50), querylabels, centerlabels, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op, genome="hg19")
  )


### ChIPseq

  setwd(wd_chip)

  queryfiles <- list.files(wd, pattern=".bed")
  querylabels <- c("all_Peak", "conservative_peak", "optimal_peak")

  op <- "metagene_profile_of_TDP43_ChIPseq_peak_Norm"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op)


  queryfiles <- list.files(wd, pattern=".bigWig")
  querylabels <- c("pValue", "FCOC")

  op <- "metagene_profile_of_TDP43_ChIPseq_signal"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)


  ## CLIP

  setwd(wd_clip)
  queryfiles <- "TDP43.merged.thUni.bam"
  querylabels <- "TDP43_AB"

  op <- "metagene_profile_of_TDP43_CLIP_signal"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)

  queryfiles <- "TDP43_AB_crosslinkregions.bed"
  querylabels <- "TDP43_AB_peaks"
  op <- "metagene_profile_of_TDP43_CLIP_peaks"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)


  beddata <- fread(queryfiles)
  dim(beddata)
  beddata <- beddata[V5 > 4]
  fwrite(beddata, "TDP43_AB_crosslinkregions_gt4.bed", sep="\t", col.names=F, quote=F)

  queryfiles <- "TDP43_AB_crosslinkregions_gt4.bed"
  querylabels <- "TDP43_AB_peaks_gt4"
  op <- "metagene_profile_of_TDP43_CLIP_peaks_gt4"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)

  ### correlation

  setwd(wd_cor)
  queryfiles <- file.path(wd_clip,"TDP43.merged.thUni.bam")
  querylabels <- "TDP43_AB_CLIP"
  centerfiles <- c(file.path(wd_chip, "GSE92026_ENCFF599JWP_conservative_idr_thresholded_peaks_hg19.bed"),
                   file.path(wd_clip, "TDP43_AB_crosslinkregions_gt4.bed"),
                   file.path(wd_r, "HEK293_R_loop.narrowPeak.bed"))

  centerlabels <- c("TDP43_ChIP_peaks", "TDP43_CLIP_peaks", "Rloop_peaks")
  op  <- "TDP43_clip_signal_around_ChIPandCLIPandRloop_peaks"
  system.time(
    plot_reference_locus(queryfiles, centerfiles, 500, c(-50,50), querylabels, centerlabels, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op)
  )



  queryfiles <- c(file.path(wd_chip,"GSE92026_ENCFF080OYX_signal_p-value_hg19.bigWig"))
  querylabels <- c("TDP43_chip")

  op  <- "TDP43_ChIP_signal_around_ChIPandCLIPandRloop_peaks"
  system.time(
    plot_reference_locus(queryfiles, centerfiles, 500, c(-50,50), querylabels, centerlabels, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op)
  )

  queryfiles <- c(file.path(wd_r,"GSM2550993_chipseq.HEK293.D210N_V5.p.bw"))
  inputfile <- file.path(wd_r,"GSM2551003_chipseq.HEK293.D210N_input.p.bw")
  querylabels <- c("Rloop_ratio_over_nput")


  op  <- "R_loop_ratioOverInput_around_ChIPandCLIPandRloop_peaks"
  system.time(
    plot_reference_locus(queryfiles, centerfiles, 500, c(-50,50), querylabels, centerlabels, ratioOverInput=T, inputfile=inputfile, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op)
  )

  queryfiles <- c(file.path(wd_chip,"GSE92026_ENCFF080OYX_signal_p-value_hg19.bigWig"),
                file.path(wd_clip,"TDP43.merged.thUni.bam"),
                file.path(wd_r,"GSM2550993_chipseq.HEK293.D210N_V5.p.bw"))
  querylabels <- c("TDP43_chip", "TDP43_AB_CLIP", "Rloop_V5")

  if(0){
    queryfiles <- queryfiles[3]
    querylabels <- querylabels[3]
    centerfiles <- centerfiles[3]
    centerlabels <- centerlabels[3]
  }

  #options(error = dump.frames)
  options(error = recover)
  op  <- "ChIPCLIPRloop_around_ChIPandCLIPandRloop_peaks_scaled"
  system.time(

    plot_reference_locus(queryfiles, centerfiles, 500, c(-50,50), querylabels, centerlabels, refPoint="center", stats=T, norm=F, scale=T, smo=T, Xlab="Peak Center", outPrefix=op)
  )

  op  <- "ChIPCLIPRloop_around_ChIPandCLIPandRloop_peaks_normalized"
  system.time(

    plot_reference_locus(queryfiles, centerfiles, 500, c(-50,50), querylabels, centerlabels, refPoint="center", stats=T, norm=T, scale=F, smo=F, Xlab="Peak Center", outPrefix=op)
  )
  dev.off()
  #load("last.dump.rda")
  #debugger()

  options(error = NULL)
}

## Nabeel Tetrahymena ChIPseq metagene plot ####
if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")
  gtffile <- "C:/GREENBLATT/Nabeel/Tetrahymena/2-Genome_GFF3.gff3"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/Tetrahymena/ChIPseq"
  setwd(wd)

  chipbed <- "Tetrahymena_ChIPseq_summits.bed"

  queryfiles <- list.files(pattern="sorted.bam$")#, "Input.SP1.A.thUni.bam")
  querylabels <- sapply(queryfiles, function(x) paste(unlist(strsplit(x, split="_"))[5:7], collapse="_"))#, "Input")

  centerfiles <- "C:/GREENBLATT/Nabeel/Tetrahymena/2-Genome_GFF3_gene.bed"
  centerlabels <- "gene"
  for(position in c("TSS", "TTS")){
    op <- paste0("individual_profile_of_H3_3_ChIPseq_at_",position)
    if(position == "TSS"){
      rp <- "start"
    }else{
      rp <- "end"
    }
    ## for individual
    plot_reference_locus(queryfiles, centerfiles, ext=c(-2000, 2000), hl=c(-2000,0), querylabels, centerlabels, inputfiles=NULL, ratioOverInput=F, fix_width=0, refPoint=rp, stats=T, norm=T, scale=F, smo=F, Xlab=position, outPrefix=op, genome="tetrahymena")
    ## for ratio
    #plot_reference_locus(queryfiles[1:2], centerfiles, c(-2000, 2000), c(-2000,0), querylabels[1:2], centerlabels, inputfiles=queryfiles[3:4], ratioOverInput=T, fix_width=0, refPoint=rp, stats=T, norm=T, scale=F, smo=F, Xlab=position, outPrefix=op, genome="tetrahymena")
  }


  op <- paste0("individual_metagene_profile_of_H3_3_ChIPseq")
  data_df <- plot_3parts_metagene(queryfiles, querylabels, txdb, meta=T, inputfiles=NULL, ratioOverInput=F, nbins=600, threeP=2000, fiveP=2000, smo=F, norm=T, CLIP_reads=F, fix_width=0, outPrefix=op, genome="tetrahymena")

  op <- paste0("individual_gene_profile_of_H3_3_ChIPseq")
  data_df <- plot_3parts_metagene(queryfiles, querylabels, txdb, meta=F, inputfiles=NULL, ratioOverInput=F, nbins=600, threeP=2000, fiveP=2000, smo=F, norm=T, CLIP_reads=F, fix_width=0, outPrefix=op, genome="tetrahymena")


  pdf("profile_of_H3_3_ChIPseq_at_gene_ends.pdf", height=6, width=8)
  plot_start_end_gtf(queryfiles, querylabels, txdb, featureName="gene", CLIP_reads=F, binsize=10, fix_width=0, norm=F, longest=T, ext=c(-2000, 1000, -1000, 2000), hl=c(-2000, 0, 0, 2000), randomize=T, stranded=T, genome="tetrahymena")
  dev.off()

  peakfile <- "Tetrahymena_ChIPseq_summits.bed"
  peaklabel <- "Tetrahymena_ChIPseq_summits"
  gtffile <- "C:/GREENBLATT/Nabeel/Tetrahymena/2-Genome_GFF3.gff3"
  genome <- "tetrahymena"
  fiveP <- 2000
  threeP <- 1000
  simple <- TRUE

  annotate_peaks(peakfile, peaklabel, gtffile, genome="tetrahymena", fiveP=2000, threeP=1000, simple=TRUE)
}

## Gio ####
if(1){
  options(error=recover)
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/GenomicPlot.R")
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

  wd <- "C:/GREENBLATT/Nabeel/Gio/bam"
  setwd(wd)

  ## bam correlation plot
  bamfiles <- list.files(pattern="^ZNF121.*.bam$")#, "Input.SP1.A.thUni.bam")
  bamlabels <- sapply(bamfiles, function(x) paste(unlist(strsplit(x, split="\\."))[1:2], collapse="_"))

  op="ZNF121_bam_correlations"
  plot_bam_correlation(bamfiles, bamlabels, binsize=1000000, outPrefix=op, genome="hg19")

  bamfiles <- list.files(pattern="^YTHDF2.*.bam$")#, "Input.SP1.A.thUni.bam")
  bamlabels <- sapply(bamfiles, function(x) paste(unlist(strsplit(x, split="\\."))[1:2], collapse="_"))

  op="YTHDF2_bam_correlations"
  plot_bam_correlation(bamfiles, bamlabels, binsize=1000000, outPrefix=op, genome="hg19")


  ## metagene plot
  queryfiles <- list.files(pattern=".*.merged.*bam$")#, "Input.SP1.A.thUni.bam")
  querylabels <- sapply(queryfiles, function(x) paste(unlist(strsplit(x, split="\\."))[1:2], collapse="_"))
  inputfiles <- list.files(pattern="^Input_.*bam$")
  inputlabels <- sapply(inputfiles, function(x) paste(unlist(strsplit(x, split="\\."))[1], collapse="_"))

  centerfiles <- c("combined_crosslink_m6AHek_top10percent.bed", "combined_crosslink_m6AHek_top20percent.bed")
  centerlabels <- c("m6AHek_top10pct", "m6AHek_top20pct")

  op <- "metagene_profile_of_ZNF121_YTHDF2_clip"
  data_df <- plot_5parts_metagene(queryfiles, querylabels, txdb, inputfiles=inputfiles, inputlabels=inputlabels, nbins=100, smo=T, norm=T,
                                  CLIP_reads=F, outPrefix=op, rm.outlier=F)

  op <- "metagene_profile_of_3parts_YTHDF2_clip"
  data_df <- plot_3parts_metagene(queryfiles=queryfiles[1], querylabels=querylabels[1], txdb=txdb, meta=T, longest=T, inputfiles=inputfiles[1],
                                  inputlabels=inputlabels[1], fiveP=1000, threeP=1000, nbins=100, smo=T, norm=T, CLIP_reads=F, fix_width=0,
                                  outPrefix=op, genome="hg19", rm.outlier=F)
  op <- "gene_profile_of_3parts_YTHDF2_clip"
  data_df <- plot_3parts_metagene(queryfiles=queryfiles[1], querylabels=querylabels[1], txdb=txdb, meta=F, longest=T, inputfiles=inputfiles[1],
                                  inputlabels=inputlabels[1], fiveP=1000, threeP=1000, nbins=100, smo=T, norm=T, CLIP_reads=F, fix_width=0,
                                  outPrefix=op, genome="hg19", rm.outlier=F)



  op <- "profile_of_ZNF121_YTHDF2_clip_at_m6A"
  data_list <- plot_reference_locus(queryfiles=queryfiles, centerfiles=centerfiles, ext=c(-500, 500), hl=c(-50,50), querylabels=querylabels,
                                    centerlabels=centerlabels, inputfiles=inputfiles, inputlabels=inputlabels, CLIP_reads=T,
                                    fix_width=0, refPoint="center", stats=T, norm=T, scale=F, smo=T, Xlab="Center", outPrefix=op,
                                    genome="hg19", rm.outlier = F)


  op <- "ratio_profile_of_ZNF121_YTHDF2_clip_at_m6A_random"
  plot_reference_locus_with_random(txdb, queryfiles, centerfiles, ext=c(-500,500), hl=c(-50,50), querylabels, centerlabels,
                                    smo=T, CLIP_reads=T, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=inputfiles,inputlabels=inputlabels, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                    genome="hg19", rm.outlier=F, n_random=1)



  wd <- "C:/GREENBLATT/Nabeel/Gio/pureclip_peaks"
  setwd(wd)

  peakfiles <- c("YTHDF2_AB_peaks_top10percent.bed", "YTHDF2_AB_peaks_top20percent.bed", "ZNF121_ABCDE_peaks_top10percent.bed", "ZNF121_ABCDE_peaks_top20percent.bed")
  peaklabels <- c("YTHDF2_top10pct", "YTHDF2_top20pct", "ZNF121_top10pct", "ZNF121_top20pct")
  names(peaklabels) <- peakfiles

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  genome <- "hg19"
  fiveP <- 1000
  threeP <- 1000
  simple <- FALSE
  RNA <- TRUE
  for(peakfile in peakfiles){
    peaklabel <- peaklabels[peakfile]
    annotate_peaks(peakfile, peaklabel, gtffile, genome, fiveP, threeP, simple, RNA)
  }


  op <- "metagene_profile_of_YTHDF2_ZNF121_peaks"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, querylabels=peaklabels, txdb, inputfiles=NULL, inputlabels=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)

  op <- "metagene_profile_of_YTHDF2_ZNF121_peaks_withIntron"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, querylabels=peaklabels, txdb, inputfiles=NULL, inputlabels=NULL, nbins=100, smo=F, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F, useIntron=T)

  op <- "gene_profile_of_YTHDF2_ZNF121_peaks"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, querylabels=peaklabels, txdb, inputfiles=NULL, inputlabels=NULL, nbins=100, smo=T, norm=T,
                                  meta=F, CLIP_reads=F, outPrefix=op, rm.outlier=F)


  op <- "ZNF121_peaks_around_YTHDF2_peaks"
  plot_reference_locus(queryfiles=peakfiles[3:4], centerfiles=peakfiles[1:2], ext=c(-500,500), hl=c(-50,50), querylabels=peaklabels[3:4], centerlabels=peaklabels[1:2], smo=T, CLIP_reads=FALSE,
                                   fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=NULL, inputlabels=NULL, stranded=TRUE,
                                   stats=T, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F)
  op <- "YTHDF2_peaks_around_ZNF121_peaks"
  plot_reference_locus(queryfiles=peakfiles[1:2], centerfiles=peakfiles[3:4], ext=c(-500,500), hl=c(-50,50), querylabels=peaklabels[1:2], centerlabels=peaklabels[3:4], smo=T, CLIP_reads=FALSE,
                       fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=NULL, inputlabels=NULL, stranded=TRUE,
                       stats=T, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F)


  op <- "ZNF121_peaks_around_YTHDF2_peaks_withRandom"
  plot_reference_locus_with_random(txdb, queryfiles=peakfiles[3:4], centerfiles=peakfiles[1:2], ext=c(-200,200), hl=c(-50,50), querylabels=peaklabels[3:4], centerlabels=peaklabels[1:2],
                                    smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op, n_random=3,
                                    genome="hg19", rm.outlier=F)
  op <- "YTHDF2_peaks_around_ZNF121_peaks_withRandom"
  plot_reference_locus_with_random(txdb, queryfiles=peakfiles[1:2], centerfiles=peakfiles[3:4], ext=c(-200,200), hl=c(-50,50), querylabels=peaklabels[1:2], centerlabels=peaklabels[3:4],
                                   smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                   inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op, n_random=3,
                                   genome="hg19", rm.outlier=F)



  centerfiles <- c("combined_crosslink_m6AHek.uniq_top10percent.bed", "combined_crosslink_m6AHek.uniq_top20percent.bed")
  centerlabels <- c("m6AHek_top10pct", "m6AHek_top20pct")
  names(centerlabels) <- centerfiles
  for(centerfile in centerfiles){
    centerlabel <- centerlabels[centerfile]
    annotate_peaks(centerfile, centerlabel, gtffile, genome, fiveP, threeP, simple, RNA)
  }


  op <- "YTHDF2_peak_over_m6A_with_random"
  plot_reference_locus_with_random(txdb, peakfiles[1:2], centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabels[1:2], centerlabels,
                                    smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                    genome="hg19", rm.outlier=F, n_random=3)


  op <- "ZNF121_peak_over_m6A_with_random"
  plot_reference_locus_with_random(txdb, peakfiles[3:4], centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabels[3:4], centerlabels,
                                    smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                    genome="hg19", rm.outlier=F, n_random=3)


  centerfiles <- c("combined_crosslink_YTHDF2_top10percent.bed", "combined_crosslink_YTHDF2_top20percent.bed")
  centerlabels <- c("YTHDF2cits_top10pct", "YTHDF2cits_top20pct")
  names(centerlabels) <- centerfiles
  for(centerfile in centerfiles){
    centerlabel <- centerlabels[centerfile]
    annotate_peaks(centerfile, centerlabel, gtffile, genome, fiveP, threeP, simple, RNA)
  }


  op <- "YTHDF2pureclip_over_YTHDF2cits_with_random"
  plot_reference_locus_with_random(txdb, peakfiles[1:2], centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabels[1:2], centerlabels,
                                   smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                   inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                   genome="hg19", rm.outlier=F, n_random=3)


  op <- "ZNF121pureclip_over_YTHDF2cits_with_random"
  plot_reference_locus_with_random(txdb, peakfiles[3:4], centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabels[3:4], centerlabels,
                                   smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                   inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                   genome="hg19", rm.outlier=F, n_random=3)


  peakfiles <- c("combined_crosslink_YTHDF2_top10percent.bed", "combined_crosslink_YTHDF2_top20percent.bed")
  peaklabels <- c("YTHDF2cits_top10pct", "YTHDF2cits_top20pct")
  centerfiles <- c("combined_crosslink_m6AHek_top10percent.bed", "combined_crosslink_m6AHek_top20percent.bed")
  centerlabels <- c("m6AHek_top10pct", "m6AHek_top20pct")

  op <- "YTHDF2cits_over_m6AHek_with_random"
  plot_reference_locus_with_random(txdb, peakfiles[1:2], centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabels[1:2], centerlabels,
                                   smo=FALSE, CLIP_reads=FALSE, fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                   inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                   genome="hg19", rm.outlier=F, n_random=3)


  ##
  setwd("C:/GREENBLATT/Nabeel/Gio/cits_peaks")

  peakfiles <- c("YTHDF2.A.thUni.tag.CITS.0.01.bed", "YTHDF2.A.thUni.tag.CIMS.0.01.bed",
                 "YTHDF2.B.thUni.tag.CITS.0.01.bed", "YTHDF2.B.thUni.tag.CIMS.0.01.bed",
                 "YTHDF2.A.thUni.crossLink_site.0.01.bed", "YTHDF2.B.thUni.crossLink_site.0.01.bed")
  peaklabels <- c("YTHDF2.A_CITS", "YTHDF2.A_CIMS", "YTHDF2.B_CITS", "YTHDF2.B_CIMS", "YTHDF2.A_CLS", "YTHDF2.B_CLS")
  names(peakfiles) <- peaklabels
  op <- "YTHDF2_replicates_peak_overlap"

  overlapBed(peakfiles, outPrefix=op, fix_width=20, fixPoint="center")

  peakfiles <- c("combined_CIMS_0.01_YTHDF2.recurring.bed", "combined_CIMS_0.01_ZNF121.recurring.bed",
                 "combined_CITS_0.01_YTHDF2.recurring.bed", "combined_CITS_0.01_ZNF121.recurring.bed",
                 "YTHDF2_AB_peaks_top10percent.bed", "ZNF121_ABCDE_peaks_top10percent.bed")
  peaklabels <- c("YTHDF2_CIMS", "ZNF121_CIMS", "YTHDF2_CITS", "ZNF121_CITS", "YTHDF2_pureclip", "ZNF121_pureclip")
  names(peakfiles) <- peaklabels
  op <- "YTHDF2_ZNF121_cits-pureclip_peak_overlap"
  overlapBed(peakfiles, outPrefix=op, fix_width=20, fixPoint="center")

  peakfiles <- c("combined_crosslink_0.01_YTHDF2.recurring.bed", "combined_crosslink_0.01_ZNF121.recurring.bed",
                 "YTHDF2_AB_peaks_top10percent.bed", "ZNF121_ABCDE_peaks_top10percent.bed")
  peaklabels <- c("YTHDF2_CLS", "ZNF121_CLS", "YTHDF2_pureclip", "ZNF121_pureclip")
  names(peakfiles) <- peaklabels
  op <- "YTHDF2_ZNF121_CLS-pureclip_peak_overlap"
  overlapBed(peakfiles, outPrefix=op, fix_width=20, fixPoint="center")


  peakfiles <- c("combined_CIMS_0.01_YTHDF2.recurring.bed", "combined_CIMS_0.01_ZNF121.recurring.bed",
                 "combined_CITS_0.01_YTHDF2.recurring.bed", "combined_CITS_0.01_ZNF121.recurring.bed",
                 "combined_crosslink_0.01_YTHDF2.recurring.bed", "combined_crosslink_0.01_ZNF121.recurring.bed")
  peaklabels <- c("YTHDF2_CIMS", "ZNF121_CIMS", "YTHDF2_CITS", "ZNF121_CITS", "YTHDF2_CLS", "ZNF121_CLS")
  names(peakfiles) <- peaklabels
  op <- "YTHDF2_ZNF121_peak_overlap"

  overlapBed(peakfiles[], outPrefix=op, fix_width=20, fixPoint="center")

  genome <- "hg19"
  fiveP <- 1000
  threeP <- 1000
  simple <- FALSE
  RNA <- TRUE

  lapply(c(5,6), function(x)annotate_peaks(peakfiles[x], gtffile, genome, fiveP, threeP, simple, RNA))

  op <- "metagene_profile_of_YTHDF2_ZNF121_peaks"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)

  peakfiles <- c("combined_CIMS_0.01_m6AHek.recurring.bed", "combined_CITS_0.01_m6AHek.recurring.bed", "combined_crosslink_0.01_m6AHek.recurring.bed",
                 "combined_CITS_0.01_m6AHek.merged.bed", "combined_CIMS_0.01_m6AHek.merged.bed")
  peaklabels <- c("m6AHek_recurring_CIMS", "m6AHek_recurring_CITS", "m6AHek_recurring_CLS", "m6AHek_CITS_all", "m6AHek_CIMS_all")
  names(peakfiles) <- peaklabels

  op <- "metagene_profile_of_m6AHek_peaks"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)
  lapply(seq_along(peakfiles), function(x)annotate_peaks(peakfiles[x], gtffile, genome, fiveP, threeP, simple, RNA))

  peakfiles <- c("combined_CIMS_0.01_YTHDF2.recurring.bed", "combined_CITS_0.01_YTHDF2.recurring.bed", "combined_crosslink_0.01_YTHDF2.recurring.bed",
                 "combined_CITS_0.01_YTHDF2.merged.bed", "combined_CIMS_0.01_YTHDF2.merged.bed")
  peaklabels <- c("YTHDF2_recurring_CIMS", "YTHDF2_recurring_CITS", "YTHDF2_recurring_CLS", "YTHDF2_CITS_all", "YTHDF2_CIMS_all")
  names(peakfiles) <- peaklabels

  op <- "metagene_profile_of_YTHDF2_peaks"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)
  lapply(seq_along(peakfiles), function(x)annotate_peaks(peakfiles[x], gtffile, genome, fiveP, threeP, simple, RNA))

  peakfiles <- c("combined_CIMS_0.01_ZNF121.recurring.bed", "combined_CITS_0.01_ZNF121.recurring.bed", "combined_crosslink_0.01_ZNF121.recurring.bed",
                 "combined_CITS_0.01_ZNF121.merged.bed", "combined_CIMS_0.01_ZNF121.merged.bed")
  peaklabels <- c("ZNF121_recurring_CIMS", "ZNF121_recurring_CITS", "ZNF121_recurring_CLS", "ZNF121_CITS_all", "ZNF121_CIMS_all")
  names(peakfiles) <- peaklabels

  op <- "metagene_profile_of_ZNF121_peaks"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)
  lapply(seq_along(peakfiles), function(x)annotate_peaks(peakfiles[x], gtffile, genome, fiveP, threeP, simple, RNA))



  peakfiles <- c("CITS95_crosslink_hotspot.bed", "CIMS95_crosslink_hotspot.bed")
  peaklabels <- c("CITS_hotspot", "CIMS_hotspot")
  names(peakfiles) <- peaklabels

  op <- "metagene_profile_of_CITS_hotspot"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)
  lapply(seq_along(peakfiles), function(x)annotate_peaks(peakfiles[x], gtffile, genome, fiveP, threeP, simple, RNA))
}

## Nabeel ZBTB48 over expression ####
if(1){
  options(error=recover)
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/misc_genomics_functions.R")
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

  wd <- "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS"
  setwd(wd)

  ## create m6Am file from merged cits file
  bedfile <- "combined_CITS_0.01_m6AHek.merged_A_slop0_out.bed"
  filter_bed_by_genomicFeature(bedfile, "utr5")


  genome <- "hg19"
  fiveP <- 1000
  threeP <- 1000
  simple <- FALSE
  RNA <- TRUE

  peakfiles <- c("combined_CITS_0.01_m6AZBTB48.recurring.bed", "combined_CITS_0.01_m6AFTO.recurring.bed", "combined_CITS_0.01_ZBTB48m6A.recurring.bed")
  peaklabels <- c("GFP", "FTO", "ZBTB48")
  names(peakfiles) <- peaklabels


  op <- "metagene_profile_of_m6A_CITS_peaks_recurring"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=F,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)


  op <- "gene_profile_of_m6A_CITS_peaks_recurring"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=F,
                                  meta=F, CLIP_reads=F, outPrefix=op, rm.outlier=F)



  op <- "CITS_peak_overlap"

  overlapBed(peakfiles, outPrefix=op, fix_width=20, fixPoint="center")
  lapply(seq_along(peakfiles), function(x)annotate_peaks(peakfiles[x], gtffile, genome, fiveP, threeP, simple, RNA))


  for(featureName in c("utr5", "utr3", "cds")){

    op <- paste0("m6A_CITS_peaks_at_", featureName)
    plot_start_end_feature(queryfiles=peakfiles, inputfiles=NULL, txdb, featureName, CLIP_reads=F, binsize=10, fix_width=0,
                           longest=T, ext=c(-200, 200, -200, 200), hl=c(-50, 50, -50, 50), randomize=F, stranded=T, norm=T, scale=F, smo=T, heatmap=F,
                           rm.outlier=F, genome="hg19", outprefix=op, useScore=F, useSizeFactor=F)

  }

  peakfiles <- c("combined_CITS_0.01_m6AZBTB48.merged.bed", "combined_CITS_0.01_m6AFTO.merged.bed", "combined_CITS_0.01_ZBTB48m6A.merged.bed")
  peaklabels <- c("GFP", "FTO", "ZBTB48")
  names(peakfiles) <- peaklabels


  op <- "metagene_profile_of_m6A_CITS_peaks_merged"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)


  op <- "gene_profile_of_m6A_CITS_peaks_merged"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=F, CLIP_reads=F, outPrefix=op, rm.outlier=F)

  wd <- "C:/GREENBLATT/Nabeel/ZBTB48/bam"
  setwd(wd)

  genome <- "hg19"
  fiveP <- 1000
  threeP <- 1000
  simple <- FALSE
  RNA <- TRUE

  peakfiles <- c("m6AZBTB48_merged.bam", "m6AFTO_merged.bam", "ZBTB48m6A_merged.bam")
  peaklabels <- c("GFP", "FTO", "ZBTB48")
  names(peakfiles) <- peaklabels


  op <- "metagene_profile_of_m6A_reads"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)


  op <- "gene_profile_of_m6A_reads"
  data_df <- plot_5parts_metagene(queryfiles=peakfiles, txdb=txdb, inputfiles=NULL, nbins=100, smo=T, norm=T,
                                  meta=F, CLIP_reads=F, outPrefix=op, rm.outlier=F)

  for(featureName in c("utr5", "utr3", "cds")){

    op <- paste0("m6A_reads_at_", featureName)
    plot_start_end_feature(queryfiles=peakfiles, inputfiles=NULL, txdb, featureName, CLIP_reads=T, binsize=10, fix_width=0,
                       longest=T, ext=c(-200, 200, -200, 200), hl=c(-50, 50, -50, 50), randomize=F, stranded=T, norm=T, scale=F, smo=T, heatmap=F,
                       rm.outlier=F, genome="hg19", outprefix=op, useScore=F, useSizeFactor=F)

  }
}

## Nabeel ZBTB48 iCLIP ####
options(error=recover)
source("C:/GREENBLATT/Rscripts/gProfilePlot/R/misc_genomics_functions.R")
gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
txdb <-  makeTxDbFromGFF(gtffile)

metaFeatures <- prepare_5parts_genomic_features(txdb, longest=TRUE, meta=TRUE, nbins=100, useIntron=FALSE, fiveP=1000, threeP=1000)
geneFeatures <- prepare_5parts_genomic_features(txdb, longest=TRUE, meta=FALSE, nbins=100, useIntron=FALSE, fiveP=1000, threeP=1000)

metaFeaturesLF <- prepare_5parts_genomic_features(txdb, longest=FALSE, meta=TRUE, nbins=100, useIntron=FALSE, fiveP=1000, threeP=1000)
geneFeaturesLF <- prepare_5parts_genomic_features(txdb, longest=FALSE, meta=FALSE, nbins=100, useIntron=FALSE, fiveP=1000, threeP=1000)

##### ZBTB48 ####
wd <- "C:/GREENBLATT/Nabeel/ZBTB48/ZBTB48_clip"
setwd(wd)

queryfiles <- c("ZBTB48.A.thUni.bam", "ZBTB48.B.thUni.bam", "ZBTB48.merged.thUni.bam")
querylabels <- c("ZBTB48_A", "ZBTB48_B", "ZBTB48_merged")
inputfiles <- "Input_ZBTB48.thUni.bam"
inputlabels <- "Input_ZBTB48"
names(inputfiles) <- inputlabels
names(queryfiles) <- querylabels


op <- "metagene_profile_of_ZBTB48_iCLIP_reads"
data_df <- plot_5parts_metagene(queryfiles, metaFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_ZBTB48_iCLIP_reads"
data_df <- plot_5parts_metagene(queryfiles, geneFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "metagene_profile_of_ZBTB48_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[3], metaFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_ZBTB48_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[3], geneFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

centerfiles <- c("C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE79577_FTO_clusters.bed",
                 "combined_CITS_0.01_FTO.recurring.bed")
centerlabels <- c("public_FTO_peaks", "recurring_FTO_peaks")
names(centerfiles) <- centerlabels

op <- "ZBTB48_iCLIP_ratioOverInput_around_FTO_peaks"
plot_reference_locus(queryfiles[3], centerfiles, ext=c(-200,200), hl=c(-50,50), smo=FALSE, CLIP_reads=FALSE, useSizeFactor=FALSE,
                                 fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=inputfiles, stranded=TRUE,
                                 stats=FALSE, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "ZBTB48_iCLIP_ratioOverInput_around_FTO_peaks_withRandom"
plot_reference_locus_with_random(queryfiles[3], centerfiles[2], txdb, ext=c(-200,200), hl=c(-0,0), useSizeFactor=T,
                                 smo=T,  CLIP_reads=FALSE,  fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                 inputfiles=inputfiles, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19",
                                 rm.outlier=F, n_random=1, stats.method="wilcox.test", useScore=FALSE)


centerfiles <- c("C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/CITS_m6A_12051.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/NIHMS870376_m6Am.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE180253_293T_m6Am_sites.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/putative_m6Am_sites_from_m6AHek_CITS.bed")
centerlabels <- c("public_m6A", "public_m6Am", "GSE180253_m6Am", "Hek_m6A", "Hek_m6Am")
names(centerfiles) <- centerlabels

op <- "ZBTB48_iCLIP_ratioOverInput_around_m6A_sites"
plot_reference_locus(queryfiles[3], centerfiles, ext=c(-200,200), hl=c(-0,0), smo=TRUE, CLIP_reads=FALSE, useSizeFactor=TRUE,
                     fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=inputfiles, stranded=TRUE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "ZBTB48_iCLIP_signal_around_m6A_sites_withRandom"
plot_reference_locus_with_random(queryfiles[3], centerfiles[4:5], txdb, ext=c(-200,200), hl=c(-50,50), shade=F, useSizeFactor=T,
                                             smo=T,  CLIP_reads=FALSE,  fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                             inputfiles=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19",
                                             rm.outlier=F, n_random=10, stats.method="wilcox.test", useScore=FALSE)


##### FTO ####
wd <- "C:/GREENBLATT/Nabeel/ZBTB48/ZBTB48_clip"
setwd(wd)

queryfiles <- c("FTO_A.thUni.bam", "FTO_B.thUni.bam", "FTO.S3.thUni.bam", "FTO.S4.thUni.bam", "FTO.merged.thUni.bam")
querylabels <- c("FTO_A", "FTO_B", "FTO_S3", "FTO_S4", "FTO_merged")
inputfiles <- "Input_FTO_A.thUni.bam"
inputlabels <- "Input_FTO"
names(inputfiles) <- inputlabels
names(queryfiles) <- querylabels


op <- "metagene_profile_of_FTO_iCLIP_reads"
data_df <- plot_5parts_metagene(queryfiles, metaFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_FTO_iCLIP_reads"
data_df <- plot_5parts_metagene(queryfiles, geneFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "metagene_profile_of_FTO_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[5], metaFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_FTO_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[5], geneFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

centerfiles <- c("C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/CITS_m6A_12051.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/NIHMS870376_m6Am.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE180253_293T_m6Am_sites.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/putative_m6Am_sites_from_m6AHek_CITS.bed")
centerlabels <- c("public_m6A", "public_m6Am", "GSE180253_m6Am", "Hek_m6A", "Hek_m6Am")

names(centerfiles) <- centerlabels

op <- "FTO_iCLIP_ratioOverInput_around_m6A_sites"
plot_reference_locus(queryfiles[5], centerfiles, ext=c(-200,200), hl=c(-50,50), smo=TRUE, CLIP_reads=FALSE, useSizeFactor=TRUE,
                     fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=inputfiles, stranded=TRUE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "FTO_iCLIP_ratioOverInput_around_m6A_sites_withRandom"
plot_reference_locus_with_random(queryfiles[5], centerfiles[4:5], txdb, ext=c(-200,200), hl=c(-0,0), useSizeFactor=T,
                                 smo=T,  CLIP_reads=FALSE,  fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                 inputfiles=inputfiles, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19",
                                 rm.outlier=F, n_random=1, stats.method="wilcox.test", useScore=FALSE)


queryfiles <- c("ZBTB48.merged.thUni.bam", "FTO.merged.thUni.bam")
querylabels <- c("ZBTB48_merged", "FTO_merged")
inputfiles <- c("Input_ZBTB48.thUni.bam", "Input_FTO_A.thUni.bam")
inputlabels <- c("Input_ZBTB48", "Input_FTO")
names(inputfiles) <- inputlabels
names(queryfiles) <- querylabels

op <- "metagene_profile_of_ZBTB48FTO_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles, metaFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_ZBTB48FTO_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles, geneFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)


### Downloaded FTO ####
wd <- "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A"
setwd(wd)
queryfiles <- c("FTO_1.bam", "FTO_2.bam", "FTO_3.bam", "FTO_merge.bam")
querylabels <- c("FTO_1", "FTO_2", "FTO_3", "FTO_merged")
inputfiles <- "FTO_Input1_2.bam"
inputlabels <- "Input_FTO"
names(inputfiles) <- inputlabels
names(queryfiles) <- querylabels


op <- "metagene_profile_of_publicFTO_iCLIP_reads"
data_df <- plot_5parts_metagene(queryfiles, metaFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_publicFTO_iCLIP_reads"
data_df <- plot_5parts_metagene(queryfiles, geneFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "metagene_profile_of_publicFTO_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[4], metaFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = T)

op <- "gene_profile_of_publicFTO_iCLIP_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[4], geneFeatures, inputfiles, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = T)

centerfiles <- c("C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/CITS_m6A_12051.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/NIHMS870376_m6Am.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE180253_293T_m6Am_sites.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/putative_m6Am_sites_from_m6AHek_CITS.bed")
centerlabels <- c("public_m6A", "public_m6Am", "GSE180253_m6Am", "Hek_m6A", "Hek_m6Am")
names(centerfiles) <- centerlabels


op <- "publicFTO_CLIP_ratioOverInput_around_m6A_sites"
plot_reference_locus(queryfiles[4], centerfiles, ext=c(-200,200), hl=c(-50,50), smo=TRUE, CLIP_reads=FALSE, useSizeFactor=TRUE,
                     fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=inputfiles, stranded=TRUE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "Hekm6A_sites_around_publicm6A_sites"
plot_reference_locus(centerfiles, centerfiles, ext=c(-100,100), hl=c(-50,50), smo=FALSE, CLIP_reads=FALSE, useSizeFactor=FALSE,
                     fix_width=0, norm=F, binsize=1, refPoint="center", Xlab="Center", inputfiles=NULL, stranded=TRUE,verbose=FALSE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "metagene_profile_of_m6A_sites"
data_df <- plot_5parts_metagene(centerfiles, metaFeatures, inputfiles=NULL, smo=F, norm=F,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = F)

op <- "gene_profile_of_m6A_sites"
data_df <- plot_5parts_metagene(centerfiles, geneFeatures, inputfiles=NULL, smo=F, norm=F,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = F)


A_TSS_count <- count_overlap_at_TSS(centerfiles, txdb, TSSflank=10, CLIP_reads=FALSE, longest=TRUE)

### FTO in siZBTB48 ####
wd <- "C:/GREENBLATT/Nabeel/ZBTB48/FTO_clip_siZBTB48"
setwd(wd)

queryfiles <- list.files(pattern="\\.bam$")
querylabels <- gsub(".thUni.bam", "", queryfiles, fixed=T)
names(queryfiles) <- querylabels


op <- "metagene_profile_of_FTO_in_siZBTB48_reads"
data_df <- plot_5parts_metagene(queryfiles[c(1,2,4,5)], metaFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "gene_profile_of_FTO_in_siZBTB48_reads"
data_df <- plot_5parts_metagene(queryfiles[c(1,2,4,5)], geneFeatures, inputfiles=NULL, smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F)

op <- "metagene_profile_of_FTO_in_siZBTB48_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[c(3,6)], metaFeatures, inputfiles=queryfiles[c(7,8)], smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = T)

op <- "gene_profile_of_FTO_in_siZBTB48_ratioOverInput"
data_df <- plot_5parts_metagene(queryfiles[c(3,6)], geneFeatures, inputfiles=queryfiles[c(7,8)], smo=T, norm=T,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = T)

centerfiles <- c("C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/CITS_m6A_12051.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/NIHMS870376_m6Am.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE180253_293T_m6Am_sites.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/putative_m6Am_sites_from_m6AHek_CITS.bed")
centerlabels <- c("public_m6A", "public_m6Am", "GSE180253_m6Am", "Hek_m6A", "Hek_m6Am")
names(centerfiles) <- centerlabels

op <- "FTO_iCLIP_in_siZBTB48_signal_around_m6A_sites"
plot_reference_locus(queryfiles[c(3,6)], centerfiles[4:5], ext=c(-200,200), hl=c(-50,50), shade=F,  smo=TRUE, CLIP_reads=FALSE, useSizeFactor=TRUE,
                     fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center", inputfiles=NULL, stranded=TRUE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "FTO_iCLIP_in_siZBTB48_signal_around_m6A_sites_withRandom"
plot_reference_locus_with_random(queryfiles[c(3,6)], centerfiles[4:5], txdb, ext=c(-200,200), hl=c(-0,0), useSizeFactor=T,
                                 smo=T,  CLIP_reads=FALSE,  fix_width=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                 inputfiles=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19",
                                 rm.outlier=F, n_random=1, stats.method="wilcox.test", useScore=FALSE)

### deseq2 for ZBTB48 #############
if(1){
   source("C:/GREENBLATT/Rscripts/gProfilePlot/R/misc_genomics_functions.R")

   gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
   txdb <-  makeTxDbFromGFF(gtffile)


   wd <- "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/bam"
   setwd(wd)

   queryfiles <- c("m6AZBTB48_S1.thUni.bam", "m6AZBTB48_S2.thUni.bam",
                   "ZBTB48m6A.S1.thUni.bam", "ZBTB48m6A.S2.thUni.bam",
                   "m6AFTO_S1.thUni.bam", "m6AFTO_S2.thUni.bam")
   querylabels <- c("GFP_S1", "GFP_S2", "ZBTB48_S1", "ZBTB48_S2", "FTO_S1", "FTO_S2")
   names(queryfiles) <- querylabels
   sampleNames <- c("GFP", "FTO", "ZBTB48")
   subsets <- list("ZBTB48" = c("GFP_S1", "GFP_S2", "ZBTB48_S1", "ZBTB48_S2"), "FTO"=c("GFP_S1", "GFP_S2", "FTO_S1", "FTO_S2"))

   plot_feature_overlap_count(queryfiles, NULL, txdb, sampleNames, longest=T, CLIP_reads=T, flank=5, norm=T, subsets=subsets, input_type="reads")


   peakfile <- c("C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/putative_m6Am_sites_from_m6AHek_CITS.bed")
   peaklabel <- c("Hek_m6A", "Hek_m6Am")
   names(peakfile) <- peaklabel

   for(i in 1:2){
      plot_feature_overlap_count(queryfiles, peakfile[i], txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=T, subsets=subsets, input_type="reads")
   }


   ## bam correlation
   op <- "m6A_bams_correlation"
   plot_bam_correlation(queryfiles, outPrefix = op)


   wd <- "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS"
   setwd(wd)
   queryfiles <- c("m6AZBTB48_S1.thUni.tag.CITS.bed", "m6AZBTB48_S2.thUni.tag.CITS.bed",
                   "ZBTB48m6A.S1.thUni.tag.CITS.bed", "ZBTB48m6A.S2.thUni.tag.CITS.bed",
                   "m6AFTO_S1.thUni.tag.CITS.bed", "m6AFTO_S2.thUni.tag.CITS.bed")

   querylabels <- c("GFP_S1", "GFP_S2", "ZBTB48_S1", "ZBTB48_S2", "FTO_S1", "FTO_S2")
   names(queryfiles) <- querylabels
   sampleNames <- c("GFP", "FTO", "ZBTB48")
   subsets <- list("ZBTB48" = c("GFP_S1", "GFP_S2", "ZBTB48_S1", "ZBTB48_S2"), "FTO"=c("GFP_S1", "GFP_S2", "FTO_S1", "FTO_S2"))

   plot_feature_overlap_count(queryfiles, NULL, txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=F,subsets=subsets, input_type="peaks_individual")

   for(i in 1:2){
      plot_feature_overlap_count(queryfiles, peakfile[i], txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=F, subsets=subsets, input_type="peaks_individual")
   }

   peakfiles <- c("combined_CITS_0.01_m6AZBTB48.merged.bed",
                  "combined_CITS_0.01_m6AFTO.merged.bed",
                  "combined_CITS_0.01_ZBTB48m6A.merged.bed")
   peaklabels <- c("GFP_peak", "FTO_peak", "ZBTB48_peak")
   names(peakfiles) <- peaklabels
   sampleNames <- c("GFP", "FTO", "ZBTB48")

   plot_feature_overlap_count(queryfiles=peakfiles, NULL, txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=F,subsets=NULL, input_type="peaks_merged")
   for(i in 1:2){
      plot_feature_overlap_count(queryfiles=peakfiles, peakfile[i], txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=F,subsets=NULL, input_type="peaks_merged")
   }

   peakfiles <- c("combined_CITS_0.01_m6AZBTB48.recurring.bed",
                  "combined_CITS_0.01_m6AFTO.recurring.bed",
                  "combined_CITS_0.01_ZBTB48m6A.recurring.bed")
   peaklabels <- c("GFP_peak", "FTO_peak", "ZBTB48_peak")
   names(peakfiles) <- peaklabels
   sampleNames <- c("GFP", "FTO", "ZBTB48")

   plot_feature_overlap_count(queryfiles=peakfiles, NULL, txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=F, subsets=NULL, input_type="peaks_recurring")
   for(i in 1:2){
      plot_feature_overlap_count(queryfiles=peakfiles, peakfile[i], txdb, sampleNames, longest=F, CLIP_reads=T, flank=5, norm=F, subsets=NULL, input_type="peaks_recurring")
   }
}


### Derive m6Am from Hekm6A ####
wd <- "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS"
setwd(wd)

tx_hg19 <- transcripts(txdb,  use.name=F)
TSS_annot_hg19 <- resize(tx_hg19, width=1, fix="start", use.names=FALSE, ignore.strand=FALSE)
TSS_annot_hg19_bed <- annoGR2DF(TSS_annot_hg19) %>%
   select(c(chr, start, end, name=tx_name, score=width, strand)) %>%
   mutate(end=as.integer(end)+1)

hg19_TSS_file <- "hg19_TSS.bed"
names(hg19_TSS_file) <- "annot_TSS"
write.table(TSS_annot_hg19_bed, hg19_TSS_file, col.names=F, row.names = F, sep="\t", quote=F)

#chainFile <- "C:/GREENBLATT/genomic_feature/hg38ToHg19.over.chain"
#hg38bed <- "refTSS_v3.1_human_coordinate.hg38.bed"

#exp_TSS_file <- liftOverBed(chainFile, hg38bed)
#names(exp_TSS_file) <- "exp_TSS"

exp_TSS_file <- "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/FANTOM5_Hek293_hg19.ctss.bed"
names(exp_TSS_file) <- "FANTOM5_TSS"
TSS_exp_hg19 <- import.bed(exp_TSS_file)

cagefile <- "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/FANTOM5_Hek293.hg19.nobarcode.bam"
names(cagefile) <- "FANTOM5_CAGE"

m6Amfile_GSE180253 <- "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE180253_293T_m6Am_sites.bed"
names(m6Amfile_GSE180253) <- "GSE180253_m6Am"

## create m6Am file from merged cits file
mergedfile <- "combined_CITS_0.01_m6AHek.merged.bed"
names(mergedfile) <- "m6AHek_CITS_merged"
bedfile <- "combined_CITS_0.01_m6AHek.merged_A_slopl2r2.bed"
names(bedfile) <- "m6AHek_CITS_merged_atA"
m6Amfile <- "putative_m6Am_sites_from_m6AHek_CITS.bed"
names(m6Amfile) <- "putative_m6Am"
m6Afile <- "combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed"
names(m6Afile) <- "putative_m6A"
#annotate_peaks(bedfile, gtffile, genome="hg19", fiveP=1000, threeP=1000, simple=FALSE, RNA=TRUE)
op <- "metagene_profile_of_CITS_merged"
data_df <- plot_5parts_metagene(c(mergedfile, bedfile), metaFeatures, inputfiles=NULL, smo=F, norm=F,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = F)

utr5_filtered <- filter_bed_by_genomicFeature(bedfile, featureName="utr5", txdb=txdb, resizeFraction=0.25, longest=TRUE)
m6A <- import.bed(m6Afile)
m6Am <- filter_by_nonoverlaps_stranded(utr5_filtered$bed, m6A)
utr5_filtered_bed <- annoGR2DF(m6Am) %>%
   select(c(chr, start, end, name, score, strand)) %>%
   mutate(start=as.integer(start)-1) %>%
   mutate(chr=as.character(chr))
dim(utr5_filtered_bed)

write.table(utr5_filtered_bed, m6Amfile, col.names=F, row.names = F, sep="\t", quote=F)
op <- "metagene_profile_of_CITS_m6Am"
data_df <- plot_5parts_metagene(m6Amfile, metaFeatures, inputfiles=NULL, smo=F, norm=F,
                                CLIP_reads=F, outPrefix=op, rm.outlier=F, useSizeFactor = F)

#annotate_peaks(m6Amfile, gtffile, genome="hg19", fiveP=1000, threeP=1000, simple=FALSE, RNA=TRUE)

peakfiles <- c(exp_TSS_file, hg19_TSS_file, cagefile)
names(peakfiles) <- c("FANTOM5_TSS", "annot_TSS", "FANTOM5_CAGE")

centerfiles <- c("C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/CITS_m6A_12051.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/NIHMS870376_m6Am.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/downloaded_FTO_M6A/GSE180253_293T_m6Am_sites.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/combined_CITS_0.01_m6AHek.merged_DRACH_slopl2r2.bed",
                 "C:/GREENBLATT/Nabeel/ZBTB48/M6A_clip_oeZBTB48/CITS/putative_m6Am_sites_from_m6AHek_CITS.bed")
centerlabels <- c("public_m6A", "public_m6Am", "GSE180253_m6Am", "Hek_m6A", "Hek_m6Am")
names(centerfiles) <- centerlabels

op <- "distance_between_putative_m6A-m6Am_and_TSS"
plot_reference_locus(queryfiles=peakfiles, centerfiles=centerfiles, ext=c(-50,50), hl=c(-10,10), smo=FALSE, CLIP_reads=FALSE, useSizeFactor=FALSE, verbose=FALSE,
                     fix_width=0, norm=FALSE, binsize=1, refPoint="center", Xlab="Center", inputfiles=NULL, stranded=TRUE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")

op <- "distance_between_TSS_and_putative_m6A-m6Am"
plot_reference_locus(queryfiles=centerfiles, centerfiles=peakfiles[2], ext=c(-50, 50), hl=c(-10,10), smo=FALSE, CLIP_reads=FALSE, useSizeFactor=FALSE, verbose=FALSE,
                     fix_width=0, norm=FALSE, binsize=1, refPoint="center", Xlab="Center", inputfiles=NULL, stranded=TRUE,
                     heatmap=FALSE, scale=FALSE, outPrefix=op, genome="hg19", rm.outlier=F, useScore=F, stats.method="wilcox.test")


op <- "venn_between_m6AHek_CITS_derived_peaks"
overlapBed(c(mergedfile, bedfile, m6Afile, m6Amfile), outPrefix=op, fix_width=10L, fixPoint="center", pairOnly=F, genome="hg19")

op <- "venn_between_m6A_and_m6Am"
overlapBed(c(m6Afile, m6Amfile), outPrefix=op, fix_width=0L, fixPoint="center", pairOnly=T, genome="hg19")


### create m6Am annotation table
peak_gene_association <- annotate_peaks(m6Amfile, gtffile, genome="hg19", fiveP=1000, threeP=1000, simple=FALSE, RNA=TRUE)

expTSS_filtered <- filter_bed_by_genomicFeature(m6Amfile, featureName = "expTSS", featureGr=TSS_exp_hg19, maxgap=0)
annotTSS_filtered <- filter_bed_by_genomicFeature(m6Amfile, featureName = "annotTSS", featureGr=TSS_annot_hg19, maxgap=0)

expTSS_filtered_bed <- annoGR2DF(expTSS_filtered$bed) %>%
   select(c(chr, start, end, name, score, strand))%>%
   mutate(start=as.integer(start)-1) %>%
   mutate(chr=as.character(chr)) %>%
   mutate(expTSS=rep("Yes", length(expTSS_filtered$bed)))

annotTSS_filtered_bed <- annoGR2DF(annotTSS_filtered$bed) %>%
   select(c(chr, start, end, name, score, strand))%>%
   mutate(start=as.integer(start)-1) %>%
   mutate(chr=as.character(chr)) %>%
   mutate(annotTSS=rep("Yes", length(annotTSS_filtered$bed)))

integrated_bed <- utr5_filtered_bed %>%
   mutate(UTR=rep("5'UTR", nrow(utr5_filtered_bed))) %>%
   left_join(expTSS_filtered_bed) %>%
   left_join(annotTSS_filtered_bed) %>%
   replace_na(list(expTSS = 'No', annotTSS = "No"))

out_table <- left_join(integrated_bed, peak_gene_association,
                   by=c("chr"="chrPeak", "start"="startPeak", "end"="endPeak", "strand"="strandPeak")) %>%
   select(-c(id, widthPeak))

write.table(out_table, "overlap_putative_m6Am_with_TSS.tab", sep="\t", row.names=F, quote=F)

dim(out_table)
tmp <- out_table %>%
   filter(annotTSS=="Yes")
dim(tmp)


## SP1 PAU gene expression ####

abundance <- read.table("C:/GREENBLATT/Jingwen/RNAseq/gene_expression_abundance.tab")
upPAU <- read.table("C:/GREENBLATT/Jingwen/GO_results_25_75/PAU_analysis_siSP1_up_list_dutr.bed")
downPAU <- read.table("C:/GREENBLATT/Jingwen/GO_results_25_75/PAU_analysis_siSP1_down_list_dutr.bed")

SP1_TPM <- abundance %>%
   mutate(geneName=rownames(abundance)) %>%
   mutate(SP1=(siSP1_a+siSP1_b)/2) %>%
   select(c(geneName, SP1))
head(SP1_TPM)
head(upPAU)

colnames(upPAU) <- colnames(downPAU) <- c("chr", "start","end", "id", "length","strand","geneId", "geneName")
up_PAU <- left_join(upPAU, SP1_TPM) %>%
   mutate(SP1 = replace_na(SP1,0)) %>%
   mutate(above1TPM = ifelse(SP1 > 1, 1, 0)) %>%
   mutate(above2TPM = ifelse(SP1 > 2, 1, 0))

down_PAU <- left_join(downPAU, SP1_TPM) %>%
   mutate(SP1 = replace_na(SP1,0)) %>%
   mutate(above1TPM = ifelse(SP1 > 1, 1, 0)) %>%
   mutate(above2TPM = ifelse(SP1 > 2, 1, 0))


summary(up_PAU)
summary(down_PAU)

summary_table <- as.data.frame(matrix(rep(0,8), nrow=2, ncol=4, dimnames=list(c("PAU_up", "PAU_down"), c("Count", "MeanTPM", "Fraction>1TPM", "Fraction>2TPM")))) %>%
   mutate(Count = c(nrow(up_PAU), nrow(down_PAU))) %>%
   mutate(MeanTPM = c(mean(up_PAU$SP1), mean(down_PAU$SP1))) %>%
   mutate(`Fraction>1TPM` = c(mean(up_PAU$above1TPM), mean(down_PAU$above1TPM))) %>%
   mutate(`Fraction>2TPM` = c(mean(up_PAU$above2TPM), mean(down_PAU$above2TPM)))

write.table(summary_table, "C:/GREENBLATT/Jingwen/RNAseq/gene_expression_abundance_in_siSP1_PAU.tab", col.names=NA, sep="\t")
