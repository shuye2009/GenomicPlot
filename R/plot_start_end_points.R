rm(list=ls(all=TRUE))

source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")

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
#bedfiles <- c("combined_crosslink_YTHDF2.uniq.bed") #"MR_Ac4c_merged.thUni.bed",
bedfiles <- c("TR_Ac4c_merged.thUni.bam", "MR_Ac4c_merged.thUni.bam", "NAT10_merged.thUni.bam")
bedlabels <- c("TR", "MR", "NAT10")
pdf("TR_MR_NAT10_merged.bam_metagene.pdf")
data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, sn=0, norm=TRUE, CLIP_reads=F)
dev.off()

featureName <- "cds"
hl <- c(-100, 0, 0, 100) ## boundaries for highlight, 1,2 for start, 3,4 for end
sn <- 100000
ran <- FALSE

pdf("TR_Ac4c_merged.thUni_intron.pdf")
plot_start_end_gtf(bedfiles, bedlabels, txdb, featureName, ext, hl, sn, ran)
dev.off()
}

## down-sampling noChange ZNF281 ####

if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")
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
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")

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

  bedfiles <- c("ZNF281.merged.bam", "NUDT21.merged.bam", "CPSF7.merged.bam", "CSTF2_GSM917676_275_+_density.bw", "CSTF2T_merged_+_density.bw", "Im68_merged_+_density.bw" )
  bedlabels <- c("ZNF281", "NUDT21", "CPSF7", "CSTF2", "CSTF2T", "CPSF6")

  if(1){
    pdf(paste("ZNF281_CPSF_CSTF_clip_read_starts_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=T, stats=F, scale=T)
    dev.off()
  }
  if(1){
    pdf(paste("ZNF281_CPSF_CSTF_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=F, stats=F, scale=T)
    dev.off()
  }

  bedfiles <- c("ZNF281.merged.bam", "NUDT21.merged.bam", "CPSF7.merged.bam", "CPSF160.merged.bam", "CPSF2.merged.bam", "FIP1.merged.bam")
  bedlabels <- c("ZNF281", "NUDT21", "CPSF7", "CPSF160", "CPSF2", "FIP1")

  if(1){
    pdf(paste("ZNF281_CPSF_clip_read_starts_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=T, stats=F, scale=T)
    dev.off()
  }
  if(1){
    pdf(paste("ZNF281_CPSF_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", norm=F, smo=T, CLIP_reads=F, stats=F, scale=T)
    dev.off()
  }

  bedfiles <- c("ZNF281.ABC_crosslinkregions_wInputwCL.bed", "NUDT21_CDE_crosslinkregions_wInputwCL.bed", "CPSF7_AB_crosslinkregions_wInputwCL.bed", "CPSF160_AB_crosslinkregions_wInputwCL.bed", "CPSF2_AB_crosslinkregions_wInputwCL.bed", "FIP1_AB_crosslinkregions_wInputwCL.bed")
  bedlabels <- c("ZNF281", "NUDT21", "CPSF7", "CPSF160", "CPSF2", "FIP1")

  if(1){
    pdf(paste("ZNF281_CPSF_Peaks_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", smo=T, norm=F, stats=F, scale=T)
    dev.off()
  }

### PAU regulated cleavage sites
  for(utr in c("putr", "dutr")){
    centerfiles_hg38 <- c(file.path(qapa_dir, paste("PAU_analysis_siZNF281_up_list_",utr,".bed", sep="")),
                     file.path(qapa_dir, paste("PAU_analysis_siZNF281_down_list_",utr,".bed", sep="")),
                      file.path(qapa_dir, paste("PAU_analysis_siZNF281_noChange_list_",utr,"_sub.bed", sep="")))
    centerfiles <- NULL
    for(bedFile in centerfiles_hg38){
      newFile <- liftOverBed(chainFile, bedFile)
      centerfiles <- c(centerfiles, newFile)
    }
    centerlabels <- c("Up", "Down", "noChange")

    hl <- c(-50, 50) ## c(-50, 50), c(-100, 0)
    hln <- paste(hl, collapse="-")
    ext <- 200
    # reads


    if(1){
      bedfiles <- c("ZNF281.merged.bam")
      bedlabels <- c("ZNF281")
      Inputfile <- "Input.SP1.A.thUni.bam"

      op <- paste("ZNF281_merged_clip_ratioOverInput_around_cleavage_sites", hln, utr, sep="_")
      plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=Inputfile, ratioOverInput=T, outPrefix=op)

    }


    if(1){
      bedfiles <- c("NUDT21.merged.bam")
      bedlabels <- c("NUDT21")
      Inputfile <- "Input.NUDT21.A.thUni.bam"

      op <- paste("NUDT21_merged_clip_ratioOverInput_around_cleavage_sites", hln, utr, sep="_")
      plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=Inputfile, ratioOverInput=T, outPrefix=op)

    }


    if(1){
      bedfiles <- c("CPSF7.merged.bam")
      bedlabels <- c("CPSF7")
      Inputfile <- "Input.CPA.A.thUni.bam"

      op <- paste("CPSF7_merged_clip_ratioOverInput_around_cleavage_sites", hln, utr, sep="_")
      plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=Inputfile, ratioOverInput=T, outPrefix=op)

    }

    if(1){
      bedFiles <- c("CSTF2_GSM917676_275_+_density.bw", "CSTF2T_merged_+_density.bw", "Im68_merged_+_density.bw")
      names(bedFiles) <- c("CSTF2" ,"CSTF2T", "CPSF6")
      for(bedlabels in names(bedFiles)){
        bedfiles <- bedFiles[bedlabels]

        op <- paste(bedlabels,"_clip_reads_around_cleavage_sites", hln, utr, sep="_")
        plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", heatmap=T, stats=T, norm=F, scale=F, smo=T, inputfile=NULL, ratioOverInput=F, outPrefix=op)

      }
    }
  }

  ### reads over peaks
  centerfiles <- c("ZNF281.ABC_crosslinkregions_wInputwCL.bed")
  centerlabels <- c("ZNF281_pureclip_peaks")


  bedfiles <- c("ZNF281.merged.bam")
  bedlabels <- c("ZNF281ABC")

  if(1){
    pdf(paste("ZNF281_merged_clip_reads_around_pureclip_peaks.pdf", sep=""), height=8, width=12)
    plot_reference_locus_with_random(txdb=txdb, bedfiles=bedfiles, centerfiles=centerfiles, ext=ext, hl=hl, bedlabels=bedlabels, centerlabels=centerlabels, CLIP_reads = F)
    dev.off()
  }


  ### non-regulated putr dutr sutr

  centerfiles_hg38 <- c(file.path(qapa_dir, "PAU_analysis_putr.bed"),
                        file.path(qapa_dir, "PAU_analysis_dutr.bed"),
                        file.path(qapa_dir, "PAU_analysis_sutr.bed"))
                        #file.path(qapa_dir, "UTR3_coordinates.bed"))

  centerfiles <- NULL
  for(bedFile in centerfiles_hg38){
    newFile <- liftOverBed(chainFile, bedFile)
    centerfiles <- c(centerfiles, newFile)
  }
  centerlabels <- c("putr", "dutr", "sutr")

  hl <- c(-50, 50) ## c(-50, 50), c(-100, 0)

  bedfiles <- c("CPSF7.merged.bam")
  bedlabels <- c("CPSF7AB")
  if(1){
    pdf(paste("CPSF7_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  bedfiles <- c("Im68_GSM917665_273_+_density.wig", "Im68B_GSM917664_274_+_density.wig")
  bedlabels <- c("Im68A", "Im68B")
  if(1){
    pdf(paste("Im68_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }


  # reads
  bedfiles <- c("ZNF281.A.thUni.bam", "ZNF281.B.thUni.bam", "ZNF281.C.thUni.bam")
  bedlabels <- c("ZNF281A", "ZNF281B", "ZNF281C")
  Inputfile <- "Input.SP1.A.thUni.bam"

  if(0){
    pdf(paste("ZNF281_individual_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }
  if(0){
    pdf(paste("ZNF281_individual_clip_ratioOverInput_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", inputfile=Inputfile, ratioOverInput = T)
    dev.off()
  }

  bedfiles <- c("ZNF281.AC.merged.bam", "ZNF281.merged.bam")
  bedlabels <- c("ZNF281AC", "ZNF281ABC")

  if(1){
    pdf(paste("ZNF281_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  if(1){
    pdf(paste("ZNF281_merged_clip_ratioOverInput_around_cleavage_site.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site", inputfile=Inputfile, ratioOverInput = T)
    dev.off()
  }

  bedfiles <- c("NUDT21.merged.bam")
  bedlabels <- c("NUDT21CDE")

  if(1){
    pdf(paste("NTDT21_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }


  bedfiles <- c("CSTF2_GSM917676_275_+_density.bw", "CSTF2T_A_GSM917677_283_+_density.bw")
  bedlabels <- c("CSTF2", "CSTF2T")

  if(1){
    pdf(paste("CSTF2_CSTF2T_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  bedfiles <- c("CSTF2T_merged_+_density.bw")
  bedlabels <- c("CSTF2TAB")

  if(1){
    pdf(paste("CSTF2T_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
    dev.off()
  }

  bedfiles <- c("Im68_merged_+_density.bw")
  bedlabels <- c("Im68AB")

  if(1){
    pdf(paste("Im68_merged_clip_reads_around_cleavage_sites.pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
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

    bedfiles <- c("ZNF281.merged.bam")#, "Input.SP1.A.thUni.bam")
    bedlabels <- c("ZNF281ABC")#, "Input")

    pdf("metagene_profile_of_ZNF281_merged_bam_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, CLIP_reads=F, smo=T)
    dev.off()

    bedfiles <- c("NUDT21.merged.bam")#, "Input.NUDT21.A.thUni.bam")
    bedlabels <- c("NUDT21CDE")#, "Input")

    pdf("metagene_profile_of_NUDT21_merged_bam_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, CLIP_reads=F, smo=T)
    dev.off()

    bedfiles <- c("CPSF7.merged.bam")#, "Input.NUDT21.A.thUni.bam")
    bedlabels <- c("CPSF7AB")#, "Input")

    pdf("metagene_profile_of_CPSF7_merged_bam_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, CLIP_reads=F, smo=T)
    dev.off()


    pdf("metagene_profile_of_ZNF281_merged_bam_500bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=500, CLIP_reads=F)
    dev.off()

    bedfiles <- c("ZNF281_over_Input_DESeq2_insignificant_crossLink_site.bed", "ZNF281_over_Input_DESeq2_significant_crossLink_site.bed")
    bedlabels <- c("insignificant_peaks", "significant_peaks")

    pdf("metagene_profile_of_ZNF281_significant_peaks_100bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F)
    dev.off()

    pdf("metagene_profile_of_ZNF281_significant_peaks_500bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=500, CLIP_reads=F)
    dev.off()

    bedfiles <- c("ZNF281.ABC_crosslinkregions_wInputwCL.bed")
    bedlabels <- c("pureclip_peaks")

    pdf("metagene_profile_of_ZNF281_pureclip_peaks_100bin_scaled.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T, scale=T)
    dev.off()

    pdf("metagene_profile_of_ZNF281_pureclip_peaks_500bin.pdf", height=8, width=12)
    data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=500, CLIP_reads=F, smo=T)
    dev.off()

    bedfiles <- c("CSTF2_GSM917676_275_+_density.bw")
    bedlabels <- c("CSTF2")

    if(1){
      pdf(paste("metagene_profile_of_CSTF2_downloaded_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    bedfiles <- c("Im68_GSM917665_273_+_density.bw", "Im68B_GSM917664_274_+_density.bw")
    bedlabels <- c("Im68A", "Im68B")

    if(1){
      pdf(paste("metagene_profile_of_CFIm68_downloaded_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    if(1){
      pdf(paste("metagene_profile_of_Input_bam_100bin.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(Inputfile, "Input", txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    bedfiles <- c("CSTF2T_merged_+_density.bw")
    bedlabels <- c("CSTF2TAB")

    if(1){
      pdf(paste("metagene_profile_of_CSTF2T_merged_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
      dev.off()
    }

    bedfiles <- c("Im68_merged_+_density.bw")
    bedlabels <- c("Im68AB")

    if(1){
      pdf(paste("metagene_profile_of_CFIm68_merged_wig.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, norm=FALSE, nbins=100, CLIP_reads=F, smo=T)
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
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

  chainFile <- "C:/GREENBLATT/genomic_feature/hg38ToHg19.over.chain"

  ext <- 200
  wd <- "C:/GREENBLATT/Nabeel/ZBTB7A"
  setwd(wd)

  qapa_dir <- "C:/GREENBLATT/Nabeel/ZBTB7A"

  bedfile <- file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_up_list_","putr",".bed", sep=""))
  for(utr in c("putr", "dutr")){
    centerfiles_hg38 <- c(file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_up_list_",utr,".bed", sep="")),
                          file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_down_list_",utr,".bed", sep="")),
                          file.path(qapa_dir, paste("PAU_analysis_ZBTB7A_noChange_list_",utr,"_sub.bed", sep="")))
    centerfiles <- NULL
    for(bedFile in centerfiles_hg38){
      newFile <- liftOverBed(chainFile, bedFile)
      centerfiles <- c(centerfiles, newFile)
    }
    centerlabels <- c("Up", "Down", "noChange")

    hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

    # reads
    bedfiles <- c("ZBTB7A.A.thUni.bam", "ZBTB7A.B.thUni.bam", "ZBTB7A.ABCD.merged.bam")
    bedlabels <- c("ZBTB7A.A", "ZBTB7A.B", "ZBTB7A.ABCD")
    pdf(paste("ZBTB7A_clip_reads_around_cleavage_sites3_",utr,".pdf", sep=""), height=8, width=12)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", Xlab="Cleavage site")
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

  bedfiles <- c("ZBTB7A.A.thUni.bam", "ZBTB7A.B.thUni.bam", "ZBTB7A.ABCD.merged.bam")
  bedlabels <- c("ZBTB7A.A", "ZBTB7A.B", "ZBTB7A.ABCD")
  for(refp in c("start", "end")){
    hl <- hls[[refp]]
    op <- paste0("ZBTB7A_clip_reads_around_regulated_intron_", refp)
    plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint=refp, Xlab=refp, stat=T, scale=F, smo=T, norm=F, outPrefix=op)
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
  bedfiles <- c("NAT10_merged.thUni.bam")
  bedlabels <- c("NAT10")

  pdf("reads_distance_to_center_AC4C_with_random_narrow.pdf", height=8, width=12)
  plot_reference_locus_with_random(txdb, bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, binsize=10, refPoint="center", CLIP_reads=FALSE, stranded=TRUE)
  dev.off()

  ## peaks
  bedfiles <- c("NAT10_AB_crosslinkregions_wInputwCL.bed")
  bedlabels <- c("NAT10")
  pdf("peaks_distance_to_center_AC4C.pdf", height=8, width=12)
  plot_reference_locus(bedfiles, centerfiles, ext, hl, sn, bedlabels, centerlabels)
  dev.off()

  pdf("peaks_distance_to_center_AC4C_with_random.pdf", height=8, width=12)
  plot_reference_locus_with_random(txdb, bedfiles, centerfiles, ext, hl, bedlabels, centerlabels)
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
  bedfiles <- c("TR_Ac4c_merged.thUni.bed", "MR_Ac4c_merged.thUni.bed")
  bedlabels <- c("TR_AC4C", "MR_AC4C")

  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  pdf("distance_to_ExonCenter_100_AC4C.pdf", height=8, width=12)
  plot_reference_locus(bedfiles, centerfiles, ext, hl, sn, bedlabels, centerlabels)
  dev.off()
}

## Nabeel m6A ####
if(1){

  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 400
  wd <- "C:/GREENBLATT/Nabeel/m6A/cits_peak"
  setwd(wd)
  bedfiles <- c("m6A_FTO_merged.thUni.crossLink_site.0.05.bed", "m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed")
  centerlabels <- c("Hek")
  centerfiles <- c("m6A_Hek_merged.thUni.crossLink_site.0.05.bed")
  bedlabels <- c("FTO", "ZBTB48")
  sn <- 100000
  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  pdf("distance_to_Hek_peaks.pdf", height=8, width=12)
  plot_reference_locus(bedfiles, centerfiles, ext, hl, sn, bedlabels, centerlabels)
  dev.off()

  ## plot for gene feature boundaries
  bedfiles <- c("m6A_FTO_merged.thUni.crossLink_site.0.05.bed",
                "m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed",
                "m6A_Hek_merged.thUni.crossLink_site.0.05.bed")
  bedlabels <- c("FTO", "ZBTB48", "Hek")

  bedfiles <- c("m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed")
  bedlabels <- c("ZBTB48")

  pdf("metagene_profile_of_peaks_longest_transcript.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, sn=0, norm=FALSE, longest=TRUE)
  dev.off()

  bedfiles <- c("m6A_FTO_merged.thUni.bed")
  bedlabels <- c("FTO_reads")
  pdf("metagene_profile_of_merged_FTO_bed.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, sn=0, norm=F)
  dev.off()

  for(featureName in c("utr5", "utr3", "cds")){

    ran <- ifelse(featureName %in% c("exon", "intron"), TRUE, FALSE)
    ## hl is boundaries for highlight, 1,2 for start, 3,4 for end
    hl <- c(0, 100, -100, 0)
    if(featureName %in% c("intron")) hl <- c(-100, 0, 0, 100)

    pdf(paste("m6A_FTO_reads_at_", gsub("'", "",featureName), ".pdf", sep=""), height=8, width=12)
    plot_start_end_gtf(bedfiles, bedlabels, txdb, featureName, ext, hl, sn, ran)
    dev.off()
  }
}

## plot metagene profile m6A #############
if(1){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/m6A/TSN"
  setwd(wd)
  bedfiles <- "m6A_FTO_merged.thUni.bam"
  bedlabels <- "FTO"
  bedfiles <- c("m6A_Hek_merged.thUni.bam", "m6A_ZBTB48_merged.thUni.bam", "m6A_FTO_merged.thUni.bam")
  bedlabels <- c("Hek", "ZBTB48", "FTO")
  pdf("metagene_profile_of_merged_bam_100bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, sn=0, CLIP_reads=F)
  dev.off()

  pdf("metagene_profile_of_merged_bam_500bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=500, sn=0, CLIP_reads=F)
  dev.off()
}

## plot metagene profile for NAT10 reads #############
if(1){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/kristie/NAT10/reads_metagene_profile"
  setwd(wd)
  bedfiles <- "NAT10_merged.thUni.bam"
  bedlabels <- "NAT10_merged"

  pdf("metagene_profile_of_NAT10_merged_bam_100bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, sn=0, CLIP_reads=F)
  dev.off()

  pdf("metagene_profile_of_NAT10_merged_bam_500bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=500, sn=0, CLIP_reads=F)
  dev.off()
}

## plot metagene profile for NAT10 peaks #############
if(1){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/kristie/NAT10/peak_metagene_profile"
  setwd(wd)
  bedfiles <- "NAT10_AB_crosslinkregions_wInputwCL.bed"
  bedlabels <- "NAT10_merged_peaks"

  pdf("metagene_profile_of_NAT10_merged_peaks_100bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, sn=0, CLIP_reads=F)
  dev.off()

  pdf("metagene_profile_of_NAT10_merged_peaks_500bin.pdf", height=8, width=12)
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=500, sn=0, CLIP_reads=F)
  dev.off()
}

## deseq2 for m6A #############
if(1){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/m6A/TSN"
  setwd(wd)

  bedfiles <- c("m6A_Hek_S1.thUni.bam", "m6A_Hek_S2.thUni.bam",
                "m6A_ZBTB48_S1.thUni.bam", "m6A_ZBTB48_S2.thUni.bam",
                "m6A_FTO_S1.thUni.bam", "m6A_FTO_S2.thUni.bam")
  bedlabels <- c("Hek_S1", "Hek_S2", "ZBTB48_S1", "ZBTB48_S2", "FTO_S1", "FTO_S2")


  feature_count_list <- count_overlap_in_features(bedfiles, bedlabels, txdb, longest=T)
  subsets <- list("ZBTB48" = c("Hek_S1", "Hek_S2", "ZBTB48_S1", "ZBTB48_S2"), "FTO"=c("Hek_S1", "Hek_S2", "FTO_S1", "FTO_S2"))

  for(flank in c(1:3)){
    data_dfs <- count_overlap_at_TSS(bedfiles, bedlabels, txdb, ext=flank, longest=T)

    pdf(paste("lfc_at_TSNflank",flank, "_and_transcript_forReadStarts.pdf",sep=""), height=8, width=12)

    for(subject in names(subsets)){
      subset <- subsets[[subject]]
      s2c <- data.frame(samples=subset, groups=c(rep(c("Hek", subject), times=c(2,2))))
      contrasts <- data.frame(groups=c(rep("groups", 1)), treat=subject, control=c(rep("Hek",1)))
      data_table_TSN <- data_dfs[, subset]
      data_table_TSN <- data_table_TSN[apply(data_table_TSN,1,sum)>0,]
      project <- paste("TSNflank", flank, sep="")

      res_TSN <- run_DESeq2(data_table_TSN, s2c, project, contrasts)

      for(featureName in c("utr5", "cds", "utr3", "exon")){

        data_table_feature <- feature_count_list[[featureName]][, subset]
        data_table_feature <- data_table_feature[apply(data_table_feature,1,sum)>0,]
        featureName <- ifelse(featureName == "exon", "transcript", featureName)

        res_feature <- run_DESeq2(data_table_feature, s2c, featureName, contrasts)

        log2FoldChange <- c(res_TSN[[1]]$log2FoldChange, res_feature[[1]]$log2FoldChange)


        feature <- rep(c("TSN", featureName), times=c(length(res_TSN[[1]]$log2FoldChange), length(res_feature[[1]]$log2FoldChange)))

        out_df <- data.frame(log2FoldChange, feature) %>%
          mutate(feature=factor(feature, levels=c("TSN", featureName)))

        median_df <- aggregate(out_df$log2FoldChange, list(out_df$feature), median)
        colnames(median_df) <- c("feature", "median")
        d <- ggplot(out_df, aes(x=log2FoldChange, color=feature)) +
          geom_density() +
          geom_vline(data=median_df, aes(xintercept=median, color=feature), linetype="dotted", size=1) +
          theme_classic()

        #print(d)

        count_df <- aggregate(out_df$log2FoldChange, list(out_df$feature), length)
        colnames(count_df) <- c("feature", "count")
        p <- ggplot(out_df, aes(x=feature, y=log2FoldChange, fill=feature)) +
          geom_boxplot(notch=TRUE) +
          theme_classic() + ggtitle(subject) +
          theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=14)) +
          theme(legend.position = "none") +
          labs(y=expression(paste(log[2], " (Overexpression/WT signal intensity)"))) +
          scale_x_discrete(limits=c("TSN", featureName)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A")))) +
          theme(axis.text.x = element_text(face="bold", size=14, color="black")) +
          geom_text(data = count_df, aes(x=feature, y = min(out_df$log2FoldChange)*1.1, label = paste("n=",count,sep="")))

        p1 <- add_pval(p, pairs = list(c(1, 2)), test='t.test')

        #print(p1)

        outp <- plot_grid(p1, d, ncol = 2, rel_widths = c(1,1))
        print(outp)
      }

    }
    dev.off()
  }
}



## did not work out
if(0){
library(Guitar)
Gtxdb <- makeGuitarTxdb(txdb=txdb, pltTxType="mrna")
bedfiles <- c("m6A_FTO_merged.thUni.crossLink_site.0.05.bed")
site <- blocks(import(bedfiles))
sitesGRanges <- samplePoints(list(site),
                             stSampleNum = 5,
                             stAmblguity = 5,
                             pltTxType = c("mrna"),
                             stSampleModle = "Equidistance",
                             mapFilterTranscript = FALSE,
                             Gtxdb)
sitesNormlize <- normalize(sitesGRanges,
                           Gtxdb,
                           txType = "mrna",
                           overlapIndex = 1,
                           siteLengthIndex = 1)

## plot for gene feature boundaries
bedfiles <- c("m6A_FTO_merged.thUni.crossLink_site.0.05.bed",
              "m6A_ZBTB48_merged.thUni.crossLink_site.0.05.bed",
              "m6A_Hek_merged.thUni.crossLink_site.0.05.bed")
stGRangelist<-list()
for (i in 1:length(bedfiles)) {
  stGRangelist[[i]] <-  blocks(import(bedfiles[[i]]))
}
names(stGRangelist) <- c("FTO", "ZBTB48", "Hek")
GuitarPlot(txTxdb = txdb,
           stGRangeLists = stGRangelist,
           pltTxType = "mrna",
           txPrimaryOnly = TRUE,
           stGroupName = c("FTO", "ZBTB48", "Hek"))
stBedFiles <- list("m6A_FTO_merged.thUni.bed")
GuitarPlot(txGenomeVer = "hg19",
           stBedFiles = stBedFiles,
           pltTxType = "mrna",
           txPrimaryOnly = TRUE,
           txGuitarTxdbSaveFile = "hg19_Gtxdb.dat"
)

}


if(0){
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)
  ext <- 400
  wd <- "C:/GREENBLATT/Nabeel/kristie/clip_reads_around_junction"
  setwd(wd)
  #bedfiles <- c("combined_crosslink_YTHDF2.uniq.bed") #"MR_Ac4c_merged.thUni.bed",
  bedfiles <- c("TR_Ac4c_merged.thUni.bed", "MR_Ac4c_merged.thUni.bed")
  bedlabels <- c("TR_AC4C", "MR_AC4C")

  featureName <- "5'utr"
  hl <- c(-100, 0, 0, 100) ## boundaries for highlight, 1,2 for start, 3,4 for end
  sn <- 100000
  ran <- TRUE

  pdf("TR_Ac4c_merged.thUni_intron.pdf")
  plot_start_end_gtf(bedfiles, txdb, featureName, ext, hl, sn, ran)
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
  bedfiles <- c("TR_Ac4c_merged.thUni.bed", "MR_Ac4c_merged.thUni.bed")
  bedlabels <- c("TR_AC4C", "MR_AC4C")
  sn <- 10000
  hl <- c(-50, 50) ## boundaries for highlight, 1,2 for start, 3,4 for end

  pdf("distance_to_center_AC4C.pdf")
  plot_reference_locus(bedfiles, centerfiles, ext, hl, sn, bedlabels, centerlabels)
  dev.off()
}


## Jingwen ChIPseq reads #############
if(1){


  TFs <- c("SP1", "CTCF", "MAZ", "OSR2", "ZBTB48", "ZNF281")

  cl <- start_parallel(length(TFs))

  parLapply(cl, TFs, function(TF){

    source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")
    gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
    txdb <-  makeTxDbFromGFF(gtffile)
    wd <- "C:/GREENBLATT/Jingwen/ChIPseq"
    qapa <- "C:/GREENBLATT/Jingwen/GO_results_25_75"
    setwd(wd)

    hl_list <- list(c(-500, -250), c(-250,0), c(0,250), c(250,500))## c(-50, 50), c(-100, 0)

    if(0){
      bedfiles <- paste(TF, ".sorted.bam", sep="")
      pdf(paste("metagene_profile_of_",TF, "_ChIPseq_100bin.pdf", sep=""), height=8, width=12)
      data_df <- plot_5parts_metagene(bedfiles, TF, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F)
      dev.off()

      if(TF %in% c("MAZ")){
        bedfiles <- paste(TF, ".merged_clip.bam", sep="")
        pdf(paste("metagene_profile_of_",TF, "_CLIP_100bin.pdf", sep=""), height=8, width=12)
        data_df <- plot_5parts_metagene(bedfiles, TF, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F)
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
          bedfiles <- paste(TF, ".sorted.bam", sep="")
          bedlabels <- paste(TF, "ChIPseq")
          outPrefix <- paste(TF,i, "_ChIPseq_reads_around_",utr,"_cleavage_sites", sep="")
          plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", stats=T, norm=F, smo=F, Xlab="Cleavage site",outPrefix=outPrefix)



          if(TF %in% c("MAZ")){
            bedfiles <- paste(TF, ".merged_clip.bam", sep="")
            bedlabels <- paste(TF, "CLIP")
            outPrefix <- paste(TF,i, "_CLIP_reads_around_",utr,"_cleavage_sites", sep="")
            plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="end", stats=T, norm=F, smo=F, Xlab="Cleavage site",outPrefix=outPrefix)


          }
        }


  ## find the transcripts that correspond to the utr3
        if(utr == "putr"){
          centerfiles <- unlist(lapply(centerfiles, function(x)utr3_to_transcript(x, txdb)))
          ext <- 1000

          if(1){
            bedfiles <- paste(TF, ".sorted.bam", sep="")
            bedlabels <- paste(TF, "ChIPseq")
            outPrefix <- paste(TF,i, "_ChIPseq_reads_around_TSS", sep="")
            plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="start", stats=T, norm=F, smo=F, Xlab="TSS",outPrefix=outPrefix)



            if(TF %in% c("MAZ")){
              bedfiles <- paste(TF, ".merged_clip.bam", sep="")
              bedlabels <- paste(TF, "CLIP")
              outPrefix <- paste(TF,i, "_CLIP_reads_around_TSS", sep="")
              plot_reference_locus(bedfiles, centerfiles, ext, hl, bedlabels, centerlabels, refPoint="start", stats=T, norm=F, smo=F, Xlab="TSS",outPrefix=outPrefix)

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
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  txdb <-  makeTxDbFromGFF(gtffile)

### R_loop
  wd_r <- "C:/GREENBLATT/Nabeel/TDP43/R_loop"
  wd_chip <- "C:/GREENBLATT/Nabeel/TDP43/TDP43_ChIP"
  wd_clip <- "C:/GREENBLATT/Nabeel/TDP43/CLIP"
  wd_cor <- "C:/GREENBLATT/Nabeel/TDP43/Correlations"
  setwd(wd_r)

  bedfiles <- list.files(wd, pattern=".p.bw")
  bedlabels <- paste(rep(c("V5", "Input"), 3), rep(c("rep1", "rep2", "rep3"),each=2), sep="_")

  op <- "metagene_profile_of_R_loop"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op)


  bedfiles <- c("GSM2550993_chipseq.HEK293.D210N_V5.p.bw", "GSM2551003_chipseq.HEK293.D210N_input.p.bw")
  bedlabels <- c("V5", "Input")

  op <- "metagene_profile_of_R_loop_merged"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op)

  centerfiles <- c("level9_R-loop_zones.bed", "level5_R-loop_zones.bed")
  centerlabels <- c("R_loop_level9", "R_loop_level5")

  op  <- "R_loop_signal_around_level9level5_peaks"
  system.time(
    plot_reference_locus(bedfiles, centerfiles, 500, c(-50,50), bedlabels, centerlabels, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op, genome="hg19")
  )


### ChIPseq

  setwd(wd_chip)

  bedfiles <- list.files(wd, pattern=".bed")
  bedlabels <- c("all_Peak", "conservative_peak", "optimal_peak")

  op <- "metagene_profile_of_TDP43_ChIPseq_peak_Norm"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op)


  bedfiles <- list.files(wd, pattern=".bigWig")
  bedlabels <- c("pValue", "FCOC")

  op <- "metagene_profile_of_TDP43_ChIPseq_signal"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)


  ## CLIP

  setwd(wd_clip)
  bedfiles <- "TDP43.merged.thUni.bam"
  bedlabels <- "TDP43_AB"

  op <- "metagene_profile_of_TDP43_CLIP_signal"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)

  bedfiles <- "TDP43_AB_crosslinkregions.bed"
  bedlabels <- "TDP43_AB_peaks"
  op <- "metagene_profile_of_TDP43_CLIP_peaks"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)


  beddata <- fread(bedfiles)
  dim(beddata)
  beddata <- beddata[V5 > 4]
  fwrite(beddata, "TDP43_AB_crosslinkregions_gt4.bed", sep="\t", col.names=F, quote=F)

  bedfiles <- "TDP43_AB_crosslinkregions_gt4.bed"
  bedlabels <- "TDP43_AB_peaks_gt4"
  op <- "metagene_profile_of_TDP43_CLIP_peaks_gt4"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, nbins=100, smo=T, norm=F, CLIP_reads=F, outPrefix=op)

  ### correlation

  setwd(wd_cor)
  bedfiles <- file.path(wd_clip,"TDP43.merged.thUni.bam")
  bedlabels <- "TDP43_AB_CLIP"
  centerfiles <- c(file.path(wd_chip, "GSE92026_ENCFF599JWP_conservative_idr_thresholded_peaks_hg19.bed"),
                   file.path(wd_clip, "TDP43_AB_crosslinkregions_gt4.bed"),
                   file.path(wd_r, "HEK293_R_loop.narrowPeak.bed"))

  centerlabels <- c("TDP43_ChIP_peaks", "TDP43_CLIP_peaks", "Rloop_peaks")
  op  <- "TDP43_clip_signal_around_ChIPandCLIPandRloop_peaks"
  system.time(
    plot_reference_locus(bedfiles, centerfiles, 500, c(-50,50), bedlabels, centerlabels, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op)
  )



  bedfiles <- c(file.path(wd_chip,"GSE92026_ENCFF080OYX_signal_p-value_hg19.bigWig"))
  bedlabels <- c("TDP43_chip")

  op  <- "TDP43_ChIP_signal_around_ChIPandCLIPandRloop_peaks"
  system.time(
    plot_reference_locus(bedfiles, centerfiles, 500, c(-50,50), bedlabels, centerlabels, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op)
  )

  bedfiles <- c(file.path(wd_r,"GSM2550993_chipseq.HEK293.D210N_V5.p.bw"))
  inputfile <- file.path(wd_r,"GSM2551003_chipseq.HEK293.D210N_input.p.bw")
  bedlabels <- c("Rloop_ratio_over_nput")


  op  <- "R_loop_ratioOverInput_around_ChIPandCLIPandRloop_peaks"
  system.time(
    plot_reference_locus(bedfiles, centerfiles, 500, c(-50,50), bedlabels, centerlabels, ratioOverInput=T, inputfile=inputfile, refPoint="center", stats=T, norm=F, smo=F, Xlab="Peak Center",outPrefix=op)
  )

  bedfiles <- c(file.path(wd_chip,"GSE92026_ENCFF080OYX_signal_p-value_hg19.bigWig"),
                file.path(wd_clip,"TDP43.merged.thUni.bam"),
                file.path(wd_r,"GSM2550993_chipseq.HEK293.D210N_V5.p.bw"))
  bedlabels <- c("TDP43_chip", "TDP43_AB_CLIP", "Rloop_V5")

  if(0){
    bedfiles <- bedfiles[3]
    bedlabels <- bedlabels[3]
    centerfiles <- centerfiles[3]
    centerlabels <- centerlabels[3]
  }

  #options(error = dump.frames)
  options(error = recover)
  op  <- "ChIPCLIPRloop_around_ChIPandCLIPandRloop_peaks_scaled"
  system.time(

    plot_reference_locus(bedfiles, centerfiles, 500, c(-50,50), bedlabels, centerlabels, refPoint="center", stats=T, norm=F, scale=T, smo=T, Xlab="Peak Center", outPrefix=op)
  )

  op  <- "ChIPCLIPRloop_around_ChIPandCLIPandRloop_peaks_normalized"
  system.time(

    plot_reference_locus(bedfiles, centerfiles, 500, c(-50,50), bedlabels, centerlabels, refPoint="center", stats=T, norm=T, scale=F, smo=F, Xlab="Peak Center", outPrefix=op)
  )
  dev.off()
  #load("last.dump.rda")
  #debugger()

  options(error = NULL)
}

## Nabeel Tetrahymena ChIPseq metagene plot ####
if(1){
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")
  gtffile <- "C:/GREENBLATT/Nabeel/Tetrahymena/2-Genome_GFF3.gff3"
  txdb <-  makeTxDbFromGFF(gtffile)
  wd <- "C:/GREENBLATT/Nabeel/Tetrahymena/ChIPseq"
  setwd(wd)

  chipbed <- "Tetrahymena_ChIPseq_summits.bed"

  bedfiles <- list.files(pattern="sorted.bam$")#, "Input.SP1.A.thUni.bam")
  bedlabels <- sapply(bedfiles, function(x) paste(unlist(strsplit(x, split="_"))[5:7], collapse="_"))#, "Input")

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
    plot_reference_locus(bedfiles, centerfiles, ext=c(-2000, 2000), hl=c(-2000,0), bedlabels, centerlabels, inputfiles=NULL, ratioOverInput=F, extend_reads=0, refPoint=rp, stats=T, norm=T, scale=F, smo=F, Xlab=position, outPrefix=op, genome="tetrahymena")
    ## for ratio
    #plot_reference_locus(bedfiles[1:2], centerfiles, c(-2000, 2000), c(-2000,0), bedlabels[1:2], centerlabels, inputfiles=bedfiles[3:4], ratioOverInput=T, extend_reads=0, refPoint=rp, stats=T, norm=T, scale=F, smo=F, Xlab=position, outPrefix=op, genome="tetrahymena")
  }


  op <- paste0("individual_metagene_profile_of_H3_3_ChIPseq")
  data_df <- plot_3parts_metagene(bedfiles, bedlabels, txdb, meta=T, inputfiles=NULL, ratioOverInput=F, nbins=600, threeP=2000, fiveP=2000, smo=F, norm=T, CLIP_reads=F, extend_reads=0, outPrefix=op, genome="tetrahymena")

  op <- paste0("individual_gene_profile_of_H3_3_ChIPseq")
  data_df <- plot_3parts_metagene(bedfiles, bedlabels, txdb, meta=F, inputfiles=NULL, ratioOverInput=F, nbins=600, threeP=2000, fiveP=2000, smo=F, norm=T, CLIP_reads=F, extend_reads=0, outPrefix=op, genome="tetrahymena")


  pdf("profile_of_H3_3_ChIPseq_at_gene_ends.pdf", height=6, width=8)
  plot_start_end_gtf(bedfiles, bedlabels, txdb, featureName="gene", CLIP_reads=F, binsize=10, extend_reads=0, norm=F, longest=T, ext=c(-2000, 1000, -1000, 2000), hl=c(-2000, 0, 0, 2000), randomize=T, stranded=T, genome="tetrahymena")
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
  source("C:/GREENBLATT/Rscripts/gProfilePlot/R/plot_start_end_points_lib.R")
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
  bedfiles <- list.files(pattern=".*.merged.*bam$")#, "Input.SP1.A.thUni.bam")
  bedlabels <- sapply(bedfiles, function(x) paste(unlist(strsplit(x, split="\\."))[1:2], collapse="_"))
  inputfiles <- list.files(pattern="^Input_.*bam$")
  inputlabels <- sapply(inputfiles, function(x) paste(unlist(strsplit(x, split="\\."))[1], collapse="_"))



  op <- "metagene_profile_of_ZNF121_YTHDF2_clip"
  data_df <- plot_5parts_metagene(bedfiles, bedlabels, txdb, inputfiles=inputfiles, inputlabels=inputlabels, nbins=100, smo=T, norm=T, CLIP_reads=F, outPrefix=op, rm.outlier=F)

  centerfiles <- "m6A_Hek_merged.thUni.crossLink_site.0.05.bed"
  centerlabels <- "m6A_Hek"

  op <- "profile_of_ZNF121_YTHDF2_clip_at_m6A"
  data_list <- plot_reference_locus(bedfiles=bedfiles, centerfiles=centerfiles, ext=c(-500, 500), hl=c(-50,50), bedlabels=bedlabels, centerlabels=centerlabels, inputfiles=inputfiles, inputlabels=inputlabels, CLIP_reads=T, extend_reads=0, refPoint="center", stats=T, norm=T, scale=F, smo=T, Xlab="Center", outPrefix=op, genome="hg19", rm.outlier = F)


  op <- "ratio_profile_of_ZNF121_YTHDF2_clip_at_m6A_random"
  plot_reference_locus_with_random1(txdb, bedfiles, centerfiles, ext=c(-200,200), hl=c(-50,50), bedlabels, centerlabels,
                                    smo=FALSE, CLIP_reads=FALSE, extend_reads=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=inputfiles,inputlabels=inputlabels, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                    genome="hg19", rm.outlier=F, n_random=1)

  peakfile <- "m6A_Hek_merged.thUni.crossLink_site.0.05.bed"
  peaklabel <- "m6A_Hek"
  gtffile <- "C:/GREENBLATT/genomic_feature/gencode.v19.annotation.gtf"
  genome <- "hg19"
  fiveP <- 1000
  threeP <- 1000
  simple <- FALSE

  annotate_peaks(peakfile, peaklabel, gtffile, genome, fiveP, threeP, simple)

  wd <- "C:/GREENBLATT/Nabeel/Gio/pureclip_peaks"
  setwd(wd)

  centerfiles <- "m6A_Hek_merged.thUni.crossLink_site.0.05.bed"
  centerlabels <- "m6A_Hek"

  peakfile <- "YTHDF2_AB_crosslinkregions_wInputwCL.bed"
  peaklabel <- "YTHDF2_AB"
  annotate_peaks(peakfile, peaklabel, gtffile, genome, fiveP, threeP, simple, RNA=T)

  op <- "YTHDF2_peak_over_m6A_with_random_mean"
  plot_reference_locus_with_random1(txdb, peakfile, centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabel, centerlabels,
                                    smo=FALSE, CLIP_reads=FALSE, extend_reads=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                    genome="hg19", rm.outlier=F)

  peakfile <- "ZNF121_ABCDE_peaks.bed"
  peaklabel <- "ZNF121_ABCDE"
  annotate_peaks(peakfile, peaklabel, gtffile, genome, fiveP, threeP, simple, RNA=T)
  op <- "ZNF121_peak_over_m6A_with_random_mean"
  plot_reference_locus_with_random1(txdb, peakfile, centerfiles, ext=c(-200,200), hl=c(-50,50), peaklabel, centerlabels,
                                    smo=FALSE, CLIP_reads=FALSE, extend_reads=0, norm=T, binsize=10, refPoint="center", Xlab="Center",
                                    inputfiles=NULL,inputlabels=NULL, stranded=TRUE, heatmap=FALSE, scale=FALSE, outPrefix=op,
                                    genome="hg19", rm.outlier=F)
}
