
## find the transcript that a utr3 belongs to
utr3_to_transcript <- function(utr3File, txdb){
   #utr3File <- centerfiles[1]

   transcriptFile <- gsub("\\.bed", "_TSS\\.bed", utr3File)

   longest_tx <- extract_longest_tx(txdb)
   utrInput <- handle_input(utr3File)
   utrGR <- utrInput$query

   w <- width(utrGR)

   feature <- transcripts(txdb, use.name=T)
   feature_longest <- feature[names(feature) %in% longest_tx$tx_name]

   feature_overlap <- subsetByOverlaps(feature_longest, utrGR, type="any", minoverlap=min(w))

   export.bed(feature_overlap, transcriptFile)

   return(transcriptFile)
}

extract_intron_VAST <- function(inputfile, length_filter=NULL){
   inputfile <- "C:/GREENBLATT/Nabeel/VAST-TOOLS/ZBTB7A/Intron_retention_regulated_events.tab"

   in_df <- read.delim(inputfile, header=TRUE)
   outfiles <- list()
   out_beds <- list()
   out_lens <- list()

   for(change in c("_up", "_down", "_unregulated")){
      #change <- "_unregulated"
      print(change)


      outfiles[[change]] <- gsub("\\.tab", paste0(change,"\\.bed"), inputfile)

      sub_df <- filter(in_df, grepl(change, GROUP))

      fullcoor <- lapply(as.list(sub_df$FullCO), function(x){unlist(strsplit(x, split=":", fixed=T))})
      fullcoor_strand <- unlist(lapply(fullcoor, function(x)x[3]))

      introncoor <- lapply(as.list(sub_df$COORD), function(x){unlist(strsplit(x, split=":|-", fixed=F))})
      introncoor_chr <- unlist(lapply(introncoor, function(x)x[1]))
      introncoor_start <- unlist(lapply(introncoor, function(x)x[2]))
      introncoor_end <- unlist(lapply(introncoor, function(x)x[3]))


      out_bed <- data.frame("chr"=introncoor_chr, "start"=as.integer(introncoor_start), "end"=as.integer(introncoor_end), "id"=sub_df$EVENT, "score"=sub_df$MV.dPsi._at_0.95, "strand"=fullcoor_strand)
      #write.table(out_bed, outputfile, row.names=F, col.names=F, sep="\t", quote=F)
      out_beds[[change]] <- out_bed
      out_lens[[change]] <- out_bed$end - out_bed$start + 1
   }

   ## test if up and down have different intron length distribution from the unregulated
   k <- ks.test(out_lens[["_unregulated"]], c(out_lens[["_up"]],out_lens[["_down"]]))

   print(k$p.value)
   print("Performing sub-sampling to match length distribution")
   sub <- sub_sample(out_lens[["_unregulated"]], c(out_lens[["_up"]],out_lens[["_down"]]), 0.9)
   k <- ks.test(sub, c(out_lens[["_up"]],out_lens[["_down"]]))
   print(k$p.value)
   print(summary(sub))

   plot(quantile(c(out_lens[["_unregulated"]]), probs = seq(0, 1, 0.01)), col="cyan")
   points(quantile(c(out_lens[["_up"]], out_lens[["_down"]]), probs = seq(0, 1, 0.01)), col="blue")
   points(quantile(sub, probs = seq(0, 1, 0.01)), col="red")

   qqplot(x=log(sub), y=log(c(out_lens[["_up"]], out_lens[["_down"]])))

   sub_index <- match(sub, out_lens[["_unregulated"]])

   out_beds[["_unregulated"]] <- data.frame(out_beds[["_unregulated"]])[sub_index,]

   ## filter final output using length_filter
   for(change in c("_up", "_down", "_unregulated")){
      if(is.integer(length_filter)){
         out_beds[[change]] <- data.frame(out_beds[[change]]) %>% filter((end-start+1) > length_filter)
      }
      write.table(out_beds[[change]], outfiles[[change]], row.names=F, col.names=F, sep="\t", quote=F)
   }

   unlist(outfiles)
}



get_grWithA_at_TSS <- function(txdb, longest=TRUE){

   if(longest){
      tl_longest_tx <- extract_longest_tx(txdb)
      tx_longest <- transcripts(txdb, filter=list(tx_name=tl_longest_tx$tx_name), use.name=T)
      tss_tx_longest <- resize(tx_longest, width=1, fix="start", use.names=TRUE, ignore.strand=FALSE)
      seq_tss_tx_longest <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, tss_tx_longest)
      seq_df_longest <- as.data.frame(seq_tss_tx_longest)
      A_seq_longest <- rownames(seq_df_longest)[seq_df_longest$x == "A"]
      A_tss_tx_longest <- tss_tx_longest[tss_tx_longest$tx_name %in% A_seq_longest]
      return(A_tss_tx_longest)
   }else{
      tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
      tl_protein_coding <- tl[tl$cds_len > 0, ]
      tx_protein_coding <- transcripts(txdb, filter=list(tx_name=tl_protein_coding$tx_name), use.name=T)
      tss_tx_protein_coding <- resize(tx_protein_coding, width=1, fix="start", use.names=TRUE, ignore.strand=FALSE)
      seq_tss_tx_protein_coding <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, tss_tx_protein_coding)
      seq_df_protein_coding <- as.data.frame(seq_tss_tx_protein_coding)
      A_seq_protein_coding <- rownames(seq_df_protein_coding)[seq_df_protein_coding$x == "A"]
      A_tss_tx_protein_coding <- tss_tx_protein_coding[tss_tx_protein_coding$tx_name %in% A_seq_protein_coding]
   }
}

count_overlap_at_TSS <- function(queryfiles, querylabels, txdb, ext=0, longest=TRUE){

   names(querylabels) <- queryfiles
   fileType <- "bed"
   if(grepl("\\.bam", queryfiles[1])){fileType <- "bam"}

   print(paste("The input type is", fileType, "file"))

   A_tss <- get_grWithA_at_TSS(txdb, longest=longest)
   gr <- flank(A_tss, width = ext, both=T)

   count_df <- NULL

   for(queryfile in queryfiles){
      #queryfile <- queryfiles[1]
      querylabel <- querylabels[queryfile]
      libsize <- NULL
      if(fileType == "bed"){
         queryRegions <- read_bed(queryfile)
         libsize <- countLines(queryfile)
      }else{

         ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
         bamGR <- readBamFileAsGRanges(queryfile, min.mapq = 10)
         l <- end(bamGR) - start(bamGR)
         summary(l)
         queryRegions =flank(bamGR, width=1, both=F, start=T, ignore.strand=FALSE)
         libsize <- length(queryRegions)

      }

      overlap_count <- countOverlaps(gr, queryRegions, type="any")

      count_df <- cbind(count_df, overlap_count)

   }

   colnames(count_df) <- querylabels

   return(count_df)
}

count_overlap_in_features <- function(queryfiles, querylabels, txdb, longest=T, CLIP_reads=T){

   names(querylabels) <- queryfiles
   fileType <- "bed"
   if(grepl("\\.bam", queryfiles[1])){fileType <- "bam"}

   print(paste("The input type is", fileType, "file"))

   input_list <- list()

   for(queryfile in queryfiles){
      #queryfile <- queryfiles[1]
      querylabel <- querylabels[queryfile]
      libsize <- NULL
      if(fileType == "bed"){
         bedGR <- read_bed(queryfile)
         queryRegions =flank(bedGR, width=1, both=F, start=T, ignore.strand=FALSE)
         queryRegions <- ifelse(CLIP_reads, queryRegions, bedGR)
         libsize <- countLines(queryfile)
      }else{
         ga <- readGAlignments(queryfile, use.names=T, param=ScanBamParam(mapqFilter=10))
         libsize <- length(ga)
         if(CLIP_reads){
            ## get the 5'-end -1 position of the reads, which is the crosslink sites for iCLIP reads
            queryRegions =flank(granges(ga), width=1, both=F, start=T, ignore.strand=FALSE)
         }else{
            queryRegions <- stack(grglist(ga))
         }

      }
      input_list[[querylabel]] <- queryRegions
   }

   feature_count_list <- list()

   for(featureName in c("utr5", "cds", "utr3", "exon")){
      print(featureName)
      #featureName <- "utr5"
      alist <- gtf_to_bed_longest_tx(txdb, featureName, longest=longest)
      original_grl <- alist$GRangesList
      if(featureName %in% c("utr5", "exon")){
         trimmed <- trimTranscripts(original_grl, start=1, end=0)
         feature_grl <- trimmed
      }else{
         feature_grl <- original_grl
      }

      featureCount_df <- NULL

      for(queryfile in queryfiles){
         #queryfile <- queryfiles[1]
         print(queryfile)
         querylabel <- querylabels[queryfile]
         queryRegions <- input_list[[querylabel]]

         overlap_count <- countOverlaps(feature_grl, queryRegions, type="any")

         featureCount_df <- cbind(featureCount_df, overlap_count)
      }

      colnames(featureCount_df) <- querylabels

      feature_count_list[[featureName]] <- featureCount_df
   }
   return(feature_count_list)
}

run_DESeq2 <- function(data_table, s2c, project, contrasts){

   ddsMatrix <- DESeqDataSetFromMatrix(round(data_table),
                                       colData = s2c,
                                       design = ~ groups)
   class(ddsMatrix)

   #sizeFactors(ddsMatrix) <- guide_stat_sum/guide_gomean
   dds <- DESeq(ddsMatrix, fitType = "local") #  , sfType = "iterate"
   #rnames <- resultsNames(dds)
   #rnames

   #normalized_counts <- counts(dds, normalized=TRUE)
   #lfcShrinked <- lfcShrink(dds=dds, coef=3, type="apeglm")
   #lfcShrinked <- lfcShrinked[order(lfcShrinked$padj), ]
   #plotMDS(normalized_counts)
   #pheatmap(cor(normalized_counts))
   print(paste(project, "size factors"))
   print(sizeFactors(dds))

   X <- try(rld <- vst(dds))
   if(class(X) != "try-error"){
      pdf(paste(project, "PCA_plot.pdf", sep="_"))
      print(plotPCA(rld, intgroup=c("groups")))
      dev.off()
   }

   res_list <- NULL

   for(i in 1:nrow(contrasts)){
      #i <- 1
      contrast <- as.vector(as.character(contrasts[i,]))
      print(paste("Processing ", paste(contrast, collapse=" ")))
      res <- results(dds, contrast=contrast, independentFiltering = TRUE, pAdjustMethod = "BH")
      treatSF <- paste(project, contrast[2], "vs", contrast[3], sep="_")

      summary(res)
      res <- na.omit(res)

      mcols(res)

      res <- res[order(-res$stat),]
      res_list[[treatSF]] <- res

      filename <- paste(treatSF, "DESeq2_results_table.tab", sep="_")
      write.table(res, filename, row.names=T, col.names=NA, quote=F, sep="\t")

      pdf(paste(treatSF, "MA_plot.pdf", sep="_"))
      #plotMA(res, ylim=c(-5,5))
      plot(log2(res$baseMean), res$log2FoldChange, ylim=c(-5,5))
      abline(h=0)
      dev.off()

   }
   return(res_list)
}


liftOverBed <- function(chainFile, queryfile){
   outFile <- gsub("\\.bed", "_liftOver\\.bed", queryfile)
   chain <- import.chain(chainFile)
   bedGr <- import.bed(queryfile)
   lifted <- stack(liftOver(bedGr, chain))
   mcols(lifted) <- mcols(lifted)[,2:3] ## remove the name column of numbers
   write_bed(lifted, outFile)
   return(outFile)
}

## v1: large vector to be sub-sampled
## v2: small vector to provide density distribution
## wpv: wilcox.test p-value to indicate similarity between sub-sample and v2
sub_sample <- function(v1, v2, pv){
   #v1 <- out_lens[["_unregulated"]]
   #v2 <- c(out_lens[["_up"]],out_lens[["_down"]])
   #pv <- 0.1


   n = 0
   vsub <- NULL
   ds <- density(v2)
   dv1 <- approx(x=ds$x, y=ds$y, xout=v1, rule=1)
   probs <- dv1$y
   probs[is.na(probs)] <- 0
   summary(probs)

   success <- FALSE
   wtest <- NULL

   ## find subset of vi that has the same distribution as v2
   while(!success){
      n <- n + 1
      vsub <- sample(v1, round(length(probs[probs>0])*0.99), replace=F, prob=probs)
      #vsub <- sample(v1, length(v2), replace=F, prob=probs)
      wtest <- wilcox.test(vsub, v2)
      ktest <- ks.test(vsub, v2)

      print(paste(n, wtest$p.value, ktest$p.value))
      if(wtest$p.value > pv){
         #if(lr < 1){
         success <- TRUE
         print("density distribution of v2")
         print(summary(v2))
         print("subsampling of v1 based on density distribution of v2")
         print(summary(vsub))
      }
      v1 <- vsub
      dv1 <- approx(x=ds$x, y=ds$y, xout=v1, rule=1)
      probs <- dv1$y
      probs[is.na(probs)] <- 0

   }

   ## produce the final subset of v1 that have the same length as v2
   vsub <- sample(v1, length(v2), replace=F, prob=probs)

   return(vsub)
}

