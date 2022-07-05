library(factR)
library(chromstaR)
library(ggpval)
library(Repitools)


source("C:/GREENBLATT/Rscripts/GenomicPlot/R/GenomicPlot.R")
#script.dir <- dirname(sys.frame(1)$ofile)
#source(file.path(script.dir, "GenomicPlot.R"))
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
      tx_longest <- transcripts(txdb, filter=list(tx_name=tl_longest_tx$tx_name), use.name=F)
      tss_tx_longest <- resize(tx_longest, width=1, fix="start", use.names=FALSE, ignore.strand=FALSE)
      seq_tss_tx_longest <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, tss_tx_longest)
      seq_df_longest <- as.data.frame(seq_tss_tx_longest)
      A_seq_longest <- rownames(seq_df_longest)[seq_df_longest$x == "A"]
      A_tss_tx_longest <- tss_tx_longest[tss_tx_longest$tx_id %in% A_seq_longest]
      return(A_tss_tx_longest)
   }else{
      tl <- transcriptLengths(txdb, with.utr5_len = T, with.cds_len = T, with.utr3_len = T)
      tl_protein_coding <- tl[tl$cds_len > 0, ]
      tx_protein_coding <- transcripts(txdb, filter=list(tx_name=tl_protein_coding$tx_name), use.name=F)
      tss_tx_protein_coding <- resize(tx_protein_coding, width=1, fix="start", use.names=FALSE, ignore.strand=FALSE)
      seq_tss_tx_protein_coding <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19, tss_tx_protein_coding)
      seq_df_protein_coding <- as.data.frame(seq_tss_tx_protein_coding)
      A_seq_protein_coding <- rownames(seq_df_protein_coding)[seq_df_protein_coding$x == "A"]
      A_tss_tx_protein_coding <- tss_tx_protein_coding[tss_tx_protein_coding$tx_id %in% A_seq_protein_coding]
      return(A_tss_tx_protein_coding)
   }
}

count_overlap_at_TSS <- function(queryfiles, txdb, TSSflank=0, CLIP_reads=FALSE, longest=TRUE){
   querylabels <- names(queryfiles)
   names(querylabels) <- queryfiles

   input_list <- handle_input(inputFiles=queryfiles, CLIP_reads=CLIP_reads, fix_width=0, useScore=F, outRle=F, useSizeFactor=F, genome="hg19")

   A_tss <- get_grWithA_at_TSS(txdb, longest=longest)
   gr <- flank(A_tss, width = TSSflank, both=T)

   count_df <- NULL

   for(queryfile in queryfiles){
      #queryfile <- queryfiles[1]
      queryRegions <- input_list[[queryfile]]$query
      overlap_count <- countOverlaps(gr, queryRegions, type="any")

      count_df <- cbind(count_df, overlap_count)

   }

   colnames(count_df) <- querylabels
   rownames(count_df) <- get_genomicCoordinates(gr)

   return(count_df)
}

count_overlap_in_features1 <- function(queryfiles, txdb, longest=T, CLIP_reads=T, TSSflank=2){
   querylabels <- names(queryfiles)
   names(querylabels) <- queryfiles

   input_list <- handle_input(inputFiles=queryfiles, CLIP_reads=CLIP_reads, fix_width=0, useScore=F, outRle=F, useSizeFactor=F, genome="hg19")

   feature_count_list <- list()

   for(featureName in c("utr5", "cds", "utr3", "transcript", "TSN")){
      #featureName <- "utr5"
      if(featureName == "TSN"){
         ## process transcription start sites with 'A'
         A_tss <- get_grWithA_at_TSS(txdb, longest=longest)
         feature_grl <- flank(A_tss, width=TSSflank, both=T)
      }else{
         alist <- gtf_to_bed_longest_tx(txdb, featureName, longest=longest)
         original_grl <- alist$GRangesList
         if(featureName %in% c("transcript")){
            trimmed <- trimTranscripts(original_grl, start=TSSflank, end=0)
            feature_grl <- trimmed
         }else{
            feature_grl <- original_grl
         }
      }

      featureCount_df <- NULL

      for(queryfile in queryfiles){
         #queryfile <- queryfiles[1]
         #print(queryfile)
         queryRegions <- input_list[[queryfile]]$query
         overlap_count <- countOverlaps(feature_grl, queryRegions, type="any")

         featureCount_df <- cbind(featureCount_df, overlap_count)
      }

      colnames(featureCount_df) <- querylabels

      feature_count_list[[featureName]] <- as.data.frame(featureCount_df)
   }

   sizes <- sapply(input_list, function(x)x$size)
   return(list("count"=feature_count_list, "size"=sizes))
}

run_DESeq2 <- function(data_table, s2c, project, contrasts, plot=FALSE, useSizeFactor=TRUE){

   ddsMatrix <- DESeqDataSetFromMatrix(round(data_table),
                                       colData = s2c,
                                       design = ~ groups)
   class(ddsMatrix)

   if(!useSizeFactor){
      sizeFactors(ddsMatrix) <- rep(1, ncol(data_table))
   }

   dds <- DESeq(ddsMatrix, fitType = "local")

   print(paste(project, "size factors"))
   print(sizeFactors(dds))

   X <- try(rld <- vst(dds))
   if(class(X) != "try-error" && plot){
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

      if(plot){
         pdf(paste(treatSF, "MA_plot.pdf", sep="_"))
         plot(log2(res$baseMean), res$log2FoldChange, ylim=c(-5,5))
         abline(h=0)
         dev.off()
      }
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

plot_feature_overlap_count <- function(queryfiles, peakfile=NULL, txdb, sampleNames, longest=T, CLIP_reads=T, flank=0, norm=T, subsets=NULL, input_type="reads"){

   peaklabel <- ifelse(is.null(peakfile), "genomic_features", names(peakfile))
   input_type <- paste0(input_type,"_in_", peaklabel)

   feature_count_list <- count_overlap_in_features(queryfiles, peakfile, txdb, longest=longest, CLIP_reads=CLIP_reads, TSSflank=flank)
   if(norm){
      rpm <- feature_count_list[["inputsize"]]/median(feature_count_list[["inputsize"]])
   }else{
      rpm <- rep(1, length(feature_count_list[["inputsize"]]))
   }
   names(rpm) <- names(queryfiles)

   print("plot counts in features")
   pdf(paste("overlap_count_of_", input_type, ".pdf",sep=""), height=8, width=12)
   for(featureName in names(feature_count_list$count)){
      #featureName <- "transcript"
      if(!is.null(peakfile)){
         count_table_feature <- feature_count_list[["count"]][[featureName]][[peakfile]]
      }else{
         count_table_feature <- feature_count_list[["count"]][[featureName]]
      }

      if(nrow(count_table_feature) == 0) next
      count_table_feature <- count_table_feature[apply(count_table_feature,1,sum)>0,]
      count_table_feature <- t(t(count_table_feature)/rpm)

      sum_table_feature <- sapply(sampleNames, function(sampleN){
         sub_df <- as.data.frame(count_table_feature[, grepl(sampleN, colnames(count_table_feature))])
         asum <- apply(sub_df, 1, sum)
      })
      sum_table_feature <- as.data.frame(log2(sum_table_feature+1))

      sum_df <- tidyr::pivot_longer(sum_table_feature, cols=colnames(sum_table_feature), names_to="Group", values_to="Count") %>%
         mutate(Group=as.factor(Group))
      median_df <- aggregate(sum_df$Count, list(sum_df$Group), median)
      colnames(median_df) <- c("Group", "median")

      if(length(levels(sum_df$Group)) < 2) next

      p <- ggplot(sum_df, aes(x=Group, y=Count, fill=Group)) +
         geom_boxplot(notch=FALSE) +
         theme_classic() + ggtitle(paste(featureName, "n =", nrow(count_table_feature))) +
         theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=14)) +
         theme(legend.position = "none") +
         labs(y=expression(paste(log[2], " (Count)"))) +
         #scale_x_discrete(limits=c("GFP", featureName)) + # labels=c("TSN" = expression(paste(m^6,"Am")), featureName = expression(paste("5'-UTR ", m^6, "A")))) +
         theme(axis.text.x = element_text(face="bold", size=14, color="black")) +
         stat_summary(fun=mean, geom="point", shape=23, size=4, fill="black") +
         geom_signif(comparisons = combn(levels(sum_df$Group), 2, simplify = F), test="wilcox.test", map_signif_level=TRUE, step_increase = 0.1)
      #geom_text(data = _df, aes(x=feature, y = min(out_df$log2FoldChange)*1.1, label = paste("n=",count,sep="")))

      d <- ggplot(sum_df, aes(x=Count, color=Group)) +
         geom_density() +
         labs(x=expression(paste(log[2], " (Count)"))) +
         geom_vline(data=median_df, aes(xintercept=median, color=Group), linetype="dotted", size=1) +
         theme_classic()

      outp <- plot_grid(p, d, ncol = 2, rel_widths = c(1,1))
      print(outp)
   }
   dev.off()


   if(!is.null(subsets)){
      print("plot lfc for feautures")
      pdf(paste("lfc_at_TSNflank",flank, "_and_transcript_for_",input_type, ".pdf",sep=""), height=8, width=12)
      for(featureName in names(feature_count_list$count)){
         lfc_list <- list()
         for (subject in names(subsets)){
            subset <- subsets[[subject]]
            s2c <- data.frame(samples=subset, groups=c(rep(c("GFP", subject), times=c(2,2))))
            contrasts <- data.frame(groups=c(rep("groups", 1)), treat=subject, control=c(rep("GFP",1)))

            project <- paste(input_type, flank, sep="")

            if(!is.null(peakfile)){
               data_table_feature <- feature_count_list[["count"]][[featureName]][[peakfile]][, subset]
            }else{
               data_table_feature <- feature_count_list[["count"]][[featureName]][, subset]
            }

            if(nrow(data_table_feature) == 0) next
            data_table_feature <- data_table_feature[apply(data_table_feature,1,sum)>0,]

            print(rpm[subset])
            data_table_feature <- round(t(t(data_table_feature)/rpm[subset]))
            print(head(data_table_feature))
            data_table_feature[is.na(data_table_feature)] <- 0

            res_feature <- NULL
            X <- try(
               res_feature <- run_DESeq2(data_table_feature, s2c, featureName, contrasts, plot=F, useSizeFactor = F)
            )
            if(class(X) != "try-error"){
               lfc_df <- data.frame(Group=rep(subject, length(res_feature[[1]]$log2FoldChange)), log2FoldChange=res_feature[[1]]$log2FoldChange)
               lfc_list[[subject]] <- lfc_df
            }
         }

         out_df <- as.data.frame(Reduce(rbind, lfc_list))
         if(nrow(out_df) == 0) next

         colnames(out_df) <- c("Group", "log2FoldChange")
         out_df <- mutate(out_df, Group=as.factor(Group))
         if(length(levels(out_df$Group)) < 2) next

         if(length(levels(out_df$Group)) > 1){
            median_df <- aggregate(out_df$log2FoldChange, list(out_df$Group), median)
            colnames(median_df) <- c("Group", "median")

            count_df <- aggregate(out_df$log2FoldChange, list(out_df$Group), length)
            colnames(count_df) <- c("Group", "count")
            p <- ggplot(out_df, aes(x=Group, y=log2FoldChange, fill=Group)) +
               geom_boxplot(notch=FALSE) +
               theme_classic() + ggtitle(featureName) +
               theme(axis.title.x = element_blank(), axis.title.y = element_text(face="bold", size=14)) +
               theme(legend.position = "none") +
               labs(y=expression(paste(log[2], " (treat/control signal intensity)"))) +
               scale_x_discrete(limits=names(subsets)) + # labels=c("TSN" = expression(paste(m^6,"Am")), GroupName = expression(paste("5'-UTR ", m^6, "A")))) +
               theme(axis.text.x = element_text(face="bold", size=14, color="black")) +
               geom_text(data = count_df, aes(x=Group, y = min(out_df$log2FoldChange)*1.1, label = paste("n=",count,sep=""))) +
               geom_signif(comparisons = combn(levels(out_df$Group), 2, simplify = F), test="t.test", map_signif_level=TRUE, step_increase = 0.1)

            #p1 <- add_pval(p, pairs = list(c(1, 2)), test='t.test')

            d <- ggplot(out_df, aes(x=log2FoldChange, color=Group)) +
               geom_density() +
               geom_vline(data=median_df, aes(xintercept=median, color=Group), linetype="dotted", size=1) +
               theme_classic()

            outp <- plot_grid(p, d, ncol = 2, rel_widths = c(1,1))
            print(outp)
         }
      }
      dev.off()
   }
}

count_overlap_in_features <- function(queryfiles, peakfile=NULL, txdb, longest=T, fix_width=0, CLIP_reads=T, TSSflank=2){
   querylabels <- names(queryfiles)
   names(querylabels) <- queryfiles
   input_list <- handle_input(inputFiles=queryfiles, CLIP_reads=CLIP_reads, fix_width=fix_width, useScore=F, outRle=F, useSizeFactor=F, genome="hg19")

   if(!is.null(peakfile)){
      peaklabel <- names(peakfile)
      names(peaklabel) <- peakfile
      peak_list <- handle_input(inputFiles=peakfile, CLIP_reads=F, fix_width=fix_width, useScore=F, outRle=F, useSizeFactor=F, genome="hg19")
   }

   feature_count_list <- list()

   for(featureName in c("utr5", "cds", "utr3", "transcript", "TSN", "unrestricted")){
      #featureName <- "utr5"
      if(featureName == "TSN"){
         ## process transcription start sites with 'A'
         A_tss <- get_grWithA_at_TSS(txdb, longest=longest)
         feature_grl <- flank(A_tss, width=TSSflank, both=T)
      }else if(featureName == "unrestricted"){
         feature_grl <- ifelse(!is.null(peakfile), peak_list[[peakfile]]$query, NULL)
      }else{
         alist <- gtf_to_bed_longest_tx(txdb, featureName, longest=longest)
         original_grl <- alist$GRangesList
         if(featureName %in% c("transcript")){
            trimmed <- trimTranscripts(original_grl, start=TSSflank, end=0)
            feature_grl <- trimmed
         }else{
            feature_grl <- original_grl
         }
      }

      if(!is.null(peakfile)){
         peakRegions <- peak_list[[peakfile]]$query
         if(unique(runValue(strand(peakRegions))) %in% c("*", ".", "")){
            feature_grl <- filter_by_overlaps(peakRegions, feature_grl) ## feature_grl becomes peaks overlapping with feature
         }else{
            feature_grl <- filter_by_overlaps_stranded(peakRegions, feature_grl)
         }
         #feature_grl <- resize(feature_grl, fix="center", width=(TSSflank*2+1))
      }

      if(!is.null(feature_grl)){
         featureCount_df <- sapply(queryfiles, function(queryfile){
            #queryfile <- queryfiles[1]
            #print(queryfile)
            queryRegions <- input_list[[queryfile]]$query
            overlap_count <- countOverlaps(feature_grl, queryRegions, type="any")
         })

         if(!is.null(peakfile)){
            feature_count_list[[featureName]][[peakfile]] <- as.data.frame(featureCount_df)
         }else{
            feature_count_list[[featureName]] <- as.data.frame(featureCount_df)
         }
      }

   }

   bamsizes <- sapply(input_list, function(x)x$size)

   if(!is.null(peakfile)){
      peaksizes <- sapply(peak_list, function(x)x$size)
      return(list("count"=feature_count_list, "inputsize"=bamsizes, "peaksize"=peaksizes))
   }else{
      return(list("count"=feature_count_list, "inputsize"=bamsizes))
   }
}

filter_bed_by_genomicFeature <- function(bedfile, featureName, txdb=NULL, resizeFraction=NULL, longest=FALSE, featureGr=NULL, featureBed=NULL, maxgap=-1L){

   if(featureName %in% c("utr3", "utr5", "cds", "intron", "exon", "transcript", "gene")){
      feature <- gtf_to_bed_longest_tx(txdb, featureName, longest=longest)
      featureGr <- feature$GRanges
   }else if(is.null(featureGr)){
      if(!is.null(featureBed)){
         featureGr <- import.bed(featureBed)
      }else{
         stop("Both featureGr and featureBed are null!")
      }
   }

   if(!is.null(resizeFraction)){
      new_w <- as.integer(width(featureGr)*resizeFraction)
      featureGr <- resize(featureGr, width=new_w, fix="start", use.names=TRUE, ignore.strand=FALSE)
   }

   bedGr <- import.bed(bedfile)

   #filtered <- filter_by_overlaps_stranded(bedGr, featureGr, maxgap=maxgap)
   overlaps <- GenomicRanges::findOverlaps(bedGr, featureGr, maxgap=maxgap)

   print(paste("from", length(unique(overlaps@from)), "to", length(unique(overlaps@to))))

   return(list(bed=bedGr[unique(overlaps@from)], feature=featureGr[unique(overlaps@to)]))
}

## bedfile is the "A" containing miCLIP peaks file
derive_m6Am_sites <- function(bedfile, m6Afile, txdb, m6Amfile){
   ## get "A" containing peaks in 5' UTR
   utr5_filtered <- filter_bed_by_genomicFeature(bedfile, featureName="utr5", txdb=txdb, resizeFraction=0.25, longest=TRUE)
   m6A <- import.bed(m6Afile)
   ## remove peaks that are overlapping with m6A sites
   m6Am <- filter_by_nonoverlaps_stranded(utr5_filtered$bed, m6A)
   ## convert to 6-column bed format
   utr5_filtered_bed <- annoGR2DF(m6Am) %>%
      select(c(chr, start, end, name, score, strand)) %>%
      mutate(start=as.integer(start)-1) %>%
      mutate(chr=as.character(chr))
   print(dim(utr5_filtered_bed))
   write.table(utr5_filtered_bed, m6Amfile, col.names=F, row.names = F, sep="\t", quote=F)
}

