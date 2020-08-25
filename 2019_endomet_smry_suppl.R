
orderSeqPlate <- function(x){
  x %>% 
    group_by(seqPlate) %>%
    dplyr::count(seqPlate) %>% 
    arrange(-n)
}

countSeqPlate <- function(x, sqPlates){
  x %>%
    group_by(bin_expr, seqPlate, ensembl_gene_id) %>%
    dplyr::count(seqPlate) %>%
    ungroup() %>%
    mutate(seqPlate = factor(seqPlate, levels=sqPlates)) %>%
    mutate(bin_expr = factor(bin_expr, c(1,0)))
}

#scripts for doing GSEA with Fisher's exact test
#R package: gSEAfish - Set enrichment analysis using Fisher's exact test
#also provides visualization tools
 
genesetEnrichmentFisher <- function(B, Targ, GS.GS, GS.names){ 
  #calculates the fisher exact t test for enrichment in a particular set of
  #target gene entrez IDs (T) among the number of background genes (B). GS.GS
  #takes the form "hGMT$genesets", where hGMT is a variable that stores the .gmt
  #file object downloaded from the msigdb website that was read in using the
  #GMT() function. GS.names takes the form "hGMT$geneset.names"
  
  #1. Initialize an empty dataframe
  t <- unlist(fisher.test(matrix(data = c(1,4,6,2), nrow = 2), conf.int = T))
  colNam <- c(names(t), "geneset.names", "MatchesBtwTargAndPathway", "NumbGenesInMsigdbPathway", "TargLen", "BackgroundLen")
  fisherDf <- data.frame(matrix(NA, nrow=1, ncol=13))
  names(fisherDf) <- colNam
  
  TargLen <- length(Targ)
  #2. Loop through each geneset and calculate pvalues
  for(i in 1:length(GS.GS)){
    pathwayEntrezIDs <- GS.GS[[i]] #these are characters
    M <- length(intersect(Targ, pathwayEntrezIDs)) #number of overlapping genes
    P <- length(pathwayEntrezIDs) #numb of genes in pathway
    conting.matrix <- matrix(c(M, TargLen-M, P-M, B-TargLen-P+M), nrow = 2, byrow = TRUE) #create contigency matrix
    fisherTest <- fisher.test(conting.matrix, conf.int = T)
    fisherDf <- rbind(fisherDf, c(unlist(fisherTest), GS.names[[i]], M, P, TargLen, B))
  }
  fisherDf = fisherDf[-1,]   #remove the row of NA's (in row 1)
  fisherDf$bhPval <- p.adjust(fisherDf$p.value, method = "BH")
  fisherDf$bonferroni <- p.adjust(fisherDf$p.value, method = "bonferroni")
  fisherDf$neg.log10.bhPval <- -log10(as.numeric(fisherDf$bhPval))
  fisherDf$dataset <-
    unlist(lapply(strsplit(fisherDf$geneset.names, "_"), function(x) x[1]))
  fisherDf <- fisherDf[order(fisherDf$bhPval),]
  #rank by bhPval pbal
  return(fisherDf)
}


forestplot <- function(d, title, xlab=expression(log[10] *"(Odds Ratio)"), ylab="Pathway", size=12){
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi)) +
    geom_pointrange(size = 1.25) +
    coord_flip() +
    #geom_hline(aes(x=0, yintercept = 1)) +
    ylab(xlab) +
    xlab(ylab) + #switch because of the coord_flip() above
    ggtitle(title) +
    scale_y_log10() +
    theme_classic() +
    geom_text(aes(label=numbMatches)) +
    guides(alpha=FALSE) + 
    theme(axis.text.x=element_text(size=size), axis.text.y = element_text(size = 7), legend.text = element_text(size = 12))
  return(p)
}

createForestDF <- function(data, title){
  d <- data.frame(x=data$geneset.names, y=as.numeric(data$`estimate.odds ratio`), ylo = as.numeric(data$conf.int1), yhi = as.numeric(data$conf.int2), pval=as.numeric(data$p.value), numbGeneInPway = data$NumbGenesInMsigdbPathway,
                  numbMatches=data$MatchesBtwTargAndPathway)
  d$x <-factor(d$x, levels=d[order(-d$y), "x"])
  return(forestplot(d, title, xlab=expression("Odds Ratio("*log[10] *" scale)")))
}

gseWrapper <- function(B, Targ, genesets, geneset.names, qval, OR){
  #dir can equal "up" or "down" for up or down regulated, respectivtly
  res_gse <- genesetEnrichmentFisher(B = B, Targ = Targ, GS.GS = genesets, GS.names=geneset.names)
  res_gse_filt <- res_gse[as.numeric(res_gse$bhPval) < qval & as.numeric(res_gse$`estimate.odds ratio`) > OR & as.numeric(res_gse$conf.int1) > OR,]
  return(res_gse_filt)
}


readHmKgRctmepathways <- function(gsea_path){
  hallmark <- GSA.read.gmt(paste0(gsea_path,"h.all.v5.1.entrez.gmt")) #50 genesets
  kegg <- GSA.read.gmt(paste0(gsea_path,"c2.cp.kegg.v5.1.entrez.gmt")) #185 genesets
  reactome <- GSA.read.gmt(paste0(gsea_path,"c2.cp.reactome.v5.1.entrez.gmt")) #673 genesets
  hallmarkKeggReactome <- list(genesets = c(hallmark$genesets, kegg$genesets, reactome$genesets), geneset.names =
                                 c(hallmark$geneset.names, kegg$geneset.names, reactome$geneset.names), geneset.descriptions =
                                 c(hallmark$geneset.descriptions, kegg$geneset.descriptions, reactome$geneset.descriptions))
  return(hallmarkKeggReactome)
}

ensgToEntrez2 <- function(ensg, geneMap){
  g <- geneMap %>% 
    dplyr::select(-ucsc) %>% 
    distinct %>% 
    filter(ensembl_gene_id %in% ensg)
  return(g)
}


geneSymbToEnsg <- function(geneSymb){
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  geneEntrez <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene"),
                      filters = "hgnc_symbol",
                      values = geneSymb,
                      mart = ensembl)
  return(geneEntrez)
}






twoSdBelowMean <- function(mn, sd){
  twosdblw <- mn - 2*sd
  return(twosdblw)
}

twoSdAboveMean <- function(mn, sd){
  twosdabv <- mn + 2*sd
  return(twosdabv)
}

threeSdBelowMean <- function(mn, sd){
  threesdblw <- mn - 3*sd
  return(threesdblw)
}

threeSdAboveMean <- function(mn, sd){
  threesdabv <- mn + 3*sd
  return(threesdabv)
}
fourSdBelowMean <- function(mn, sd){
  fsdblw <- mn - 4*sd
  return(fsdblw)
}

fourSdAboveMean <- function(mn, sd){
  fsdabv <- mn + 4*sd
  return(fsdabv)
}

poundToKg <- function(lb){
  kg_per_lb = 2.2046226218
  return(lb/kg_per_lb)
}

inToCm <- function(inches){
  cm_per_in = 2.54
  return(cm_per_in*inches)
}




plotSurvAndReturnStat <- function(timeToEvent, alive0dead1, classif, surv.col, gnNm=""){
  fit <- survfit(Surv(timeToEvent, alive0dead1) ~ classif)
  fit.stat <- survdiff(Surv(timeToEvent, alive0dead1) ~ classif)
  print(fit.stat)
  #650x250
  p <- GGally::ggsurv(fit, back.white = T, surv.col = surv.col, size.est = 1.25) + 
    ylim(c(0,1)) + theme_classic() + 
    ggtitle(gnNm) +
    theme(aspect.ratio = 1, legend.title=element_blank(), 
          text = element_text(size=42, color = "black"), 
          axis.title=element_blank(),
          plot.title =element_text(size = 54, margin=margin(0,0,-5,0), face="italic")) +
  #add the number of patients
  annotate("text", x = 1, y = .35, label = paste0(names(fit.stat$n)[1], ", N = ",fit.stat$n[1])) +
  annotate("text", x = 1, y = .15, label = paste0(names(fit.stat$n)[2], ", N = ",fit.stat$n[2])) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

  print(p)
  return(fit.stat)
}

plotSurvAndReturnStat2 <- function(timeToEvent, alive0dead1, classif, surv.col, gnNm=""){
  fit <- survfit(Surv(timeToEvent, alive0dead1) ~ classif)
  fit.stat <- survdiff(Surv(timeToEvent, alive0dead1) ~ classif)
  print(fit.stat)
  #650x250
  p <- GGally::ggsurv(fit, back.white = T, surv.col = surv.col, size.est = 1.25) + 
    ylim(c(0,1)) + theme_classic() + 
    ggtitle(gnNm) +
    theme(aspect.ratio = 1, legend.title=element_blank(), 
          text = element_text(size=25, color = "black"), 
          axis.title=element_blank(),
          plot.title =element_text(face="italic")) +
    #add the number of patients
    annotate("text", x = 1, y = .35, label = paste0(names(fit.stat$n)[1], ", N = ",fit.stat$n[1])) +
    annotate("text", x = 1, y = .15, label = paste0(names(fit.stat$n)[2], ", N = ",fit.stat$n[2])) +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  print(p)
  return(list(fit.stat = fit.stat, fit = fit))
}

plotSurvAndReturnStat2_likrat <- function(timeToEvent, alive0dead1, classif, surv.col, gnNm=""){
  fit <- survfit(Surv(timeToEvent, alive0dead1) ~ classif)
  fit.stat <- survdiff(Surv(timeToEvent, alive0dead1) ~ classif)
  print(fit.stat)
  #650x250
  p <- GGally::ggsurv(fit, back.white = T, surv.col = surv.col, size.est = 1.25) + 
    ylim(c(0,1)) + theme_classic() + 
    ggtitle(gnNm) +
    theme(aspect.ratio = 1, legend.title=element_blank(), 
          text = element_text(size=25, color = "black"), 
          axis.title=element_blank(),
          plot.title =element_text(face="italic")) +
    #add the number of patients
    annotate("text", x = 1, y = .35, label = paste0(names(fit.stat$n)[1], ", N = ",fit.stat$n[1])) +
    annotate("text", x = 1, y = .15, label = paste0(names(fit.stat$n)[2], ", N = ",fit.stat$n[2])) +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  print(p)
  return(list(fit.stat = fit.stat, fit = fit))
}


returnDeGenesLimma <- function(expr, design, colOfInt, pval=1, logFC=0, retAll=F){
  fit_t <- lmFit(object = assay(expr), design = design)
  ebfit_t <- eBayes(fit_t)
  if(retAll == T){
    topTableRes_t <- topTable(ebfit_t, coef = colOfInt, number = Inf) 
    return(topTableRes_t)
  }
  
  topTableRes_t <- topTable(ebfit_t, coef = colOfInt, number = nrow(assay(expr)), adjust.method = "BH", p.value = pval, lfc = logFC) 
  return(topTableRes_t)
}

bindClassToColData <- function(cdat, obj){
  #if the columns exist, replace them
  #if the column does not exist, add it
  cno <- colnames(obj)
  for(col in 1:ncol(obj)){
    cdat$temp1 <- obj[,col] #add a temporary column from obj
    if(cno[col] %in% colnames(cdat)){
      #the column does exist, so replace it
      cdat[, cno[col]] <- NULL
    }
    colnames(cdat)[which(colnames(cdat) == "temp1")] <- cno[col]
  }
  return(cdat)
}


plotGeneHist2 <- function(exprNml, exprTum, isof, title, clrs){

  tidyDf <- as.data.frame(cbind(as.numeric(c(exprTum[isof,], exprNml[isof,])),
                                as.factor(c(rep("tumor",ncol(exprTum)), rep("normal",ncol(exprNml))))),
                          stringsAsFactors=FALSE)
  colnames(tidyDf) <- c("expr", "type")
  print(mean(exprNml[isof,]))
  print(twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  expr <- type <- ..density.. <- NULL # Setting the variables to NULL first
  p1 <- ggplot(tidyDf, aes(x=expr, color=as.factor(type),
                           fill=as.factor(type),
                           group=as.factor(type))) +
    theme_classic() +
    geom_histogram(data=subset(tidyDf,type == 2),
                   alpha=0.2, aes(y=..density..), binwidth = .2, 
                   color=clrs[2], fill=clrs[2]) +
    geom_rug(alpha=0.3, show.legend=FALSE) +
    scale_fill_manual(values=clrs) + 
    scale_color_manual(values=clrs) + 
    theme(axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          legend.position="none",
          text = element_text(size=35, color="black"),
          plot.title=element_text(face="italic")) + coord_fixed(ratio=15) + 
    xlim(c(-0.1, 10)) +
    ylim(c(0,.63))+
    ggtitle(title) + 
    stat_function(fun="dnorm", colour=clrs[1],
                  args=list(mean=mean(exprNml[isof,]),
                            sd=sd(exprNml[isof,])), linetype="dashed", size=1.5) +
    geom_vline(size=2, alpha=0.8, color=clrs[1], xintercept = twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
    print(p1)
}



plotGene_bothHist <- function(exprNml, exprTum, isof, title, clrs){

  tidyDf <- as.data.frame(cbind(as.numeric(c(exprTum[isof,], exprNml[isof,])),
                                as.factor(c(rep("tumor",ncol(exprTum)), rep("normal",ncol(exprNml))))),
                          stringsAsFactors=FALSE)
  colnames(tidyDf) <- c("expr", "type")
  print(mean(exprNml[isof,]))
  print(twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  expr <- type <- ..density.. <- NULL # Setting the variables to NULL first
  p1 <- ggplot(tidyDf, aes(x=expr, color=as.factor(type),
                           fill=as.factor(type),
                           group=as.factor(type))) +
    theme_classic() +
    geom_histogram(data=subset(tidyDf,type == 1),
                   alpha=0.2, aes(y=..density..), binwidth = .2, 
                   color=clrs[1], fill=clrs[1]) +
    geom_histogram(data=subset(tidyDf,type == 2),
                   alpha=0.2, aes(y=..density..), binwidth = .2, 
                   color=clrs[2], fill=clrs[2]) +
    geom_rug(alpha=0.3, show.legend=FALSE) +
    scale_fill_manual(values=clrs) + 
    scale_color_manual(values=clrs) + 
    theme(axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          legend.position="none",
          text = element_text(size=35, color="black"),
          plot.title=element_text(face="italic")) + coord_fixed(ratio=15) + 
    xlim(c(-0.1, 10)) +
    ggtitle(title) + 
    geom_vline(size=2, alpha=0.8, color=clrs[1], xintercept = twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  print(p1)
}


plotGene_oneHist <- function(exprNml, exprTum, isof, title, clrs, vert_line = TRUE){
  
  tidyDf <- as.data.frame(cbind(as.numeric(c(exprTum[isof,], exprNml[isof,])),
                                as.factor(c(rep("tumor",ncol(exprTum)), rep("normal",ncol(exprNml))))),
                          stringsAsFactors=FALSE)
  colnames(tidyDf) <- c("expr", "type")
  print(mean(exprNml[isof,]))
  print(twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  expr <- type <- ..density.. <- NULL # Setting the variables to NULL first
  p1 <- ggplot(tidyDf, aes(x=expr, color=as.factor(type),
                           fill=as.factor(type),
                           group=as.factor(type))) +
    theme_classic() +
    geom_histogram(data=subset(tidyDf,type == 1),
                   alpha=0.2, aes(y=..density..), binwidth = .2, 
                   color=clrs[1], fill=clrs[1]) +
     scale_fill_manual(values=clrs) + 
    scale_color_manual(values=clrs) + 
    theme(axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          legend.position="none",
          text = element_text(size=35, color="black"),
          plot.title=element_text(face="italic")) + coord_fixed(ratio=15) + 
    xlim(c(-0.1, 10)) +
    ggtitle(title) #+ 
  if(vert_line == TRUE){
    p1 <- p1 + geom_vline(size=2, alpha=0.8, color=clrs[1], xintercept = twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  } 
  print(p1)
}


plotGene_oneHistNmlCurve <- function(exprNml, exprTum, isof, title, clrs, vert_line = TRUE, histogram = TRUE){
  
  tidyDf <- as.data.frame(cbind(as.numeric(c(exprTum[isof,], exprNml[isof,])),
                                as.factor(c(rep("tumor",ncol(exprTum)), rep("normal",ncol(exprNml))))),
                          stringsAsFactors=FALSE)
  colnames(tidyDf) <- c("expr", "type")
  print(mean(exprNml[isof,]))
  print(twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  expr <- type <- ..density.. <- NULL # Setting the variables to NULL first
  p1 <- ggplot(tidyDf, aes(x=expr, color=as.factor(type),
                           fill=as.factor(type),
                           group=as.factor(type))) +
    theme_classic() +
    scale_fill_manual(values=clrs) + 
    scale_color_manual(values=clrs) + 
    theme(axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          legend.position="none",
          text = element_text(size=35, color="black"),
          plot.title=element_text(face="italic")) + coord_fixed(ratio=15) + 
    xlim(c(-0.1, 10)) +
    ggtitle(title) + 
    stat_function(fun="dnorm", colour=clrs[1],
                  args=list(mean=mean(exprNml[isof,]),
                            sd=sd(exprNml[isof,])), linetype="dashed", size=1.5)
    
  if(histogram == TRUE){
    p1 <- p1 + geom_histogram(data=subset(tidyDf,type == 1),
                              alpha=0.2, aes(y=..density..), binwidth = .2, 
                              color=clrs[1], fill=clrs[1])
  }
  if(vert_line == TRUE){
    p1 <- p1 + geom_vline(size=2, alpha=0.8, color=clrs[1], xintercept = twoSdBelowMean(mn=mean(exprNml[isof,]),sd = sd(exprNml[isof,])))
  }

  print(p1)
}




pvalFxn <- function(survDifObj){
  p.val <- 1 - pchisq(survDifObj$chisq, length(survDifObj$n) - 1)
  return(p.val)
}



sig2star = function(s, breaks=c(-Inf, 0.0001,0.001,0.01,1), labels=c("***","**","*","")) {
  r <- s
  r[] <- as.character(cut(s, breaks=breaks, labels)) 
  r[is.na(r)] <- ""
  return(r)
}



#' Cluster and Reorder Rows and Columns
#'
#' Clusters the rows and columns of a numeric matrix using hclust
#'
#' @param matr matrix to be clustered
#' @param dimen Indicates the dimension that the clustering is performed in. One of either "row", "column", or "both". 
#'
#' @return returns a matrix of clustered values by either row, column, or both, 

clusterAndReorder <- function(matr, dimen="both"){
  #dimen can be one of row, column or both
  if(dimen == "row"){
    hcl <- hclust(dist(matr))
    matr.new <- matr[hcl$order,]
    return(matr.new)
  } else if(dimen == "column"){
    hcl <- hclust(dist(t(matr)))
    matr.new <- matr[,hcl$order]
    return(matr.new)
  } else if(dimen == "both"){
    hcl.r <- hclust(dist(matr))
    hcl.c <- hclust(dist(t(matr)))
    matr.new <- matr[hcl.r$order, hcl.c$order]
    return(matr.new)
  }
}



set2 <- RColorBrewer::brewer.pal(n = 8, "Set2")
bw <- RColorBrewer::brewer.pal(n = 9, "Greys")

ggBar <- list(
  theme_classic(),
  theme(legend.position="top",
        axis.ticks = element_blank(),
        axis.text = element_text(size=18, color="black"),
        axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5),
        axis.title = element_text(size=22)),
  scale_y_continuous(expand = c(0,0)))

ggHist <- list(
  theme_classic(),
  theme(legend.position="top",
        axis.text = element_text(size=24)),
  scale_y_continuous(expand = c(0,0)))


ggSeqPlate = list(
  geom_bar(stat="identity"),
  facet_grid(.~ensembl_gene_id),
  theme_classic(), 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5),
        axis.text.y = element_text(size=16)),
  scale_fill_manual(values=set2[2:1]))


clean_dataset <- function(e){
  
  #returns adj normal and tumor dataframes
  print(dim(e))
  
  #remove ffpe samples
  e <- e[,!e$is_ffpe]
  print("Dimensions after FFPE removed:")
  print(dim(e)) 
  
  #add technical variables
  colData(e)$hosp <- substr(rownames(colData(e)), 6,7)
  colData(e)$seqPlate <- substr(rownames(colData(e)), 22,25)
  
  ###remove all metaData columns that are NA
  colData(e) <- colData(e)[,colSums(is.na(colData(e)))<nrow(colData(e))]
  
  ##TPM normalization, then log2(tpm+1) transform
  eu <- assay(e)
  eu <- scale(eu, center=FALSE, scale=colSums(eu)) * 1000000
  eu <- log(eu+1, 2)
  assay(e) <- eu
  
  #separate tumor and adjacent normal samples
  et <- e[,substr(colnames(e),14,15) == "01"] 
  en <- e[,substr(colnames(e),14,15) == "11"]
  print("Dimensions of tumor and adj nml dataframes:")
  print(dim(et)); print(dim(en))
  
  #remove duplicates
  et <- et[,!duplicated(substr(colnames(et),1,12))] 
  en <- en[,!duplicated(substr(colnames(en),1,12))] 
  print("Dimensions of tumor and adj nml dataframes (duplicates removed):")
  print(dim(et)); print(dim(en))

  #remove patients with NA values
  et <- et[,!apply(assay(et), 2, anyNA)] 
  en <- en[,!apply(assay(en), 2, anyNA)] 
  print("Dimensions of tumor and adj nml dataframes (duplicates removed + pts with NA values removed):")
  print(dim(et)); print(dim(en))
  
  #Filter out genes where <20% of patients have a non zero expression value
  et <- et[rowSums(assay(et)==0)<=ncol(et)*.80,] #select genes w/ less than 80% zeros
  en <- en[rowSums(assay(en)==0)<=ncol(en)*.80,] #select genes w/ less than 80% zeros
  print("Dimensions of tumor and adj nml dataframes (duplicates removed + pts with NA values removed + 20% filter):")
  print(dim(et)); print(dim(en))
  
  #get the genes in common between the two
  genes_in_common <- intersect(rownames(et), rownames(en))
  et <- et[genes_in_common,]
  en <- en[genes_in_common,]
  print("Dimensions of tumor and adj nml dataframes (duplicates removed + pts with NA values removed + 20% filter + intersection):")
  print(dim(et)); print(dim(en))
  
  return(list(et=et, en=en))
}


centroid_plot_elem <- function(){
  list(ggplot2::theme_classic(), 
         geom_point(size=3),
         geom_text(aes(label=Type)))
}

returnDeMirnaLimma <- function(expr, design, colOfInt, pval=1, logFC=0, retAll=F){
  fit_t <- lmFit(object = expr, design = design)
  ebfit_t <- eBayes(fit_t)
  if(retAll == T){
    topTableRes_t <- topTable(ebfit_t, coef = colOfInt, number = Inf) 
    return(topTableRes_t)
  }
  
  topTableRes_t <- topTable(ebfit_t, coef = colOfInt, number = nrow(assay(expr)), adjust.method = "BH", p.value = pval, lfc = logFC) 
  #topTableRes_t2 <-topTableRes_t[{topTableRes_t$logFC > logFC} | {topTableRes_t$logFC < logFC},] 
  return(topTableRes_t)
}

assay_to_df <- function(sum_expt, tum_type, gene){
  data.frame(assay(sum_expt[gene,])) %>% 
    t %>% 
    as_tibble() %>%
    mutate(tum_type=tum_type)
}   


returnClassif <- function(tum_matr, tum_res, gene){
  res <- {assay(tum_matr[gene,]) > {tum_res %>% filter(geneNm == gene) %>% pull(lowerBound)}} %>%
    ifelse(1,0) %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% rename(!!paste0("classif_", gene) := gene)
  return(res)
}

calc_days_last_fu <- function(days_to_death, days_to_last_follow_up){
  days_last_fu <- ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death)
  days_last_fu <- ifelse(days_last_fu < 1, 1, days_last_fu)
  days_last_fu <- ifelse(days_last_fu > 365.25 * 5, 365.25 * 5, days_last_fu)
  days_last_fu <- as.data.frame(days_last_fu)
  return(days_last_fu)
}

calc_vital_stat <- function(vital_stat, days_last_fu ){
  vital_stat <- ifelse(vital_stat == "alive", 0, 1)
  vital_stat <- ifelse(days_last_fu >= 365.25 * 5, 0, vital_stat)
  vital_stat <- as.data.frame(vital_stat)
  return(vital_stat)
}
