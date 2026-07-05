#THE FUNCTION RECEIVES A GENE LIST AND COMPARES ITS ENRICHMENT/DEPLETION COMPARED TO THE BACKGROUND (GENOMIC OR USER-SUPPLIED BACKGROUND LIST)

#genelist - a gene set of interest to analyze - can be a vector or a list 
#DBversion - 12, 14, HPO 2022-04-14 release, HPO 2024-04-26 Release. Recommended: 2024-04-26 Release
#curationLevel - confident or confident+tentative. confident: observed in HPO. confident+tentative: observed in HPO or DisGeNET (which includes also mouse and rat phenotypes). HPO 2022-04-14 release, HPO 2024-04-26 Release take only confident

#OPTIONAL INPUTS:
#background - a larger set of genes, must include all genes in genelist or empty
#enrORdep - the direction of the test, can be enr for enrichment (higher probability to see term in genelist), dep (lower probability) or both (a non directional test) 
#minGenes - minimum genes found per term. Terms with fewer genes will not be FDR-corrected. If the user chose
#outpath - a path for saving the results file
#locORser - running locally ("loc") or on the server ("ser") - affects on DBdir. If you run locally - set working directory to: /Users/<username>/Dropbox (Weizmann Institute)/
#comp - a string that describes our test – can be “active vs. all” for example. 
#uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately
#FreqOfPheno - typical or all. typical: observed in >50% of patients. all: observed at any frequency, HPO 2022-04-14 release, HPO 2024-04-26 Release take only all.
#onlySkel - T or F. T - analyze only skeleton-related organs.

ORGANizer_Organize <- function(genelist,background,outpath,curationLevel,FreqOfPheno,DBversion,minGenes,onlySkel,enrORdep,locORser,comp,uniq)
{
  if (missing(background)) {background <- c()}
  if (missing(minGenes)) {minGenes <- 0}
  if (missing(onlySkel)) {onlySkel <- F}
  if (missing(curationLevel)) {curationLevel <- "confident"}
  if (missing(FreqOfPheno)) {FreqOfPheno <- "all"} #
  if (missing(DBversion)) {DBversion <- "HPO 2024-04-26 Release"}
  if (missing(enrORdep)) {enrORdep <- "both"}
  if(missing(outpath)) {outpath <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(locORser)) {locORser <- "ser"}
  if (missing(comp)) {comp <- "diff_activity"}
  if (missing(uniq)) {uniq <- T}
  
  # Input validation - turns lists to vectors
  
  if (class(background)=="list"){
    background=unlist(background)
  }
  if (class(genelist)=="list"){
    genelist=unlist(genelist)
  }
  # Input validation - checks that each gene that appears in genelist appears in background at least as many times 
  
  if (length(background) > 0){
    count_genelist=as.data.frame(table(genelist))
    count_background=table(background)
    gene_validation<-function(x){
      gene=unname(x['genelist'])
      if (gene %in% names(count_background)){
        bg=as.integer(count_background[gene])
        gl=as.integer(x['Freq'])
        if (gl<=bg){
          return(T)
        } else{return(F)}
      } else{ return(F)}
    }
    valid_vec=apply(count_genelist,1,FUN=gene_validation)
    if (sum(valid_vec)<nrow(count_genelist)){
      invalid_genes=count_genelist[which(!valid_vec),'genelist']
      return(paste("Invalid input:",paste(invalid_genes,collapse = ", ")))
    }
  }
  
  if (uniq){
    genelist=unique(genelist)
    background=unique(background)
  }
  #version
  DBversion=as.character(DBversion)
  
  #curationLevel and freq of pheno
  if (tolower(curationLevel) == "confident") {
    if (tolower(FreqOfPheno) == "typical") {
      appendix <- "HPO_TYP"
    } else if (tolower(FreqOfPheno) == "all") {
      appendix <- "HPO_ALL"
    } else {
      print("wrong frequency of phenotype entered")
    }
  } else if (tolower(curationLevel) == "confident+tentative") {
    appendix <- "DisHPO_v"
  } else {
    print("wrong curation level entered")
  }
  
  #load ORGANizer DB #these files are just excel versions of the ORGANizer_World mat files
  library("readxl")
  if (locORser == "ser"){
    ORGANizer_World_path <- paste("/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Annotation DBs/Gene ORGANizer/",DBversion,"/ORGANizer_World_",appendix,"_",DBversion,".xlsx",sep="")
  } else{  ORGANizer_World_path <- paste("Gokhman lab general info/USEFUL_DATASETS/Annotation DBs/Gene ORGANizer/",DBversion,"/ORGANizer_World_",appendix,"_",DBversion,".xlsx",sep="")
}
  Organs <- read_excel(ORGANizer_World_path,sheet="Organs",col_names=T)
  Systems <- read_excel(ORGANizer_World_path,sheet="Systems",col_names=T)
  GermLayers <- read_excel(ORGANizer_World_path,sheet="GermLayers",col_names=T)
  Regions <- read_excel(ORGANizer_World_path,sheet="Regions",col_names=T)
  
  #Take only skeletal organs
  if (onlySkel) {
    #this option takes genes that are associated with the skeleton
    idx <- which(Systems[,"skeleton"] == 1)
    Organs <- Organs[idx,]
    Regions <- Regions[idx,]
    GermLayers <- GermLayers[idx,]
    Systems <- Systems[idx,]
    #this option takes only skeletal organs
    if (locORser == "ser"){
      Organs_skel <-read_excel("/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Annotation DBs/Gene ORGANizer/ORGANizer_skeletalBodyParts.xlsx",skip=0,col_names=T)
    } else{      Organs_skel <- read_excel("Gokhman lab general info/USEFUL_DATASETS/Annotation DBs/Gene ORGANizer/ORGANizer_Systems2BodyParts.xlsx",skip=0,col_names=T)
    }
    Organs_skel <- unlist(Organs_skel[which(Organs_skel[,2]==1),1])
    idx <- c(1,2,match(Organs_skel,colnames(Organs)))
    Organs <- Organs[,idx]
    #Systems <-Systems[,'skeleton']
  }
  
  
  #intersect background with ORGANizer DB
  if (length(background) > 0) {
    #intersectedBackGenes <- intersect(Organs$gene_symbol,unlist(background))
    intersectedBackGenes<-background[background %in% Organs$gene_symbol]
    idxBack <- match(intersectedBackGenes,Organs$gene_symbol)
    Organs <- Organs[idxBack,]
    Systems <- Systems[idxBack,]
    GermLayers <- GermLayers[idxBack,]
    Regions <- Regions[idxBack,]
  }
  
  #intersect genelist with ORGANizer DB
  uniqueGenes <- unique(unlist(genelist))
  #intersectedGenes <- intersect(Organs$gene_symbol,unlist(genelist))
  intersectedGenes <- genelist[genelist %in% Organs$gene_symbol]
  idxGenes <- match(intersectedGenes,Organs$gene_symbol)
  Organs_genelist <- Organs[idxGenes,]
  Systems_genelist <- Systems[idxGenes,]
  GermLayers_genelist <- GermLayers[idxGenes,]
  Regions_genelist <- Regions[idxGenes,]
  
  colnames <- c("Unique IDs entered","Genes from genelist in DB","Type","Body Part","Enrichment","Observed","Expected","Genes","P-value","FDR")
  outTable <- as.data.frame(matrix(NA,nrow=(ncol(Organs)-2+ncol(Systems)-2+ncol(Regions)-2+ncol(GermLayers)-2),ncol=length(colnames)))
  colnames(outTable) <- colnames
  
  outTable[,"Type"] <- c(rep("Organ",times=ncol(Organs)-2),rep("System",times=ncol(Systems)-2),rep("Region",times=ncol(Regions)-2),rep("Germ layer",times=ncol(GermLayers)-2))
  outTable[,"Body Part"] <- c(colnames(Organs)[3:ncol(Organs)],colnames(Systems)[3:ncol(Systems)],colnames(Regions)[3:ncol(Regions)],colnames(GermLayers)[3:ncol(GermLayers)])
  #go over body parts
  ll <- 1
  ll_strt <- ll
  if (length(idxGenes) > 0) {
    for (oo in 3:ncol(Organs)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(Organs_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(Organs_genelist[which(Organs_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(Organs_genelist[,oo])/nrow(Organs[,oo])) * sum(Organs[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_org <- matrix(c(outTable[ll,"Observed"], sum(Organs[,oo]) - outTable[ll,"Observed"],
                          nrow(Organs_genelist[,oo]) - outTable[ll,"Observed"],
                          ((nrow(Organs[,oo])-sum(Organs[,oo]))) - (nrow(Organs_genelist[,oo]) - outTable[ll,"Observed"])),
                        nrow=2, byrow=T)
      row_sums=rowSums(mat_org)
      col_sums=colSums(mat_org)
      tot_sum=sum(mat_org)
      exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                       row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
      bool_mat=exp_mat<mat_org
      bool_mat=bool_mat*1
      mat_org_corrected=mat_org+bool_mat
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_org_corrected, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_org_corrected, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_org_corrected, alternative = "t")$p.value
      }

      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
    
    ll_strt <- ll
    for (oo in 3:ncol(Systems)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(Systems_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(Systems_genelist[which(Systems_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(Systems_genelist[,oo])/nrow(Systems[,oo])) * sum(Systems[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_sys <- matrix(c(outTable[ll,"Observed"], sum(Systems[,oo]) - outTable[ll,"Observed"],
                          nrow(Systems_genelist[,oo]) - outTable[ll,"Observed"],
                          ((nrow(Systems[,oo])-sum(Systems[,oo]))) - (nrow(Systems_genelist[,oo]) - outTable[ll,"Observed"])),
                        nrow=2, byrow=T)
      row_sums=rowSums(mat_sys)
      col_sums=colSums(mat_sys)
      tot_sum=sum(mat_sys)
      exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                       row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
      bool_mat=exp_mat<mat_sys
      bool_mat=bool_mat*1
      mat_sys_corrected=mat_sys+bool_mat
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_sys_corrected, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_sys_corrected, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_sys_corrected, alternative = "t")$p.value
      }
      

      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
    
    ll_strt <- ll
    for (oo in 3:ncol(Regions)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(Regions_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(Regions_genelist[which(Regions_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(Regions_genelist[,oo])/nrow(Regions[,oo])) * sum(Regions[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_reg <- matrix(c(outTable[ll,"Observed"], sum(Regions[,oo]) - outTable[ll,"Observed"],
               nrow(Regions_genelist[,oo]) - outTable[ll,"Observed"],
               ((nrow(Regions[,oo])-sum(Regions[,oo]))) - (nrow(Regions_genelist[,oo]) - outTable[ll,"Observed"])),
             nrow=2, byrow=T)
      row_sums=rowSums(mat_reg)
      col_sums=colSums(mat_reg)
      tot_sum=sum(mat_reg)
      exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                       row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
      bool_mat=exp_mat<mat_reg
      bool_mat=bool_mat*1
      mat_reg_corrected=mat_reg+bool_mat
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_reg_corrected, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_reg_corrected, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_reg_corrected, alternative = "t")$p.value
      }
      
      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
    
    ll_strt <- ll
    for (oo in 3:ncol(GermLayers)) {
      outTable[ll,"Unique IDs entered"] <- length(uniqueGenes)
      outTable[ll,"Genes from genelist in DB"] <- nrow(Organs_genelist[,oo])
      outTable[ll,"Observed"] <- sum(GermLayers_genelist[,oo])
      outTable[ll,"Genes"] <- paste(unlist(GermLayers_genelist[which(GermLayers_genelist[,oo] == 1),"gene_symbol"]),collapse=",")
      outTable[ll,"Expected"] <- (nrow(GermLayers_genelist[,oo])/nrow(GermLayers[,oo])) * sum(GermLayers[,oo])
      outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
      
      mat_gel <- matrix(c(outTable[ll,"Observed"], sum(GermLayers[,oo]) - outTable[ll,"Observed"],
                          nrow(GermLayers_genelist[,oo]) - outTable[ll,"Observed"],
                          ((nrow(GermLayers[,oo])-sum(GermLayers[,oo]))) - (nrow(GermLayers_genelist[,oo]) - outTable[ll,"Observed"])),
                        nrow=2, byrow=T)
      row_sums=rowSums(mat_gel)
      col_sums=colSums(mat_gel)
      tot_sum=sum(mat_gel)
      exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                       row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
      bool_mat=exp_mat<mat_gel
      bool_mat=bool_mat*1
      mat_gel_corrected=mat_gel+bool_mat
      if (enrORdep == "enr") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_gel_corrected, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_gel_corrected, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        outTable[ll,"P-value"] <- stats::fisher.test(mat_gel_corrected, alternative = "t")$p.value
      }
      
      
      ll <- ll+1
    }
    #fdr
    pvals_temp <- outTable[ll_strt:(ll-1),"P-value"]
    idxminGenes <- which((outTable[ll_strt:(ll-1),"Enrichment"] > 1 & outTable[ll_strt:(ll-1),"Observed"] <= minGenes) | (outTable[ll_strt:(ll-1),"Enrichment"] < 1 & outTable[ll_strt:(ll-1),"Expected"] <= minGenes))
    pvals_temp[idxminGenes] <- NA
    outTable[ll_strt:(ll-1),"FDR"] <- p.adjust(pvals_temp,method="BH")
  }
  outTable=outTable[order(outTable$FDR),]
  #WRITE FILE
  if (substr(outpath,nchar(outpath)-1,nchar(outpath)) != "/") {outpath <- paste(outpath,"/",sep="")}
  write.table(outTable, file=paste0(outpath,"GeneORGANizer_ORGANize","_v",
                                    DBversion,"_","minGenes_",minGenes,"_",enrORdep,"_",
                                    FreqOfPheno,"_",comp,"_",onlySkel,".txt"), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(outTable)
}
