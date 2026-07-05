#THE FUNCTION RECEIVES TWO GENELISTS AND COMPUTES THE RATIO OF ENRICHMENT BETWEEN THE TWO

#genelist - genelist - a gene set of interest to analyze - can be a vector or a list 
#DBversion - 12, 14, HPO 2022-04-14 release, HPO 2024-04-26 Release. Recommended: 2024-04-26 Release
#curationLevel - confident or confident+tentative. confident: observed in HPO. confident+tentative: observed in HPO or DisGeNET (which includes also mouse and rat phenotypes). HPO 2022-04-14 release, HPO 2024-04-26 Release take only confident

#OPTIONAL INPUTS:
#background - a larger set of genes, must include all genes in genelist or empty
#enrORdep - the direction of the test, can be enr for enrichment (higher probability to see term in genelist), dep (lower probability) or both (a non directional test) 
#minGenes - minimum genes found per term. Terms with fewer genes will not be FDR-corrected. If the user chose
#outpath - a path for saving the results file
#locORser - running locally ("loc") or on the server ("ser") - affects on DBdir. If you run locally - set working directory to: /Users/<username>/Dropbox (Weizmann Institute)/
#comp_main - a string that describes the comparison
#comp1 - a string that describes the first enrichment test
#comp2 - a string that describes the second enrichment test#uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately
#FreqOfPheno - typical or all. typical: observed in >50% of patients. all: observed at any frequency, HPO 2022-04-14 release, HPO 2024-04-26 Release take only all.
#onlySkel - T or F. T - analyze only skeleton-related organs.

ORGANizer_compare <- function(genelist1,genelist2,outpath,background,curationLevel,FreqOfPheno,DBversion,minGenes,onlySkel,enrORdep, locORser, comp_main, comp1,comp2, uniq)
{
  source("/home/labs/davidgo/Collaboration/backup/Lab_tools/Gokhman_Enrich/final_scripts/ORGANizer_ORGANize.R")
  #library("qvalue")
  
  #DEFAULT PARAMETERS
  if (missing(background)) {background <- c()}
  if (missing(minGenes)) {minGenes <- 0}
  if (missing(onlySkel)) {onlySkel <- F}
  if (missing(curationLevel)) {curationLevel <- "confident"}
  if (missing(FreqOfPheno)) {FreqOfPheno <- "all"} #
  if (missing(DBversion)) {DBversion <- "HPO 2024-04-26 Release"}
  if (missing(enrORdep)) {enrORdep <- "both"}
  if(missing(outpath)) {outpath <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(locORser)) {locORser <- "ser"}
  if (missing(comp_main)) {comp_main <- "up_vs_down"}
  if (missing(comp1)) {comp1 <- "up_reg"}
  if (missing(comp2)) {comp2 <- "down_reg"}
  if (missing(uniq)) {uniq <- T}
  
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
  
  #ORGANize
  ORGANizer1 <- ORGANizer_Organize(genelist1,curationLevel=curationLevel,FreqOfPheno=FreqOfPheno,background=background,DBversion=DBversion,minGenes=0,onlySkel=onlySkel, outpath=outpath, enrORdep=enrORdep, locORser=locORser, comp=comp1, uniq=uniq)
  ORGANizer2 <- ORGANizer_Organize(genelist2,curationLevel=curationLevel,FreqOfPheno=FreqOfPheno,background=background,DBversion=DBversion,minGenes=0,onlySkel=onlySkel, outpath=outpath, enrORdep=enrORdep, locORser=locORser, comp=comp2, uniq=uniq)
  
  if (!is.data.frame(ORGANizer1)) {
    if (!is.data.frame(ORGANizer2)) {
      return(c(ORGANizer1, ORGANizer2))
    } else {
      return(ORGANizer1)
    }
  } else if (!is.data.frame(ORGANizer2)) {
    return(ORGANizer2)
  }
  #remove organs linked to fewer that minGenes
  rowmax <- apply(cbind(ORGANizer1[,"Observed"],ORGANizer2[,"Observed"]), 1, max)
  idx <- which(rowmax >= minGenes)
  ORGANizer1 <- ORGANizer1[idx,]
  ORGANizer2 <- ORGANizer2[idx,]
  
  colnames <- c("Type","Body Part","Ratio","Linked1","Total in DB1","Linked2","Total in DB2","Genes1","Genes2","P-value","FDR")
  outputCompare <- as.data.frame(matrix(,nrow=nrow(ORGANizer1),ncol=length(colnames)))
  colnames(outputCompare) <- colnames
  
  if (sum(ORGANizer1$Observed,na.rm=T) > 0 & sum(ORGANizer2$Observed,na.rm=T)) {
    #FISHER'S EXACT TEST
    for (ll in 1:nrow(ORGANizer1)) {
      outputCompare[ll,"Type"] <- ORGANizer1[ll,"Type"]
      outputCompare[ll,"Body Part"] <- ORGANizer1[ll,"Body Part"]
      outputCompare[ll,"Linked1"] <- ORGANizer1[ll,"Observed"]
      outputCompare[ll,"Total in DB1"] <- ORGANizer1[ll,"Genes from genelist in DB"]
      ll2=which(ORGANizer2['Body Part']==ORGANizer1[ll,'Body Part'] & ORGANizer2['Type']==ORGANizer1[ll,'Type'])
      outputCompare[ll,"Linked2"] <- ORGANizer2[ll2,"Observed"]
      outputCompare[ll,"Total in DB2"] <- ORGANizer2[ll2,"Genes from genelist in DB"]   
      outputCompare[ll,"Ratio"] <- (outputCompare[ll,"Linked1"]/ORGANizer1[ll,"Genes from genelist in DB"])/(outputCompare[ll,"Linked2"]/ORGANizer2[ll2,"Genes from genelist in DB"])
      outputCompare[ll,"Genes1"] <- ORGANizer1[ll,"Genes"]
      outputCompare[ll,"Genes2"] <- ORGANizer2[ll2,"Genes"]
      
      if (!is.na(outputCompare[ll,"Ratio"])) {
        alternative <- "two.sided"
        testMat <- rbind(c(outputCompare[ll,"Linked1"],outputCompare[ll,"Linked2"]),c(outputCompare[ll,"Total in DB1"]-outputCompare[ll,"Linked1"],outputCompare[ll,"Total in DB2"]-outputCompare[ll,"Linked2"]))
        a<-fisher.test(testMat,alternative=alternative)
        outputCompare[ll,"P-value"] <- a$p.value
      } else {
        outputCompare[ll,"P-value"] <- NA
      }

    }
    
    #FDR
    idx <- which(ORGANizer1[,"Type"] == "Organ" & (outputCompare[,"Linked1"] >= minGenes | outputCompare[,"Linked2"] >= minGenes))
    if (length(idx) > 0) {outputCompare[idx,"FDR"] <- p.adjust(outputCompare[idx,"P-value"], method = "BH")}
    idx <- which(ORGANizer1[,"Type"] == "System" & (outputCompare[,"Linked1"] >= minGenes | outputCompare[,"Linked2"] >= minGenes))
    if (length(idx) > 0) {outputCompare[idx,"FDR"] <- p.adjust(outputCompare[idx,"P-value"], method = "BH")}
    idx <- which(ORGANizer1[,"Type"] == "Region" & (outputCompare[,"Linked1"] >= minGenes | outputCompare[,"Linked2"] >= minGenes))
    if (length(idx) > 0) {outputCompare[idx,"FDR"] <- p.adjust(outputCompare[idx,"P-value"], method = "BH")}
    idx <- which(ORGANizer1[,"Type"] == "Germ layer" & (outputCompare[,"Linked1"] >= minGenes | outputCompare[,"Linked2"] >= minGenes))
    if (length(idx) > 0) {outputCompare[idx,"FDR"] <- p.adjust(outputCompare[idx,"P-value"], method = "BH")}
  }
  outputCompare=outputCompare[order(outputCompare$FDR),]
  if (substr(outpath,nchar(outpath)-1,nchar(outpath)) != "/") {outpath <- paste(outpath,"/",sep="")}
  write.table(outputCompare, 
              file=paste0(outpath,"ORGANizer_compare_",comp_main,".txt"), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(outputCompare)
}