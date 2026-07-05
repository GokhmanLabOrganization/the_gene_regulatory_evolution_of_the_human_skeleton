#THE FUNCTION RECEIVES TWO genelistS AND COMPUTES THE RATIO OF ENRICHMENT BETWEEN THE TWO

# genelist1 and genelist2 - genes to analyze (enter gene symbols)
#OPTIONAL INPUTS:
#outpath - a path for saving the results file
# background - a larger set of genes, must include all genes in genelist or empty
#minGenes - minimum genes found per term. Terms with fewer genes will not be FDR-corrected.
#TYPofALL - typical or all. typical: observed in >50% of patients. all: observed at any frequency. TYP option is only for 14
# HPO_dbVer - The HPO DB version. Options: "HPO 2024-04-26 Release", "HPO 2022-04-14 release", 14
# bodyparts - a list of bodyparts. When the user supplies a list of bodyparts, the function will only analyze phenotypes associated with these body parts. Any body part that appears in Gene ORGANizer can be used.
# systems - a list of systems. When the user supplies a list of systems, the function will only analyze phenotypes associated with these systems. Any system that appears in Gene ORGANizer can be used.
# directional - T or F. T: analyze only directional phenotypes (i.e., phenotypes that can be described on a scale)
#enrORdep - the direction of the test, can be enr for enrichment (higher probability to see term in genelist), dep (lower probability) or both (a non directional test) 
#locORser - running locally ("loc") or on the server ("ser") - affects on DBdir. If you run locally - set working directory to: /Users/<username>/Dropbox (Weizmann Institute)/
#comp_main - a string that describes the comparison
#comp1 - a string that describes the first enrichment test
#comp2 - a string that describes the second enrichment test
#uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately

HPO_compare <- function(genelist1, genelist2, outpath, background, minGenes, TYPorALL, HPO_dbVer, systems, bodyparts, directional, enrORdep,locORser,comp_main,comp1,comp2,uniq)
{
  source("/home/labs/davidgo/Collaboration/backup/Lab_tools/Gokhman_Enrich/final_scripts/HPO_enrich.R")

  #DEFAULT PARAMETERS
  if (missing(background)) {background <- c()}
  if (missing(minGenes)) {minGenes <- 5}
  if (missing(bodyparts)) {bodyparts <- "all"} else {bodyparts <- tolower(bodyparts)}
  if (missing(systems)) {systems = "all"}
  if (missing(TYPorALL)) {TYPorALL <- "ALL"} else {TYPorALL <- toupper(TYPorALL)}
  if (missing(HPO_dbVer)) {HPO_dbVer <- "HPO 2024-04-26 Release"}
  if (missing(directional)) {directional <- F}
  if (missing(outpath)) {outpath <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(enrORdep)) {enrORdep <- "both"}
  if (missing(locORser)) {locORser <- "ser"}
  if (missing(comp_main)) {comp_main <- "up_vs_down"}
  if (missing(comp1)) {comp1 <- "up_reg"}
  if (missing(comp2)) {comp2 <- "down_reg"}
  if (missing(uniq)) {uniq <- T}
  
  #HPO enrich
  HPO1 <- HPO_enrich(genelist=genelist1, outpath=outpath, background=background, minGenes=minGenes, TYPorALL=TYPorALL, HPO_dbVer=HPO_dbVer, systems=systems, bodyparts=bodyparts, directional=directional, enrORdep=enrORdep, comp=comp1, uniq=uniq)
  HPO2 <- HPO_enrich(genelist=genelist2, outpath=outpath, background=background, minGenes=minGenes, TYPorALL=TYPorALL, HPO_dbVer=HPO_dbVer, systems=systems, bodyparts=bodyparts, directional=directional, enrORdep=enrORdep, comp=comp2, uniq=uniq)
  
  # output validation
  if (!is.data.frame(HPO1)) {
    if (!is.data.frame(HPO2)) {
      return(c(HPO1, HPO2))
    } else {
      return(HPO1)
    }
  } else if (!is.data.frame(HPO2)) {
    return(HPO2)
  }
  #FISHER'S EXACT TEST
  colnames <- c("Pheno","HPO_ID","Ratio","Linked1","Total in DB1","Linked2","Total in DB2","Genes1","Genes2","P-value","FDR")
  outputCompare <- as.data.frame(matrix(,nrow=nrow(HPO1),ncol=length(colnames)))
  colnames(outputCompare) <- colnames
  for (ll in 1:nrow(HPO1)) {
    outputCompare[ll,"Pheno"] <- HPO1[ll,"Pheno"]
    outputCompare[ll,"HPO_ID"] <- HPO1[ll,"HPO_ID"]
    outputCompare[ll,"Linked1"] <- HPO1[ll,"Observed (n_entered)"]
    outputCompare[ll,"Total in DB1"] <- HPO1[ll,"N_entered"]
    ll2=which(HPO2$Pheno==HPO1[ll,"Pheno"])
    
    outputCompare[ll,"Linked2"] <- HPO2[ll2,"Observed (n_entered)"]
    outputCompare[ll,"Total in DB2"] <- HPO2[ll2,"N_entered"]
    outputCompare[ll,"Ratio"] <- (outputCompare[ll,"Linked1"]/HPO1[ll,"N_entered"])/(outputCompare[ll,"Linked2"]/HPO2[ll2,"N_entered"])
    outputCompare[ll,"Genes1"] <- HPO1[ll,"Genes"]
    outputCompare[ll,"Genes2"] <- HPO2[ll2,"Genes"]
    
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
  idx <- which(outputCompare[,"Linked1"] >= minGenes | outputCompare[,"Linked2"] >= minGenes)
  if (length(idx) > 0) {outputCompare[idx,"FDR"] <- p.adjust(outputCompare[idx,"P-value"], method = "BH")}
  outputCompare=outputCompare[order(outputCompare$FDR),]
  if (substr(outpath,nchar(outpath)-1,nchar(outpath)) != "/") {outpath <- paste(outpath,"/",sep="")}
  write.table(outputCompare, 
              file=paste(outpath,"HPO_StatMat_compare_",comp_main,"_minGenes",minGenes,".txt",sep=""), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(outputCompare)
}