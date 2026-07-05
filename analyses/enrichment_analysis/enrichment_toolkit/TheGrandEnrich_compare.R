#THE FUNCTION RECEIVES TWO GENELISTS AND COMPUTES THE RATIO OF ENRICHMENT BETWEEN THE TWO

#db: GO_BP, GO_MF, GO_CC, KEGG, OMIM, REACTOME, GAD_Disease
#genelist1 and genelist2 - the genelists to compare enrichment between.
#OPTIONAL INPUTS:
#background - a larger set of genes, must include all genes in genelist or empty
#enrORdep - the direction of the test, can be enr for enrichment (higher probability to see term in genelist), dep (lower probability) or both (a non directional test) 
#minGenes - minimum genes found per term. Terms with fewer genes will not be FDR-corrected. If the user chose
#outdir - a path for saving the results file
#locORser - running locally ("loc") or on the server ("ser") - affects on DBdir. If you run locally - set working directory to: /Users/<username>/Dropbox (Weizmann Institute)/
#comp_main - a string that describes the comparison
#comp1 - a string that describes the first enrichment test
#comp2 - a string that describes the second enrichment test
#uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately


TheGrandEnrich_compare <- function(db,genelist1,genelist2,background,minGenes,outdir,enrORdep,locORser,comp_main,comp1,comp2,uniq)
{
  source("/home/labs/davidgo/Collaboration/backup/Lab_tools/Gokhman_Enrich/final_scripts/TheGrandEnrich_enrich.R")
  
  #DEFAULT PARAMETERS
  if (missing(db)) {db <- "GO_BP"}
  if (missing(background)) {background <- c()}
  if (missing(minGenes)) {minGenes <- 0}
  if (missing(outdir)) {outdir <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(enrORdep)) {enrORdep <- "both"}
  if (missing(locORser)) {locORser <- "ser"}
  if (missing(comp_main)) {comp_main <- "up_vs_down"}
  if (missing(comp1)) {comp1 <- "up_reg"}
  if (missing(comp2)) {comp2 <- "down_reg"}
  if (missing(uniq)) {uniq <- T}
  #ORGANize
  Enrich1 <- TheGrandEnrich_enrich(db=db,genelist1,background=background,minGenes=minGenes,outdir=outdir,enrORdep=enrORdep,locORser=locORser,comp=comp1,uniq=uniq)
  Enrich2 <- TheGrandEnrich_enrich(db=db,genelist2,background=background,minGenes=minGenes,outdir=outdir,enrORdep=enrORdep,locORser=locORser,comp=comp2,uniq=uniq)
  
  # output validation
  if (!is.data.frame(Enrich1)) {
    if (!is.data.frame(Enrich2)) {
      return(c(Enrich1, Enrich2))
    } else {
      return(Enrich1)
    }
  } else if (!is.data.frame(Enrich2)) {
    return(Enrich2)
  }
  
  #FISHER'S EXACT TEST
  colnames <- c("DB","Term","Ratio","Linked1","Total in DB1","Linked2","Total in DB2","Genes1","Genes2","P-value","FDR")
  outputCompare <- as.data.frame(matrix(,nrow=nrow(Enrich1),ncol=length(colnames)))
  colnames(outputCompare) <- colnames
  for (ll in 1:nrow(Enrich1)) {
    outputCompare[ll,"DB"] <- db
    outputCompare[ll,"Term"] <- Enrich1[ll,"Term"]
    outputCompare[ll,"Linked1"] <- Enrich1[ll,"Observed genes associated with term"]
    outputCompare[ll,"Total in DB1"] <- Enrich1[ll,"Genes in genelist that are in DB"]
    ll2=which(Enrich2$Term==Enrich1[ll,"Term"])
    outputCompare[ll,"Linked2"] <- Enrich2[ll2,"Observed genes associated with term"]
    outputCompare[ll,"Total in DB2"] <- Enrich2[ll2,"Genes in genelist that are in DB"]   
    outputCompare[ll,"Ratio"] <- (outputCompare[ll,"Linked1"]/Enrich1[ll,"Genes in genelist that are in DB"])/(outputCompare[ll,"Linked2"]/Enrich2[ll2,"Genes in genelist that are in DB"])
    outputCompare[ll,"Genes1"] <- Enrich1[ll,"Genes associated with term"]
    outputCompare[ll,"Genes2"] <- Enrich2[ll2,"Genes associated with term"]
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
  if (substr(outdir,nchar(outdir),nchar(outdir)) != "/") {outdir <- paste(outdir,"/",sep="")}
  write.table(outputCompare, 
              file=paste(outdir,db,"_TheGrandEnrich_compare_",comp_main,".txt",sep=""), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(outputCompare)
}