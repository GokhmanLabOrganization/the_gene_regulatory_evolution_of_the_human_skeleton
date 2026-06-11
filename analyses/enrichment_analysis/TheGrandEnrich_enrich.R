#THE FUNCTION RECEIVES A GENE LIST AND COMPARES ITS ENRICHMENT/DEPLETION COMPARED TO THE BACKGROUND (GENOMIC OR USER-SUPPLIED BACKGROUND LIST)

#db: GO_BP, GO_MF, GO_CC, KEGG, OMIM, REACTOME, GAD_Disease
#genelist - a gene set of interest to analyze - can be a vector or a list 

#OPTIONAL INPUTS:
#background - a larger set of genes, must include all genes in genelist or empty
#enrORdep - the direction of the test, can be enr for enrichment (higher probability to see term in genelist), dep (lower probability) or both (a non directional test) 
#minGenes - minimum genes found per term. Terms with fewer genes will not be FDR-corrected. If the user chose
#outdir - a path for saving the results file
#locORser - running locally ("loc") or on the server ("ser") - affects on DBdir. If you run locally - set working directory to: /Users/<username>/Dropbox (Weizmann Institute)/
#comp - a string that describes our test – can be “active vs. all” for example. 
#uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately



TheGrandEnrich_enrich <- function(db,genelist,background,minGenes,outdir,enrORdep,locORser,comp,uniq)
{
  if (missing(db)) {db <- "GO_BP"}
  if (missing(background)) {background <- c()} 
  if (missing(minGenes)) {minGenes <- 0}
  if (missing(outdir)) {outdir <- "/Users/dgokhman/HUJI drive/R/R outputs/"}
  if (missing(enrORdep)) {enrORdep <- "both"}
  if (missing(locORser)) {locORser <- "ser"}
  if (missing(comp)) {comp <- "diff_activity"}
  if (missing(uniq)) {uniq <- T}
  
  #load DB
  geneCol <- 1
  catCol <- 2
  if (locORser=="ser"){
    DBdir <- '/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Annotation DBs/'
    
  } else {
    DBdir <- 'Gokhman lab general info/USEFUL_DATASETS/Annotation DBs/'
  }
  if (db == "GO_BP") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt",sep="")
  } else if (db == "GO_MF") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GOTERM_MF_FAT.txt",sep="")
  } else if (db == "GO_CC") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GOTERM_CC_FAT.txt",sep="")
  } else if (db == "REACTOME-DAVID") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2REACTOME_PATHWAY.txt",sep="")
  } else if (db == "GAD_Disease") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GAD_DISEASE.txt",sep="")
  } else if (db == "GAD_Disease_class") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2GAD_DISEASE_CLASS.txt",sep="")
  } else if (db == "OMIM") {
    input_path_gene2category <- paste(DBdir,"DAVIDKnowledgebase/DAVIDKnowledgebase_current/OFFICIAL_GENE_SYMBOL2OMIM_DISEASE.txt",sep="")
  } else if (db == "KEGG") {
    input_path_gene2category <- paste(DBdir,"KEGG/current/KEGG_pathway_NEW.txt",sep="")
    geneCol <- 2
    catCol <- 5
  }
  if (db %in% c("GO_BP","GO_MF","GO_CC")) {
    gene2cat_World <- read.table(input_path_gene2category, header=F, sep='\t', stringsAsFactors=F, quote="")
  } else {
    gene2cat_World <- read.table(input_path_gene2category, header=T, sep='\t', stringsAsFactors=F, quote="")
  }
  gene2cat_World_unique=gene2cat_World[match(unique(gene2cat_World[,geneCol]), gene2cat_World[,geneCol]),]
  
  # Input validation - checks that each gene that appears in genelist appears in background
  
  if (class(background)=="list"){
    background=unlist(background)
  }
  if (class(genelist)=="list"){
    genelist=unlist(genelist)
  }
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
        } else{
          return(F)}
      } else{ return(F)}
    }
    valid_vec=apply(count_genelist,1,FUN=gene_validation)
    if (sum(valid_vec)<nrow(count_genelist)){
      invalid_genes=count_genelist[which(!valid_vec),'genelist']
      return(paste("Invalid input:",paste(invalid_genes,collapse = ", ")))
    }
  }
  
  
  #intersect DB with background entered
  if (length(background) > 0) {
    idxBack <- which(toupper(gene2cat_World[,geneCol]) %in% toupper(unlist(background)))
    countsTable_world=table(background)
    gene2cat_World <- gene2cat_World[idxBack,]
    gene2cat_World$counts=unname(sapply(gene2cat_World[,geneCol],FUN=function(x){(countsTable_world[x])})) #count the number of occurences of each gene in the list
    gene2cat_World_unique=gene2cat_World[match(unique(gene2cat_World[,geneCol]), gene2cat_World[,geneCol]),]
  }
  
  #intersect genelist with DB
  uniqueGenes <- unique(unlist(genelist))
  idxGenes <- which(unlist(gene2cat_World[,geneCol]) %in% uniqueGenes)
  countsTable_inter=table(genelist) 
  gene2cat_intersectedWithGenelist <- gene2cat_World[idxGenes,]
  gene2cat_intersectedWithGenelist$counts=unname(sapply(gene2cat_intersectedWithGenelist[,geneCol],FUN=function(x){(countsTable_inter[x])})) #count the number of occurences of each gene in the list
  intersectedGenes <- intersect(unlist(gene2cat_World[,geneCol]),unlist(genelist))
  idxGenes <- match(intersectedGenes,gene2cat_intersectedWithGenelist[,geneCol])
  gene2cat_intersectedWithGenelist_unique <- gene2cat_intersectedWithGenelist[idxGenes,]
  
  if(length(idxGenes)==0){
    return("No genes from genelist are in DB")
  }

  colnames <- c("DB","Term","Genes in genelist","Genes in background",
                "Genes in genelist that are in DB",
                "Genes in background that are in DB",
                "Observed genes associated with term",
                "Observed genes associated with term, background","Expected genes associated with term",
                "FoldChange","Genes associated with term","P-value","FDR")
  uniqueTerms <- unique(unlist(gene2cat_World[,catCol]))
  uniqueGenes_inDB <- unique(unlist(gene2cat_World[,geneCol]))
  numredundantGenes_inDB <- sum(gene2cat_World[match(unique(gene2cat_World[,geneCol]),gene2cat_World[,geneCol]),'counts']) #sum the number of occurences for all genes in background db intersect
  outTable <- as.data.frame(matrix(NA,nrow=length(uniqueTerms),ncol=length(colnames)))
  colnames(outTable) <- colnames
  
  #go over terms - run Fisher's exact on each term
  # columns with identical values
  outTable["DB"] <- db
  if(uniq){
    outTable["Genes in genelist"]<-length(uniqueGenes)
    outTable["Genes in background"]<-length(unique(background))
    outTable["Genes in genelist that are in DB"] <- nrow(gene2cat_intersectedWithGenelist_unique)
    outTable["Genes in background that are in DB"] <-nrow(gene2cat_World_unique)
  }
  else{
    outTable["Genes in genelist"]<-length(genelist)
    outTable["Genes in background"]<-length(background)
    outTable["Genes in genelist that are in DB"] <- sum(gene2cat_intersectedWithGenelist_unique$counts)
    outTable["Genes in background that are in DB"] <- sum(gene2cat_World_unique$counts)
  }
  ListRatio <- outTable[1,"Genes in genelist that are in DB"]/outTable[1,"Genes in background that are in DB"][1]
  for (ll in 1:length(uniqueTerms)) {
    term <- uniqueTerms[ll]
    outTable[ll,"Term"] <- term
    idxW <- which(gene2cat_World[,catCol] == term) # num of unique genes linked to term in the DB
    numidxW_redundant=sum(gene2cat_World[gene2cat_World[,catCol] == term,]$counts) # num of redundant genes linked to term in the DB
    idxO <- which(gene2cat_intersectedWithGenelist[,catCol] == term) # num of genes linked to term in the user-entered list
    numidxO_redundant=sum(gene2cat_intersectedWithGenelist[gene2cat_intersectedWithGenelist[,catCol] == term,]$counts)
    outTable[ll,"Genes associated with term"] <- paste(gene2cat_intersectedWithGenelist[idxO,geneCol],collapse=",") # genes linked to term in the user-entered list
    numGenes4term_DB <- length(idxW)
    numGenes4term_DB_redundant <- numidxW_redundant
    if (length(idxW) > 0) {
        #outTable[ll,"Enrichment"] <- outTable[ll,"Observed"]/outTable[ll,"Expected"]
        #if (!is.na(outTable[ll,"Enrichment"])) {
          if (uniq){
          outTable[ll,"Observed genes associated with term"] <- length(idxO)
          outTable[ll,"Observed genes associated with term, background"] <- numGenes4term_DB
          mat_grand <- matrix(c(outTable[ll,"Observed genes associated with term"], outTable[ll,"Observed genes associated with term, background"] - outTable[ll,"Observed genes associated with term"],
                                outTable[ll,"Genes in genelist that are in DB"] - outTable[ll,"Observed genes associated with term"],
                                (length(uniqueGenes_inDB)-
                                   outTable[ll,"Observed genes associated with term, background"]) -
                                  (outTable[ll,"Genes in genelist that are in DB"] - outTable[ll,"Observed genes associated with term"])),
                              nrow=2, byrow=T)
      }else{
        outTable[ll,"Observed genes associated with term"] <- numidxO_redundant
        outTable[ll,"Observed genes associated with term, background"] <- numGenes4term_DB_redundant 
          mat_grand <- matrix(c(outTable[ll,"Observed genes associated with term"],
                                outTable[ll,"Observed genes associated with term, background"] - outTable[ll,"Observed genes associated with term"],
                                outTable[ll,"Genes in genelist that are in DB"] - outTable[ll,"Observed genes associated with term"],
                                (numredundantGenes_inDB-outTable[ll,"Observed genes associated with term, background"]) - (outTable[ll,"Genes in genelist that are in DB"] - outTable[ll,"Observed genes associated with term"])),
                              nrow=2, byrow=T)
      }
        outTable[ll,"Expected genes associated with term"] <- ListRatio * outTable[ll,"Observed genes associated with term, background"]
        outTable[ll,"FoldChange"] <- outTable[ll,"Observed genes associated with term"]/outTable[ll,"Expected genes associated with term"]
        #print(outTable[ll,])
        #print("matgrand")
        #print(mat_grand)
        row_sums=rowSums(mat_grand)
        col_sums=colSums(mat_grand)
        tot_sum=sum(mat_grand)
        exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                         row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
        bool_mat=exp_mat<mat_grand
        bool_mat=bool_mat*1
        mat_grand_corrected=mat_grand+bool_mat
        #print(mat_grand_corrected)
        if (enrORdep == "enr") {
          outTable[ll,"P-value"] <- stats::fisher.test(mat_grand_corrected, alternative = "g")$p.value
        } else if (enrORdep == "dep") {
          outTable[ll,"P-value"] <- stats::fisher.test(mat_grand_corrected, alternative = "l")$p.value
        } else if (enrORdep == "both") {
          outTable[ll,"P-value"] <- stats::fisher.test(mat_grand_corrected, alternative = "t")$p.value
        }
        
        }
    #}
      
  }
  
  #fdr
  
  pvals_temp <- outTable[,"P-value"]
  if (enrORdep == "enr") {
    idxminGenes <- which(outTable[,"Observed genes associated with term"] <= minGenes)
    } 
  else if (enrORdep == "dep") {
    idxminGenes <- which(outTable[,"Expected genes associated with term"] <= minGenes)
    } 
  else if (enrORdep == "both"){
    idxminGenes <- which((outTable$FoldChange > 1 & outTable["Observed genes associated with term"] <= minGenes) | (outTable$FoldChange < 1 & outTable$Expected <= minGenes))
    }
  pvals_temp[idxminGenes] <- NA
  outTable[,"FDR"] <- p.adjust(pvals_temp,method="BH")
  outTable=outTable[order(outTable$FDR),]
  
  #WRITE FILE
  
  if (substr(outdir,nchar(outdir)-1,nchar(outdir)) != "/") {outdir <- paste(outdir,"/",sep="")}
  write.table(outTable, file=paste(outdir,db,"_","minGenes_",minGenes,"_",enrORdep,"_",comp,
                                   "_uniq_",uniq,"_TheGrandEnrich_enrich.txt"
                                   ,sep=""), sep="\t", quote = F,na = "",
              row.names = F, col.names = T)
  
  return(outTable)
}

