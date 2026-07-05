# This function receives a gene list and returns the enrichment and p-value of each HPO associated with the list

#genelist - a gene set of interest to analyze - can be a vector or a list 
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
#comp - a string that describes our test – can be “active vs. all” for example. 
#uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately



HPO_enrich <- function(genelist, outpath, background, minGenes, TYPorALL, HPO_dbVer, systems, bodyparts, directional, enrORdep,locORser,comp,uniq)
{
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
  if (missing(comp)) {comp <- "diff_activity"}
  if (missing(uniq)) {uniq <- T}
  #genelist <- unique(unlist(genelist))
  #background <- unique(unlist(background))
  
  #input validation
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
        } else{return(F)}
      } else{ return(F)}
    }
    valid_vec=apply(count_genelist,1,FUN=gene_validation)
    if (sum(valid_vec)<nrow(count_genelist)){
      invalid_genes=count_genelist[which(!valid_vec),'genelist']
      return(paste("Invalid input:",paste(invalid_genes,collapse = ", ")))
    }
  }
  
  # load DB  
  if (locORser=="ser"){
    load(paste("/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Annotation DBs/HPO/HPO_world/HPO_World_",TYPorALL,"_",HPO_dbVer,".RData",sep=""))
  } else{  load(paste("Gokhman lab general info/USEFUL_DATASETS/Annotation DBs/HPO/HPO_world/HPO_World_",TYPorALL,"_",HPO_dbVer,".RData",sep=""))
    }
  if (uniq){
    genelist=unique(genelist)
    background=unique(background)
  }
  
  #CONCATENATE BODY PARTS
  bodyparts_conc <- paste0(bodyparts,collapse = ", ")
  systems_conc <- paste0(systems,collapse = ", ")
  
  #REMOVE UNINFORMATIVE HPO TERMS
  hpos2remove <- c("Autosomal recessive inheritance","Autosomal dominant inheritance","X-linked dominant inheritance","X-linked recessive inheritance","Somatic mutation","X-linked inheritance","Polygenic inheritance","Multifactorial inheritance")
  idx2remove <- c()
  for (hpo in hpos2remove) {
    idx <- which(HPO_Phenos == hpo)
    if (length(idx) > 0) {idx2remove <- c(idx2remove,idx)}
  }
  HPO_World <- HPO_World[-idx2remove,]
  HPO_IDs <- HPO_IDs[-idx2remove]
  HPO_Phenos <- HPO_Phenos[-idx2remove]
  
  #LOOK FOR HPOs LINKED TO BODY PARTS
  library("readxl")
  if (HPO_dbVer == 14) {sheet <- "Pheno - build1268 - v14"
  } else if (HPO_dbVer == 13) {sheet <- "Pheno - build115 - v13"
  } else if (HPO_dbVer == 12) {sheet <- "Pheno - v12"
  } else {sheet <- "HPO"}
  #if (locORser=="ser"){
  #  HPO_ORGANizer_anno <- read_excel("/home/labs/davidgo/Collaboration/USEFUL DATASETS/Annotation DBs/HPO/specific phenotype ontology.xlsx",skip=0,sheet=sheet,col_names=T)
  #} else{  HPO_ORGANizer_anno <- read_excel("Gokhman lab general info/USEFUL DATASETS/Annotation DBs/HPO/specific phenotype ontology.xlsx",skip=0,sheet=sheet,col_names=T)
  #}
  if (locORser=="ser"){
    HPO_ORGANizer_anno <- read_excel("/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Phenotypes/Neanderthal Chimp Phenotype annotation/Reconstruction.xlsx",skip=1,sheet=sheet,col_names=T)
  } else{  HPO_ORGANizer_anno <- read_excel("Gokhman lab general info/USEFUL_DATASETS/Phenotypes/Neanderthal Chimp Phenotype annotation/Reconstruction.xlsx",skip=1,sheet=sheet,col_names=T)
  }
  HPO_ORGANizer_anno <- HPO_ORGANizer_anno[,c("HPO-Term-ID","system","organ / body part")]
  
  #FILTER TO TAKE ONLY THE USER-ENTERED BODY PART
  hpos2keep_sys <- c()
  if ("all" %in% systems) {
    hpos2keep_sys <- HPO_IDs; print("not filtering based on system")
  } else {
    print(paste0("taking only phenotypes associated with the systems: ",systems))
    for (system in systems) {
      for (row in 1:nrow(HPO_ORGANizer_anno)) {
        pos <- gregexpr(system,HPO_ORGANizer_anno[row,"system"],fixed=T)[[1]]
        if (!is.na(pos) & pos[1] != -1) {
          hpos2keep_sys <- c(hpos2keep_sys,HPO_ORGANizer_anno[row,"HPO-Term-ID"])
        }
      }
    }
    hpos2keep_sys = unlist(hpos2keep_sys)
  }
  
  hpos2keep <- c()
  if ("all" %in% bodyparts) {
    print("not filtering based on body parts")
    hpos2keep <- HPO_IDs
  }  else {
    print(paste0("taking phenotypes associated with these bodyparts: ",bodyparts_conc))
    print(paste0("taking phenotypes associated with these systems: ",systems_conc))
    for (row in 1:nrow(HPO_ORGANizer_anno)) {
      for (bodypart in bodyparts) {
        pos <- regexpr(bodypart, HPO_ORGANizer_anno[row,"organ / body part"],fixed=T)[[1]]
        if (!is.na(pos) & pos[1] != -1) {
          hpos2keep <- c(hpos2keep,HPO_ORGANizer_anno[row,"HPO-Term-ID"])
        }
      }
    }
    hpos2keep=unlist(hpos2keep)
  }
  
  hpo2keep_directional <- HPO_IDs
  if (directional) {
    load(paste("DerivedTraits_",TYPorALL,"_",HPO_dbVer,".RData",sep=""))
    hpo2keep_directional <- as.character(unique(gene2hpo2trait[,"HPO ID"]))
  }
  
  hpos2keep <- intersect(intersect(hpos2keep,hpos2keep_sys),hpo2keep_directional)
  bodypartHPOs <- intersect(HPO_IDs,hpos2keep)
  idxhpos2keep <- match(bodypartHPOs,HPO_IDs)
  HPO_IDs <- HPO_IDs[idxhpos2keep]
  HPO_Phenos <- HPO_Phenos[idxhpos2keep]
  HPO_World <- HPO_World[idxhpos2keep,]
  hposPerGene <- colSums(HPO_World)
  idx <- which(colSums(HPO_World) == 0)
  if (length(idx) > 0) {
    HPO_World <- HPO_World[,-idx]
    HPO_Genes <- HPO_Genes[-idx]
    HPO_ENTREZ_IDs <- HPO_ENTREZ_IDs[-idx]
  }
  
  # CREATE BACKGROUND
  if (length(background) > 0) {
    idx <- c()
    for (gg in background) {
      gene <- toupper(gg)
      ii <- which(HPO_Genes == gene)
      if (length(ii) > 0) {
        #idx <- rbind(idx,ii)
        idx <- c(idx,ii)
      }
    }
    idx <- unlist(idx)
    HPO_World <- HPO_World[,idx]
    HPO_Genes <- HPO_Genes[idx]
    HPO_ENTREZ_IDs <- HPO_ENTREZ_IDs[idx]
  }
  
  # CREATE USER-ENTERED SUB-TABLE OF HPO_WORLD WITH GENE-HPO ASSOCIATIONS BASED ON ENTERED GENE LIST
  idx <- c()
  for (gg in genelist) {
    gene <- toupper(gg)
    ii <- which(HPO_Genes == gene)
    if (length(ii) > 0) {
      ind=ii[1] #in case of redundancy - add only one row each time
      idx <- c(idx,ind)
    }
  }
  idx <- unlist(idx)
  HPO_Entered <- HPO_World[,idx]
  HPO_Entered_Genes <- HPO_Genes[idx]
  HPO_Entered_ENTREZ_IDs <- HPO_ENTREZ_IDs[idx]
  
  #GO OVER EACH HPO AND COMPUTE P-VALS
  # N_world - number of genes found both within background gene list (active genes) and HPO's gene list
  # N_entered - number of genes found both within test gene list (differentially active) and HPO's gene list
  # n_world - number of genes from background gene list found to be associated with a certain phenotype
  # Observed (n_entered) - number of genes from test gene list found to be associated with a certain phenotype
  
  cols <- c("Pheno","HPO_ID","Enrichment","N_world","n_world","N_entered","Observed (n_entered)","Expected","Genes","P-value","FDR")
  StatMat <- matrix(NA,nrow=nrow(HPO_World),ncol=length(cols))
  colnames(StatMat) <- cols
  StatMat <- as.data.frame(StatMat)
  
  N_world <- ncol(HPO_World)
  N_entered <- ncol(HPO_Entered)
  StatMat$N_world=N_world
  StatMat$N_entered<- N_entered
  if (length(N_entered) > 0) {
    for (hpo in 1:nrow(HPO_World)) {
      #StatMat$N_world[hpo] <- N_world      
      #StatMat$N_entered[hpo] <- N_entered
      StatMat$Pheno[hpo] <- HPO_Phenos[hpo]
      StatMat$HPO_ID[hpo] <- HPO_IDs[hpo]
      n_world <- sum(HPO_World[hpo,])
      StatMat$n_world[hpo] <- n_world
      n_entered <- sum(HPO_Entered[hpo,])
      StatMat$Enrichment[hpo] <- (n_entered/N_entered)/(n_world/N_world)
      StatMat[hpo,"Observed (n_entered)"] <- n_entered
      idx <- which(HPO_Entered[hpo,] == 1)
      StatMat[hpo,"Genes"] <- paste(HPO_Entered_Genes[idx],collapse=", ") # ASK SIMON: unique or not?
      StatMat$Expected[hpo] <- (n_world/N_world)*N_entered  #number of active genes associated to phenotype times all diff active genes over all active genes
      
      mat_hpo <- matrix(c(n_entered, n_world - n_entered, #number of genes associated to the phenotype that are active but not diff active
                          N_entered - n_entered,          #number of diff active genes that aren't associated to the phenotype
                          (N_world-n_world) - (N_entered - n_entered)),  #number of active, non-diff active genes that aren't associated to the phenotype
                        nrow=2, byrow=T)
      row_sums=rowSums(mat_hpo)
      col_sums=colSums(mat_hpo)
      tot_sum=sum(mat_hpo)
      exp_mat=matrix(c(row_sums[1]*col_sums[1]/tot_sum,row_sums[2]*col_sums[1]/tot_sum,
                       row_sums[1]*col_sums[2]/tot_sum,row_sums[2]*col_sums[2]/tot_sum),nrow=2)
      bool_mat=exp_mat<mat_hpo
      bool_mat=bool_mat*1
      mat_hpo_corrected=mat_hpo+bool_mat
      if (enrORdep == "enr") {
        StatMat$'P-value'[hpo] <- stats::fisher.test(mat_hpo_corrected, alternative = "g")$p.value
      } else if (enrORdep == "dep") {
        StatMat$'P-value'[hpo] <- stats::fisher.test(mat_hpo_corrected, alternative = "l")$p.value
      } else if (enrORdep == "both") {
        StatMat$'P-value'[hpo] <- stats::fisher.test(mat_hpo_corrected, alternative = "t")$p.value
      }
      
    }
  }
  #FDR
  idx <- which(StatMat[,"Observed (n_entered)"] >= minGenes)
  if (length(idx) > 0) {StatMat$FDR[idx] <- p.adjust(StatMat$'P-value'[idx], method = "BH")}
  
  idx <- which(StatMat$FDR < 0.05)
  print(paste("Found ",length(idx)," significant HPO terms",sep=""))
  StatMat=StatMat[order(StatMat$FDR),]
  
  #SAVE OUTPUT
  if (substr(outpath,nchar(outpath)-1,nchar(outpath)) != "/") {outpath <- paste(outpath,"/",sep="")}
  write.table(StatMat, 
              file=paste(outpath,"HPO_enrich_","v",HPO_dbVer,"_","minGenes_",
                         minGenes,"_",systems,"_","bp_",bodyparts,"_","freq_",
                         TYPorALL,"_",enrORdep,"_",comp,".txt",sep=""), 
              sep="\t", quote = F,na = "", row.names = F, col.names = T)
  
  return(StatMat)
}
