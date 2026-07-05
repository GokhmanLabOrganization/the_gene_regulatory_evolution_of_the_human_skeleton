# This is an example of how to run the TheGrandEnrich_enrich.R script but the flow is similar in all other scripts

# make sure you have the readxl library
library("readxl")

# load script
grandEnrichDir <- "/home/labs/davidgo/Collaboration/backup/Lab_tools/Gokhman_Enrich/final_scripts/TheGrandEnrich_enrich.R"
source(grandEnrichDir)

# set arguments - notice that the needed arguments vary between scripts
# don't forget to change the lists and directory
db="GO_BP"
locORser = "ser"
genelist = c("BRCA1","BRCA2","PALB2")
background = c("BRCA2","PALB2","CDH1","BRCA1","EIF4A3","BRCA1","PALB2")
enrORdep="both"
minGenes = 1
outdir = "/home/labs/davidgo/Collaboration/Lab_Tools/Enrichment"
comp="TestRun"
uniq=T

# run script
result=TheGrandEnrich_enrich(db=db,genelist=genelist,background=background,enrORdep=enrORdep,minGenes=minGenes,outdir=outdir,comp=comp,uniq=uniq)
head(result)
