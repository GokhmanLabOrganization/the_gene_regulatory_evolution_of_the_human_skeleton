library(readxl)
library(MPRAnalyze)
library(glue)

args = commandArgs(trailingOnly=TRUE)

lib = args[4]
adaptor = args[5]
cells = args[6]
input = args[7]

# read in full data
mpra_colannot <- read.table(
	glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_comparative/{cells}_{lib}{adaptor}_col_annot_comparative_improved.txt"),
	header=T, sep='\t', stringsAsFactors=F, quote="", row.names=1)
mpra_data_DNA <- read.table(
	glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_comparative/{cells}_{lib}{adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt"),
	header=T, sep='\t', stringsAsFactors=F, quote="", row.names=1)
mpra_data_RNA <- read.table(
	glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_comparative/{cells}_{lib}{adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt"),
	header=T, sep='\t', stringsAsFactors=F, quote="", row.names=1)

print('done reading in data for full model')

# convert to matrix from data.frame
ce.dnaCounts <- as.matrix(as.data.frame(mpra_data_DNA))
ce.rnaCounts <- as.matrix(as.data.frame(mpra_data_RNA))

# assign column annotations, can stay as data.frame
ce.colAnnot <- mpra_colannot


### WITHOUT CONTROLS ###
# initialize MpraObject without controls
obj <- MpraObject(dnaCounts = ce.dnaCounts, rnaCounts = ce.rnaCounts, 
                  colAnnot = ce.colAnnot)
                  
print('made big MPRAobject')
# get depth factors
obj <- estimateDepthFactors(obj, lib.factor = "replicate")

d_depth <- obj@dnaDepth
r_depth <- obj@rnaDepth
print('calculated depth factors')

# read the smaller dataframe
dnaname = args[1]
rnaname = args[2]
name = args[3]

mpra_data_RNA_small <- read.table(glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_comparative/rna_chunks_{input}/{rnaname}"), header=F, sep='\t', stringsAsFactors=F, quote="", row.names=1)
mpra_data_DNA_small <- read.table(glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_comparative/dna_chunks_{input}/{dnaname}"), header=F, sep='\t', stringsAsFactors=F, quote="", row.names=1)
print('done reading in small dataframes')

d = as.matrix(as.data.frame(mpra_data_DNA_small))
r = as.matrix(as.data.frame(mpra_data_RNA_small))

obj2 <- MpraObject(dnaCounts = d, rnaCounts = r, colAnnot = ce.colAnnot)
print('made small MPRAobject')
obj2 <- setDepthFactors(obj2, dnaDepth = d_depth, rnaDepth = r_depth)
# obj2 <- analyzeComparative(obj2, 
                          # dnaDesign = ~ barcode_allele + replicate, 
                          # rnaDesign = ~ allele + replicate,
                          # correctControls = FALSE,
                          # reducedDesign = ~ replicate)
# 14.09.23 - changed the design so that replicates are not modelled in the RNA design
obj2 <- analyzeComparative(obj2, 
                          dnaDesign = ~ barcode_allele + replicate, 
                          rnaDesign = ~ allele,
                          correctControls = FALSE,
                          reducedDesign = ~1)
print('analyze comparative done')

# test for differential activity
res <- testLrt(obj2)

# save results
write.table(res, file=glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_comparative/res_{input}/{name}.txt"), 
            sep="\t", quote = F,na = "", row.names = T, col.names = T)
