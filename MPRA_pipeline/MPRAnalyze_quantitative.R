# install.packages('devtools')

# Load packages

library("devtools")

# install the correct version of MPRAnalyze
devtools::install_github("YosefLab/MPRAnalyze", force=TRUE)

library("MPRAnalyze")
library("glue")

sessionInfo()
packageVersion("MPRAnalyze") #should be commented

args = commandArgs(trailingOnly=TRUE)
print(args)
lib <- args[1]
adaptor <- args[2]
cells <- args[3]



# Read in data
 
annot = read.csv(
	glue("./{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_col_annot_quantitative.txt"), 
	sep = '\t', 
	row.names = 1
	)

annot$replicate<-as.factor(annot$replicate)
annot$barcode<-as.factor(annot$barcode)

dcounts_table = read.table(
	glue("./{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_DNA_MPRAnalyze_quantitative.txt"), 
	header=T, 
	sep='\t', 
	stringsAsFactors=F, 
	quote="", 
	row.names=1
	)
rcounts_table = read.table(
	glue("./{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_RNA_MPRAnalyze_quantitative.txt"), 
	header=T, 
	sep='\t', 
	stringsAsFactors=F, 
	quote="", 
	row.names=1
	)

# Coerce to matrices for later

dcounts = as.matrix(as.data.frame(dcounts_table))
rcounts = as.matrix(as.data.frame(rcounts_table))
neg_ctrls = grepl("scrambled", row.names(dcounts))   # could be done on the dcounts or the rcounts, they're the same
non_neg_ctrls = !grepl("scrambled", row.names(dcounts))  # all the rows that aren't negative controls

# Estimate and extract depth factors on entire dataset

obj <- MpraObject(dnaCounts = dcounts, rnaCounts = rcounts, colAnnot = annot, controls=neg_ctrls)
obj <- estimateDepthFactors(obj, lib.factor = "replicate", which.lib="both")
obj <- analyzeQuantification(
	obj,
	rnaDesign = ~1,
	dnaDesign = ~ replicate
	)
alpha_res <- getAlpha(obj)
res <- testEmpirical(obj, useControls=TRUE)

# Add alphas to res
identical(rownames(alpha_res), rownames(res)) # should be TRUE
res$alpha <- cbind(alpha_res$alpha)

# save(list = c("res", "alpha_res","obj"), file = "data") #KL 30.1.23 - not sure why we this is done. not really needed

# added by SF: this line was previously run before the mpranalyze model (on dcounts), before running mpranalyze, which ommits several enahncers. Now moved here and run on res
non_neg_ctrls = !grepl("scrambled", row.names(res))  # all the rows that aren't negative controls

# Fix FDRs to correct values excluding the negative controls from that calculation (we don't need to include them in the tests)
# res$fdr.mad[non_neg_ctrls] <- p.adjust(res$pval.mad[non_neg_ctrls], method='BH')
# KL 30.01.23: at the moment: include all oligos in fdr correction, later might remove scrambled and negative controls
res$fdr.mad <- p.adjust(res$pval.mad, method='BH')

# Write results to out

write.table(
			alpha_res, 
			      file=glue("./{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_alpha_full_quantitative_nobc_annot_factor_new_design.txt"), 
			      sep="\t", 
            quote = F, 
            na = "", 
            row.names = T, 
            col.names = T
            )

write.table(
			res, 
			      file=glue("./{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt"),
            sep="\t", 
            quote = F, 
            na = "", 
            row.names = T, 
            col.names = T
            )

# # divide up results into positive, negative controls, and ancestral/derived oligos (easier to input into VARs_anno this way)
# pos_ctrls = grepl("ctrl", row.names(res))
# neg_ctrls = grepl("scrambled", row.names(res))
# ancestral_oligos = grepl("Ancestral", row.names(res))
# derived_oligos = grepl("Derived", row.names(res))

# pos = res[pos_ctrls, ]
# neg = res[neg_ctrls, ]
# ancestral = res[ancestral_oligos, ]
# derived = res[derived_oligos, ]

# # write them all to tables
# write.table(
  # pos, 
  # file=glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_results_full_quantitative_nobc_fdr_pos_only.txt"), 
  # sep="\t", 
  # quote = F, 
  # na = "", 
  # row.names = T, 
  # col.names = T
# )

# write.table(
  # neg, 
  # file=glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_results_full_quantitative_nobc_fdr_neg_only.txt"), 
  # sep="\t", 
  # quote = F, 
  # na = "", 
  # row.names = T, 
  # col.names = T
# )

# write.table(
  # ancestral, 
  # file=glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_results_full_quantitative_nobc_fdr_ancestral.txt"), 
  # sep="\t", 
  # quote = F, 
  # na = "", 
  # row.names = T, 
  # col.names = T
# )

# write.table(
  # derived, 
  # file=glue("/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{lib}{adaptor}/output/mpranalyze_quantitative/{cells}_{lib}{adaptor}_results_full_quantitative_nobc_fdr_derived.txt"), 
  # sep="\t", 
  # quote = F, 
  # na = "", 
  # row.names = T, 
  # col.names = T
# )
