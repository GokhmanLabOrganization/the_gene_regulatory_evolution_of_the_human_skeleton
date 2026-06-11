# run R script aseByPosition for all given samples
# input for R is: ase_path, outpath and window size

# combos=(
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC1-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC1-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC2-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC2-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC3-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC3-3_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/ASE/human/hg1-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/ASE/human/hg1-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/ASE/human/hg2-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/ASE/human/hg2-3_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL28-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL28-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL31-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|20"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL31-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|20"
# )

# combos=(
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC1-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC2-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/ASE/human/HC3-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/ASE/human/hg1-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/ASE/human/hg2-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL28-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL31-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|50"
# "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/Cartilage_Organoids/backup_counts_and_deseq2/S2-HC-HL1-30_ase_by_reads_merged_humanref.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/Cartilage_Organoids/aneuploidy_plots/|50"
# )

# for combo in "${combos[@]}"; do
#   IFS="|" read -r aseFile outpath window <<< "$combo"

#   echo "Running:"
#   echo "  aseFile = $aseFile"
#   echo "  outpath = $outpath"
#   echo "  window  = $window" 
#   echo

#   Rscript ./aseByPosition.R "$aseFile" "$outpath" "$window"
# done


GTF="/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/Cartilage_Organoids/raw_data/Cartilage_Organoids_Hybrids_processed/references/human/human.gtf"
PYTHON="python"
ASE_SCRIPT="/home/labs/davidgo/galbo/backup/hybrids-data-scripts/counts/aneuploidy_analysis.py"

# List of jobs: each entry is "counts|output|window"
JOBS=(
#"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/ASE/human/HC1-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/aneuploidy_plots/|50"
#"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/ASE/human/HC2-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/aneuploidy_plots/|50"
#"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/ASE/human/HC3-2_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/aneuploidy_plots/|50"
"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/ASE/human/hg1-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/aneuploidy_plots/|50"
"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/ASE/human/hg2-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/ExpLBM/outputs_05Jan2026_humanMPRA_draft1/aneuploidy_plots/|50"
#"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL28-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|50"
#"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/ASE/human/HGL31-1_ase_by_reads_merged.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_gorilla/iPSCs/outputs/aneuploidy_plots/|50"
#"/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/Cartilage_Organoids/backup_counts_and_deseq2/S2-HC-HL1-30_ase_by_reads_merged_humanref.txt|/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids/human_chimp/Cartilage_Organoids/aneuploidy_plots/|50"
)

# ---- RUN ----

for job in "${JOBS[@]}"; do
  IFS="|" read -r COUNTS OUTFILE WINDOW <<< "$job"

  echo "Running aneuploidy analysis:"
  echo "  counts : $COUNTS"
  echo "  out    : $OUTFILE"
  echo "  window : $WINDOW"
  echo

  "$PYTHON" "$ASE_SCRIPT" \
    --gff "$GTF" \
    --counts "$COUNTS" \
    --window "$WINDOW" \
    --out "$OUTFILE"
done

echo "All jobs finished."

