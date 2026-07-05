# MPRALegNet prediction

MPRALegNet is a deep-learning model that predicts regulatory (MPRA) activity from
sequence. In the paper it was trained on the chondrocyte MPRA data and used to (i)
assess predicted-vs-measured activity (Extended Data, quality control) and (ii)
perform *in silico* saturation mutagenesis on the *EVC* differentially active cCRE
(Methods, *"MPRALegNet training and in silico saturation mutagenesis"*; Fig. 5c).

## What's here

`legnet_prediction.sh` — a SLURM driver that runs MPRALegNet inference over one or
more FASTA files and ensembles the per-checkpoint predictions into a single score
table. For each input FASTA it:
1. splits the FASTA into chunks of `CHUNK` sequences,
2. runs `human_legnet/fasta_predict.py` for every model checkpoint (`$MODEL/*.ckpt`)
   on every chunk,
3. concatenates and merges the results so the output has one column of sequence
   names followed by one score column per checkpoint.

Edit the config block at the top of the script (marked `==== EDIT ====`) to set the
model, input FASTA(s), and output paths.

## Dependencies (external)

This is a **driver only**. The model and inference code come from the MPRALegNet
repository, which is **not bundled here**:

* **[autosome-ru/human_legnet](https://github.com/autosome-ru/human_legnet)** — clone
  locally and point `HUMAN_LEGNET` at it. It provides `fasta_predict.py` and the
  trained model (`models/HC/best_models/*.ckpt`, `models/HC/config.json`). The paper
  trained MPRALegNet on the chondrocyte MPRA "with minor modifications" to this repo.
* A conda environment named `legnet` with the model's dependencies (PyTorch + GPU).
* A GPU (the script requests `--gpus 1` on the `gpu` partition).

> Note: model **training** and the *EVC* **in silico saturation mutagenesis** (Fig. 5c)
> were run within the human_legnet environment and are not included in this repository.

## Output

A single TSV (`$OUT_DIR/$FNAME`): column 1 = sequence name, subsequent columns =
predicted activity per model checkpoint (average across columns for the ensemble
prediction).
