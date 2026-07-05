#!/bin/bash
#SBATCH --mem=4G
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --gpus 1
#SBATCH -o stdout/%x.o%j       # adjust to your own log directory (must exist before submitting)
#SBATCH -e stdout/%x.e%j
#SBATCH -p gpu

# Run MPRALegNet (https://github.com/autosome-ru/human_legnet) over one or more FASTA
# files and ensemble the per-checkpoint predictions into a single score table.
# Requires a local clone of human_legnet with a trained model, and a conda env named
# 'legnet' (see README.md).

source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate legnet

# ==== EDIT: paths & inputs ====
HUMAN_LEGNET=~/human_legnet                          # local clone of autosome-ru/human_legnet
MODEL=$HUMAN_LEGNET/models/HC/best_models            # directory of ensemble checkpoints (*.ckpt)
CONFIG=$HUMAN_LEGNET/models/HC/config.json           # model config
TMP=~/tmp/HC                                          # scratch dir for splits + intermediate predictions
FASTA=(~/Osteoarthritis/Hatzikotoulas_2025_mut.fa)   # input FASTA(s) to score
OUT_DIR=~/Osteoarthritis                             # output directory
FNAME=Hatzikotoulas_2025_mut.legnet.tsv              # output filename
BATCH_SIZE=1
CHUNK=100000
# ==============================

mkdir -p $TMP
mkdir -p $OUT_DIR

for fasta in ${FASTA[@]}
do
  SEQ=`basename -s .fa $fasta`
  OUT=`dirname $fasta`
  mkdir -p $OUT
  rm -r $TMP/${SEQ}
  mkdir -p $TMP/$SEQ

  # FASTA分割
  SPLIT_DIR=$TMP/${SEQ}/splits
  mkdir -p $SPLIT_DIR
  awk -v chunk=$CHUNK -v outprefix=$SPLIT_DIR/chunk_ '
    BEGIN { i=0; c=0; out=sprintf("%s%05d.fa", outprefix, i) }
    /^>/ {
      if (c >= chunk) {
        close(out)
        i++
        out=sprintf("%s%05d.fa", outprefix, i)
        c=0
      }
      c++
    }
    { print >> out }
  ' "$fasta"

  # 推論実行
  for ckpt in $MODEL/*.ckpt
  do
    ckpt_label=`basename -s .ckpt $ckpt`
    for chunk_fa in $SPLIT_DIR/*.fa
    do
      chunk_base=`basename -s .fa $chunk_fa`
      python $HUMAN_LEGNET/fasta_predict.py \
        --config $CONFIG \
        --model $ckpt \
        --fasta $chunk_fa \
        --out_path ${SPLIT_DIR}/${chunk_base}_${ckpt_label}.tsv \
        --device 0 \
        --batch_size $BATCH_SIZE
    done
    cat ${SPLIT_DIR}/chunk_*_${ckpt_label}.tsv|sort -k1,1 > $TMP/${SEQ}/${ckpt_label}.tsv
  done

  # 出力結合
  # 1列目（seq名）の取得
  cat $TMP/${SEQ}/${ckpt_label}.tsv | cut -f 1 > $OUT/$FNAME

  # 各モデルのスコア列を追加
  for j in $TMP/${SEQ}/best_model*.tsv
  do
    cut -f 2 $j |paste $OUT/$FNAME - > $OUT/tmp
    mv $OUT/tmp $OUT/$FNAME
  done
done
