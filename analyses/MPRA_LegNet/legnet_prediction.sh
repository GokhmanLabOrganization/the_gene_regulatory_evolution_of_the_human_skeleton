#!/bin/bash
#SBATCH --mem=4G
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --gpus 1
#SBATCH -o /home/zicong_zhang/stdout/%x.o%j
#SBATCH -e /home/zicong_zhang/stdout/%x.e%j
#SBATCH -p gpu

source ~/miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate legnet

MODEL=~/human_legnet/models/HC/best_models
CONFIG=~/human_legnet/models/HC/config.json
TMP=~/tmp/HC
FASTA=(~/Osteoarthritis/Hatzikotoulas_2025_mut.fa)
OUT_DIR=~/Osteoarthritis
FNAME=Hatzikotoulas_2025_mut.legnet.tsv
BATCH_SIZE=1
CHUNK=100000

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
      python ~/human_legnet/fasta_predict.py \
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
