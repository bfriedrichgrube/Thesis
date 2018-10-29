
DIR=/cluster/scratch/fbetty/IPF_20180820/pyprophet
cd ${DIR}

bsub -J tric -R "rusage[mem=500000,scratch=500000]" -W 48:00 \
feature_alignment.py \
--in TCGA-*.tsv \
--out feature_alignment.csv \
--out_matrix feature_alignment_matrix.csv \
--file_format openswath \
--fdr_cutoff 0.01 \
--max_fdr_quality 0.2 \
--mst:useRTCorrection True \
--mst:Stdev_multiplier 3.0 \
--method LocalMST \
--max_rt_diff 60 \
--alignment_score 0.01 \
--frac_selected 0 \
--realign_method lowess_cython \
--disable_isotopic_grouping


