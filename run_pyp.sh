# Run Pyprophet completely

OUT_DIR=/cluster/scratch/fbetty/IPF_20180820/pyprophet
cd ${OUT_DIR}

# merge and score input files
bsub -J pypMerge_for_scoring -R "rusage[mem=300000,scratch=300000]" -W 120:00 ./run_pyp_merge_for_scoring.sh

# merge all files
bsub -J pypMerge -w "done(pypMerge_for_scoring)" -R "rusage[mem=300000,scratch=300000]" -W 48:00 ./run_pyp_merge_fast.sh

# run scoring on full merged file
bsub -J scoring -w "done(pypMerge)" -R "rusage[mem=500000,scratch=500000]" -W 120:00 ./run_pyp_scoring.sh
bsub -J ipf -w "done(scoring)" -R "rusage[mem=500000,scratch=500000]" -W 120:00 ./run_pyp_ipf.sh

# export report
bsub -J export -w "done(ipf)" -R "rusage[mem=500000,scratch=500000]" -W 24:00 ./run_pyp_export.sh
bsub -J export_matrix -w "done(export)" -R "rusage[mem=500000,scratch=500000]" -W 24:00 ./run_pyp_export_matrix.sh
