# load data
cp merged_ipf.osw $TMPDIR
mv $TMPDIR/merged_ipf.osw $TMPDIR/merged_export_matrix.osw

pyprophet export --max_rs_peakgroup_pep=1.0 --format=matrix --in=$TMPDIR/merged_export_matrix.osw

cp $TMPDIR/merged_export_matrix* .


