# load data
cp merged_ipf.osw $TMPDIR
mv $TMPDIR/merged_ipf.osw $TMPDIR/merged_export.osw

pyprophet export --max_rs_peakgroup_pep=1.0 --format=legacy_split --in=$TMPDIR/merged_export.osw

cp $TMPDIR/merged_export* .



