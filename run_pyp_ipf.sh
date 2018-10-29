# load data
cp merged_final.osw $TMPDIR
mv $TMPDIR/merged_final.osw $TMPDIR/merged_ipf.osw

pyprophet ipf --in=$TMPDIR/merged_ipf.osw

cp $TMPDIR/merged_ipf* .
