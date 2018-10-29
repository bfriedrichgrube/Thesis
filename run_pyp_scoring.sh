# load data from merge and scoring
cp merged_for_scoring.osw $TMPDIR
cp merged.osw $TMPDIR
mv $TMPDIR/merged.osw $TMPDIR/merged_final.osw

pyprophet score --in=$TMPDIR/merged_final.osw --apply_weights=$TMPDIR/merged_for_scoring.osw --level=ms2
pyprophet score --in=$TMPDIR/merged_final.osw --apply_weights=$TMPDIR/merged_for_scoring.osw --level=ms1
pyprophet score --in=$TMPDIR/merged_final.osw --apply_weights=$TMPDIR/merged_for_scoring.osw --level=transition

cp $TMPDIR/merged_final* .
