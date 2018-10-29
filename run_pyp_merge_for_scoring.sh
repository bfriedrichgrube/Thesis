# load input files for scoring model
cp ../TCGA-*.osw $TMPDIR

pyprophet merge --out=$TMPDIR/merged_for_scoring.osw \
--subsample_ratio=0.01 $TMPDIR/*.osw

pyprophet score --in=$TMPDIR/merged_for_scoring.osw --level=ms2
pyprophet score --in=$TMPDIR/merged_for_scoring.osw --level=ms1
pyprophet score --in=$TMPDIR/merged_for_scoring.osw --level=transition

pyprophet export --format=score_plots --in=$TMPDIR/merged_for_scoring.osw

cp $TMPDIR/merged* .
