cd $SUBJECTS_DIR
sub_list=$(ls | grep -v fs)
for sub in $(echo $sub_list)
do 
    for hemi in lh rh
    do 
        echo " "
        echo "=========================== Processing $sub $hemi ==========================="
        mris_anatomical_stats -a "$SUBJECTS_DIR"/"$sub"/label/"$hemi".glasser-360_mics.annot -b "$sub" "$hemi" > "$SUBJECTS_DIR"/"$sub"/stats/"$hemi".glasser.stats
    done
done
