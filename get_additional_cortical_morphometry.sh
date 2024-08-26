cd $SUBJECTS_DIR
sub_list=$(ls | grep -v fs)
for sub in $(echo $sub_list)
do 
    for hemi in lh rh
    do 
        echo " "
        echo "=========================== Processing $sub $hemi ==========================="
        mris_anatomical_stats -a "$SUBJECTS_DIR"/"$sub"/label/"$hemi".glasser-360_mics.annot -b "$sub" "$hemi" > "$SUBJECTS_DIR"/"$sub"/stats/"$hemi".glasser.stats
        mris_anatomical_stats -a "$SUBJECTS_DIR"/"$sub"/label/"$hemi".glasser-360_mics.annot -t "$SUBJECTS_DIR"/"$sub"/surf/"$hemi".sulc -b "$sub" "$hemi"  > "$SUBJECTS_DIR"/"$sub"/stats/"$hemi".glasser_sulc.stats
    done
done


cd $SUBJECTS_DIR
sub_list=$(ls | grep -v fs)
for sub in $(echo $sub_list)
do 
    for morph in sulc area volume 
    do 
        for hemi in lh rh
        do 
        echo "=========================================== Processing sub $sub $lh $morph"
        mri_convert $SUBJECTS_DIR/"$sub"/surf/"$hemi"."$morph" $SUBJECTS_DIR/"$sub"/surf/"$hemi"."$morph"_fsnative.mgh
        mri_surf2surf --s "$sub" --sval $SUBJECTS_DIR/"$sub"/surf/"$hemi"."$morph"_fsnative.mgh --trgsubject fsaverage5 --tval $SUBJECTS_DIR/$sub/surf/"$hemi"."$morph"_fsaverage5.mgh --hemi "$hemi"
        done
    done
done
