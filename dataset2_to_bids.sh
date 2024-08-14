for i in {001..081}; 
do 
	find /media/koba/MULTIBOOT/blindness_gradients/datasets/baltimore/share_data_101623 -type f  | grep $i | grep gz
done > nii_list.txt

for file in $(cat nii_list.txt);
do 
	sub="$(echo $file | cut  -d "/" -f 4 | cut -d "_" -f2)"

	session="$(echo $file | cut  -d "/" -f 3 | cut -d "n" -f2)"

	type="$(echo $file | cut  -d "/" -f 2 | cut -d "_" -f1)"

	#echo $sub  $session $type

	if [[ "$type" == "structure" ]]; then
		path="bids/sub-"$sub"/ses-"00$session"/anat"
		mkdir -p $path
		cp $file "$path/sub-"$sub"_ses-"00$session"_T1w.nii.gz"
		echo "$file --------------------> $path/sub-"$sub"_ses-"00$session"_T1w.nii.gz"
	elif [[ "$type" == "function" ]]; then
		path="bids/sub-"$sub"/ses-00"$session"/func"
		mkdir -p $path
		cp $file "$path/sub-"$sub"_ses-"00$session"_task-rest_bold.nii.gz"
		echo "$file --------------------> $path/sub-"$sub"_ses-"00$session"_task-rest_bold.nii.gz"
	fi
done > bids.log