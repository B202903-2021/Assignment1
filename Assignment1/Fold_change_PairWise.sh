#!/bin/bash
#Compares all samples with control and between samples (relative to control, relative to clones)
#Fold_change_PairWise.sh
count=0
for loop0 in $(cat variable_files.txt | cut -d '_' -f1 | sort | uniq);
do
count=$count
for loop in $(cat variable_files.txt | grep "$loop0");
do 
name=$(echo $loop |cut -d '_' -f1 )
name_size=${#name}
name_for_awk=${loop::-11}
#Name convention Control
name_control=${loop:$name_size} #cut control name only
# Detects files, with similar name, removes itself (loop) 
cat variable_files.txt | grep -v $loop |grep  "$name_control" > temp_name_file.txt # Removes itself then greps similar searches 
(echo -e "$loop" && cat temp_name_file.txt) > temp_output.txt && mv temp_output.txt file_list.txt


for loop2 in $(cat file_list.txt | grep -v "$loop0");
do 
count=$((count + 1))
#Name conevention samples 
name_sample=$(echo $loop2 | cut -d '_' -f1) #cuts sample names only
control_sample=$name"_"$name_sample

paste $loop $loop2 > $control_sample"_mean_combined.csv"

awk_input=$control_sample"_mean_combined.csv"
awk_output=$name_for_awk"_"$name_sample
awk -F '\t' '{ 
for (count=2; count<=NF; count++); if ($1!=0) print ($2-$1)/$1; else if ($1==0) print $2-1 ; 
}' $awk_input > $awk_output".fc.csv"
header_short=$name_for_awk"_relative_"$name_sample
echo_input=$awk_output".fc.csv"
(echo -e "$header_short" && cat $echo_input) > file1.txt && mv file1.txt $awk_output.headed.csv #adds header to each individual column 

control_sample_short=$awk_output"_fold_change.csv"
remove_centre=$(echo $name_for_awk |cut -d '_' -f2,3)
remove_complete_name=$name_sample"_"$remove_centre"_"$loop0"_fold_change.csv"

#In order to ensure duplicates will not occur, we have to check and ensure the reverse comparison is not printed
if test ! -f $remove_complete_name;
then
paste Gene_Descriptor_Field.csv $awk_output.headed.csv > $control_sample_short 
cat $control_sample_short | grep -v "Gene_ID" | sort -t$'\t' -nrk 3 > temp.csv # adds reference to individual clone comparison with WT control
sorted_fold_change=$awk_output"_sorted.csv"
(head -n1 $control_sample_short && cat temp.csv) > temp2.csv && mv temp2.csv $sorted_fold_change
fi


echo $count
done

done

return_comparison_name=$loop0
done
