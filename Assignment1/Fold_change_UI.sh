#!/bin/bash
# remove 0 uninduced & searches for all uninduce into list : Search between induced and uninduced at specific times 
#Fold_change_UI.sh
count=0
for loop0 in $(cat variable_files.txt | cut -d '_' -f1 | sort | uniq); # Loops the sample type
do
for loop4 in $(cat variable_files.txt | grep $loop0); # loops the file type
do
UorI=$(cat variable_files.txt | grep $loop4 | cut -d '_' -f3 | sort | uniq) # Detects if sample is induced or uninduced

for UorIloop in $(cat variable_files.txt | grep $loop0 | grep $UorI | grep -v "0") ;
do 
for Inner_Uor_loop in $(cat variable_files.txt | grep $loop0 | grep $UorI | grep -v $UorIloop);
do
count=$((count + 1))
#Name conevention samples 

TimeC_F=$(echo $UorIloop | cut -d '_' -f2)
TimeC_T=$(echo $Inner_Uor_loop | cut -d '_' -f2)
control_sample=$TimeC_F"_"$UorI"_"$TimeC_T
paste $UorIloop $Inner_Uor_loop > $control_sample"_mean_combined.csv"

#awk input and output files 
awk_input=$control_sample"_mean_combined.csv"
awk_output=$TimeC_F"_"$UorI"_"$loop0"_"$TimeC_T

#Calculates fold change between two mean files 
awk -F '\t' '{ 
for (count=2; count<=NF; count++); if ($1!=0) print ($2-$1)/$1; else if ($1==0) print $2-1 ; 
}' $awk_input > $awk_output".fc.csv"
header_short=$name_for_awk"_relative_"$loop0
echo_input=$awk_output".fc.csv"
(echo -e "$header_short" && cat $echo_input) > file1.txt && mv file1.txt $awk_output.headed.csv #adds header to each individual column 
#Name convention for prevening duplicate
control_sample_short=$awk_output"_fold_change.csv"
remove_centre=$(echo $name_for_awk |cut -d '_' -f2,3)
remove_complete_name=$TimeC_T"_"$UorI"_"$loop0"_"$TimeC_F"_sorted_UI.csv"

#In order to ensure duplicates will not occur, we have to check and ensure the reverse comparison is not printed
if test ! -f $remove_complete_name; # If the file is not present it may create
then
#Combines fold change field to sample Id and description, ensuring correct orientation from large to small
paste Gene_Descriptor_Field.csv $awk_output.headed.csv > $control_sample_short 
cat $control_sample_short | grep -v "Gene_ID" | sort -t$'\t' -nrk 3 > temp.csv # adds reference to individual clone comparison with WT control
sorted_fold_change=$awk_output"_sorted_UI.csv"
(head -n1 $control_sample_short && cat temp.csv) > temp2.csv && mv temp2.csv $sorted_fold_change
fi


echo $count
done 
done
done
done

#There are 12 possible combinations without duplicates
