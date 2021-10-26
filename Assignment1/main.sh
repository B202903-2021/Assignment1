#!bin/bash
########Stage 1: Quality check to assess the quality of the raw sequence data #########
fastqc -o . /localdisk/data/BPSM/AY21/fastq/*.fq.gz

#-------------------------------------------Stage1

#Checking which sequences pass which control tests
for fastqc_loop in $(ls -1 *fastqc.zip)
do #Initial for zip loop list
unzip $fastqc_loop # unzip the file in order to acess the output files
#Create parameters to create files of sequences that failed specific quality controls
Edited_Name=${fastqc_loop::-4}
QC_name_input=$Edited_Name"/summary.txt" #removes the .zip on the end to enable the acess of summary file within the unzipped file 
# We want to loop the quality control to check if the 1st column is a pass or not
count=$0
qcount=$0 # used to count what QC test we are on
for quality_control_loop in $(cat $QC_name_input | cut -f1) # loops through each quality control for that specific fastqc file
#-f1 was used as it gives the correct number of QC stages, -f2 instead results in counting each individual word. leads to Wrong number of loops.
# also can act as Pass or Fail variable
do #QC do
qcount=$((qcount + 1)) # allows use to name files based on the stage of QC

QC_control=$(cat $QC_name_input | cut  -f2) # This sets the variable to the QC value (eg. per base sequence quality)
# it automatically cuts by tab delimited, thus not specification for what divides the fields

# This loop assesses whether a file that contains the QC summary exists, if not one can be created
if test ! -f QC_comparison_summary.txt ;
then
cat $QC_name_input | cut -f2 > QC_comparison_summary.txt
fi #Comparison of P,F,W


qfile_PorForW=$"Files_that_"$quality_control_loop"_""QC_"$qcount"_S1.txt" # create files for passing or failing the QC
if [ $quality_control_loop == PASS ] ; 
then 
echo "$fastqc_loop" >> $qfile_PorForW
fi # PASS

if [ $quality_control_loop == FAIL ] ; 
then 
echo "$fastqc_loop" >> $qfile_PorForW
fi # PASS

if [ $quality_control_loop == WARN ] ; # Creates a file containing sequences that are flaged as WARNING
then 
echo "$fastqc_loop" >> $qfile_PorForW
fi # PASS
# $QC_control


#loop check
echo $quality_control_loop
count=$((count + 1))

# Assess the number of PASS, FAIL, WARN 
qnumber_PorForW=$"Numbers_of_"$quality_control_loop"_""QC_"$qcount"_S1.txt"
number=$(cat $qfile_PorForW | wc -l) 
echo $number > $qnumber_PorForW

done # QC done
echo $count

done #ending of loop for zip loop list


#Create a report that summarises the numbers 

# Using the field QC plus the compC, it will keep track of the Quality Control Variable and their respective sequence grades
compC=$0 

#Loops the list of quality control checks
while read QC_list;
do
compC=$((compC + 1))
echo $QC_list > temp_QC_F.txt #Ensures the Quality Variables will be pasted along side, PASS, FAIL, WARN fields
#Re-arranging the name, for input into paste function
QC_N=$"QC_"$compC"_"
QC_F=$"_C1.txt"

#Matches the number of PASS, FAIL, WARN with quality check (based on compC)
F_N=$(ls -1 Number*.txt | grep $QC_N | grep "FAIL") 
P_N=$(ls -1 Number*.txt | grep $QC_N | grep "PASS")
W_N=$(ls -1 Number*.txt | grep $QC_N | grep "WARN")

#Combines the fields to produce a comparison of the grades
paste temp_QC_F.txt $P_N $F_N $W_N >> Comparison_Table_Unheaded.txt

done < QC_comparison_summary.txt # Ends the loop that goes through quality check variables
# The result of the paste, is unheaded, thus using echo, we can add a series of headers, seperated by tab delimited
# They are added to a file first, followed by the Comparison_Table_Unheaded.txt. We use the temporary as a intermediate to change its name to Comparison_Table.txt
(echo -e "Quality_Check\tPass\tFail\tWarn" && cat Comparison_Table_Unheaded.txt) > file1.txt && mv file1.txt Comparison_Table.txt

rm -f -R *.fastqc

#------------------------------------------------Stage2: Align Read pairs to Trypanosoma congolense genome 

# We achieve this using bowtie2

## Unzip the Genome file : Using gunzip, we can unzip the file to a fasta, -c ensures that the source file is left unchanged. 

gunzip -c /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz >Tcongo_genome.fasta

## We will index this : This enables easier placement of read sequences against the reference genome, matching at positions in which the reads are only slightly missmatched or peferecly matched
# This reduces the amount of time needed to align the reads to multiple locations across the T. congolense genome
# Genomes are only needed to be indexed a single time, regardless of the quantity of samples needed to be mapped.
bowtie2-build Tcongo_genome.fasta Tcongo # Tcongo is an argument

## We need to map the reads


# For every line in ID of 100k file, make a variable of ID + the additional tag
# The 100k.fqfiles contains alot of usefull information on the reads generated from the RNA sequencing. 
# Cat writes this index file to the pipeline, we use awk to print the first field, ID, only show lines from the second forward as we need to get rid of the header, then the output is put in a new file 

cat /localdisk/data/BPSM/AY21/fastq/100k.fqfiles | awk '{print $1}' | tail -n +2  > Tcongo_ID.txt

#Before we map the the reads, we must first assess whether the reads we are going to map are either paired or unpaired.
# Though looking 100k.fqfiles we can see that there is only two samples, if larger set of data is given, it is important to take into account unpaired reads as well. As they are inputed differently in the bowtie2 command line.
#In order to determine where it is paired or unpaired we can count the number of files, that share the sample sample name. If 2 then its paired, if 1 then unpaired. 
# The pipeline below makes it presentable into the bowtie2 command line as m1 and m2 for paired or as -U for unpaired.


# For every line in ID of 100k file, make a variable of ID (The loop will go through a list of these variables)
# The 100k.fqfiles contains alot of usefull information on the reads generated from the RNA sequencing. 
# Cat writes this index file to the pipeline, we use awk to print the first field, ID, only show lines from the second forward as we need to get rid of the header, then the output is put in a new file 

for reads in $(cat /localdisk/data/BPSM/AY21/fastq/100k.fqfiles | awk '{print $1}' | tail -n +2 ) ; #Loop ID list

do

#First we must assess whether it is paired or unpaired (filters for selected sample, then counts the number of read files)
Pair_Check=$(ls -l /localdisk/data/BPSM/AY21/fastq/*.fq.gz | grep -c $reads)

#The files will be outputed as -S .sam 
#Output variable to create distinct sam files dependent on the variables in the loop (ID)
read_sam=$"100k."$reads".sam" 


if [ $Pair_Check == 2 ] ; # If the number of files is 2 then execute paired bowtie2
then
#Create read variables to input into bowtie2 command line for paired mapping
read1=$"/localdisk/data/BPSM/AY21/fastq/100k."$reads"_1.fq.gz" 
read2=$"/localdisk/data/BPSM/AY21/fastq/100k."$reads"_2.fq.gz" 
#-x Tcongo, refers to the index argument we declared when we built the reference genome
fi
bowtie2 -x Tcongo -1 $read1 -2 $read2 -p 10 -S $read_sam # -p refers to threads, this process can be very demmanding, thus allocating additional threads may be beneficial 
if [ $Pair_Check == 1 ] ;
then 
 # Otherwise, execute unpaired bowtie2 mapping
readU="/localdisk/data/BPSM/AY21/fastq/100k."$reads"_U.fq.gz" #I assume that the the unpaired read would be declared in some manner, here I declare it as U
bowtie2 -x Tcongo -U $readU -p 10 -S $read_sam 
fi

done # End ID list lop

#The output file can take up a large amount of storage depending on the size and number of reads needed to be mapped
#Thus, for future reference, if the data being dealt with is extremely large, it is recommended  to view the data using grep or head as the file viewers will crash



## Converting the sam files to bam (http://quinlanlab.org/tutorials/samtools/samtools.html)
# To obtain data that can be analysed from the alignment, it has to be converted into a binary equivalent 
# We do this to make it easier for programes to compute. This however, comes at a cost that the files are difficult for us to understand at a glance
# The convertion of sam files to bam is carried out using Samtools

for sam_ID in $(ls . | grep ".sam"); # Loop of all sam files 
 do
sam_b=${sam_ID::-4}".bam" #Removing .sam to replace with .bam ::-4 removes 4 characters from the end 

echo $sam_ID is being converted into .bam #Notification that the converstion is taking place

# samtool view is the command used to convert 
# -S is used to specify .sam input as the command line is defaulted to expect .bam input
# It is also neccesary to state the output as .bam using -b
# -@ is a sematical key that adds threads when declared. As the process takes some time, we have alocated 10 
samtools view -@ 10 -S -b $sam_ID > $sam_b
done 
echo "Conversion to Bam complete" 

## Samtools sort 
#The alignment of the FASTQ read sequences to the genome are generated in randomized order, but in relation to the reference T. congolense reference genome
#The SAM and the resulting BAM file remains the order in which the reads are inputed in the FASTQ files. 
#In order to carry out analysis on the output data, it is neccesary to maniplute the .bam file further. 
#This time ordered in a manner that the reads alignments match the locations on the reference genome
#This process is carried out using samtools sort 

for bam_ID in $(ls . | grep ".bam$"); # loops using a list of .bam files : as we need to sort them
do
bam_b=${bam_ID::-4}"_sorted.bam" #Naming the sorted .bam output files 

echo "$bam_ID is being ordered"

#-o states the output, with -@ stating 10 threads (makes the process faster) 
samtools sort $bam_ID -o $bam_b -@ 10
done # ends .bam loop 
echo "Ordering of .bam files complete" 

## Indexing the sorted.bam files
#The genome being indexed enables us to rapidly retrieve areas of overlapping alignments found between T. congolense and our reads
#Additional benefit to indexing and neccesity, is in cases you wish to use a genome viewer for more efficient visualization, indexing is a requirment
#We use samtools index 

for sorted_bam in $(ls . | grep ".sorted.bam$"); #loops for every sorted bam file in directory 
do
echo "Indexing $sorted_bam" 
samtools index $sorted_bam #Indexing it
echo "Indexing complete" 
done # ends sorted bam loop

#--------------------------------------------------Stage3: Generating count data

## Generating count data: the number of reads that align to the regions of the genome that code for genes; bed tools
# bedtools multicov was used to quantify the alighments from the bam files we have generated, reporting the regions of alignment overlap 
# As the bai file generated is an index file, we will not use it in the following step. Instead we will use the sorted bam file
# The bam file is what contains all the sequence alignments. The bai file acts as a associate, containing similar suffix to the bam file
# Acting as a reference table, enabling the computer programe to locate specific locations of the bam files without the need to investigate the complete sequences
# Thus, bai file is useless without the bam file

#First we want to generate the count data for each replicate of each sample in relation to the T. congolense reference genome 
 
for count_data in $(ls *sorted.bam); # Loop all sorted.bam files 
 do
count_data_output=${count_data::-11}".csv" # Remove _sorted.bam to be replaced with .csv 

echo "Counting number of overlaps between $count_data and T. congolense at specific locals"

# Inputing the sorted bam files into the bedtool multicov command line following -bams  in reference to -bed T. congolense 
bedtools multicov  -bams $count_data -bed  /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed > $count_data_output
done # end of count_data loop

echo "Counting complete"

## Now the mean of count per gene (expression levels) was calculated using bedtools 

# First we need to average the replicates of each other (sort == -k4,4)

#Exports 100k.fqfiles to pipeline, show all lines after headers, 1st sorted by time, then replicate and then treatment, extracting the the time list 
#Then awk is used to create files containing the replicates of that specifc sample. 3 files of number counts per gene 
cat /localdisk/data/BPSM/AY21/fastq/100k.fqfiles | tail -n +2 |sort -k4,4 -k2,2 -k3,3 | awk '{
FS="\t"; Grouping_of_names=$2"_"$4"_"$5"_G.txt"; 
printf "%s.%s.%s\n", "100k", $1, "csv" >> Grouping_of_names;}' 


#Now we must compare means of each individual list inside a nested loop

ls . | grep ".G.txt$" > sample_grouped_files.txt

#out loop loops the ordered and group file

while read sample_grouped; 
do
SG_E=${sample_grouped::-4}
cat $sample_grouped > SG_inner.txt


#inner loop loops the files inside of the ordered grouped file, allowing for their last field (count per gene) and makes it easier to mean
count=0
while read SG_I;
do
SG_I_name=${SG_I::-4} ; 

cat $SG_I | awk 'FS="\t" {print $NF}' > $SG_E"_"$s_g_name"_"$count"o.txt"
count=$((count +1 ))
done <SG_inner.txt

###Combines reference Bed file from database with last field (count per gene) of each csv file within each group,

paste $SG_E*o.txt > $SG_E"m1.txt"

# We need to then calculate the mean of the last 3 fields in the new 

mean_output=$SG_E"_mean.txt"
 
 
field_value=$(awk -F '\t' '{print NF; exit}' $SG_E"m1.txt")

awk -F '\t' -v FN=$field_value -v output_value=$mean_output '{
OFS="\t"; total =0; 
for (count=1; count<=FN;count++) total=total+$count; mean=total/FN ; print mean;
}' $SG_E"m1.txt" > $mean_output

#Creates header for mean field then conjugates them together, to be pasted below
(echo -e "$SG_E" && cat $mean_output) > temp_file.txt && mv  temp_file.txt $SG_E"_HMean.csv"

#This adds gene ID and Gene descriptors to their respective columns
(echo -e "Gene_ID\tGene_Descriptor" && (cat /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed | awk -F '\t' '{OFS="\t" ; print $4,$5}'))  > temp_file.txt && mv  temp_file.txt Gene_Descriptor_Field.csv

#This creates a table with ID, Descripter and mean
paste Gene_Descriptor_Field.csv $SG_E"_HMean.csv" > $SG_E"_Complete_TRUE.csv"


done <sample_grouped_files.txt 


#There are two methods of calculating "fold change" Find out the difference 1. (Y-X)/X 



# Creating "fold change" data #
#Creat list to keep track of variables
ls . | grep "mean.txt$" > variable_files.txt
#Over time comparison realtive to WT controls 
 
source Fold_change_PairWise.sh

source Fold_change_UI.sh

#


