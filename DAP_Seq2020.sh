#/bin/bash

#######################
# Handle Command Line #
#######################
# this is used for Bowtie to make sure I use enough CPUs
threads=$1
if [ "$1" == "" ]
then
        threads=5
fi


#######################
# Modified 07.20.2017 #
#######################

# I wanted to try different folder names so I didn't have to do so much organizing at the end

#####################
# Setup Environment #
#####################
genomesPath=/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/GENOMES/
rawdataPath=$2
processingPath=$3
processedPath=$4
scriptsPath=/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/scripts
gemJar=gem.v2.5.jar

echo "Hello Ying! Today is $(date). I'm excited to analyze your data for you!"
echo "My name is New_revolution and I was made to process DAP-Seq data"
echo "Your raw data came from $rawdataPath"
echo "Your genomes came from $genomesPath"
echo "Your output will go to $processingPath"


################
# Load Modules #
################
module load trim_galore
module load Bowtie/2.3.0
module load SAMtools/1.9
module load Java/11
module load Python/2.7.13
module load R/3.4.2
# module load picard/2.18.3


# Plant to genome map
# use for analysis
declare -A plant_map
plant_map['At']="Arabidopsis_TAIRv10"
plant_map['Cr']="Crubella_JGIv1.0"
plant_map['EsJ']="Esalsugineum_JGIv1.0"
plant_map['Si']="Sirio_CoGev2.51"
plant_map['Sp']="Sparvulav_2.1"

#At.chr27445.gtf 
#Es.chr26351.gtf 
#Sp.chr26847.gtf
#Si.chr27367.gtf 
#Cr.chr26521.gtf

# Plant to genome map
# no chr
# use for paper 
#declare -A plant_map
#plant_map['At']="Arabidopsis_TAIRv10paper"
#plant_map['Cr']="Crubella_JGIv1.0paper"
#plant_map['EsJ']="Esalsugineum_JGIv1.0paper"
#plant_map['Si']="Sirio_CoGev2.51paper"
#plant_map['Sp']="Sparvulav_2.1paper"

#At.27445.gtf
#Es.26351.gtf
#Sp.26847.gtf
#Si.27367.gtf
#Cr.26521.gtf


# Plant to chromosome file
declare -A chromosome_file
chromosome_file['At']="At_genome.chrom.sizes"
chromosome_file['Cr']="Cr_genome.chrom.sizes"
chromosome_file['EsJ']="EsJ_genome.chrom.sizes"
chromosome_file['Si']="Si_genome.chrom.sizes"
chromosome_file['Sp']="Sp_genome.chrom.sizes"

# Plant to genome_size
# These numbers were from counts.list file generated from Bowtie_build_genome.txt then I used excel to sum up the second column
declare -A genome_size
genome_size['At']="119146348"
genome_size['Br']="284855373"
genome_size['Cr']="134834574"
genome_size['EsJ']="243117811"
genome_size['EsC']="233653057"
genome_size['Si']="259494581"
genome_size['Sp']="120016289"

# Plant to madebam_file
# This should give me all the bam_files generated
declare -A madebam_file
madebam_file['At']="bam_files_At"
madebam_file['Cr']="bam_files_Cr"
madebam_file['EsJ']="bam_files_EsJ"
madebam_file['Si']="bam_files_Si"
madebam_file['Sp']="bam_files_Sp"

# Plant to madenarrow_file
# This should give me all the bam_files generated
declare -A madenarrow_file
madenarrow_file['At']="narrowpeak_GEM_At"
madenarrow_file['Cr']="narrowpeak_GEM_Cr"
madenarrow_file['EsJ']="narrowpeak_GEM_EsJ"
madenarrow_file['Si']="narrowpeak_GEM_Si"
madenarrow_file['Sp']="narrowpeak_GEM_Sp"

# Plant to madesorted_file
# This should give me all the bam_files generated
declare -A madesorted_file
madesorted_file['At']="sorted_narrowpeak_GEM_At"
madesorted_file['Cr']="sorted_narrowpeak_GEM_Cr"
madesorted_file['EsJ']="sorted_narrowpeak_GEM_EsJ"
madesorted_file['Si']="sorted_narrowpeak_GEM_Si"
madesorted_file['Sp']="sorted_narrowpeak_GEM_Sp"


# Move to the processing folder
pushd $processingPath >/dev/null 2>&1

# Get Plants
for plant in $(ls)
do
        plantIndexPath="$genomesPath/${plant_map[$plant]}/index/${plant_map[$plant]}"
        readdistributionPath="$genomesPath/Read_distribution_txt/${chromosome_file[$plant]}"
        plantGenomePath="$genomesPath/${plant_map[$plant]}/"
        plantbamfilePath="$processedPath/${madebam_file[$plant]}/"
        plantnarrowPath="$processedPath/${madenarrow_file[$plant]}/"
        plantsortednarrowPath="$processedPath/${madesorted_file[$plant]}/"

        cd $plant
        echo $plant
        for protein in $(ls)
        do
                cd $protein
                for tissue in $(ls)
                do
                        cd $tissue
                        for samplenum in $(ls)
                        do
                                cd $samplenum/R1
                                pwd
                                r1file=${rawdataPath}/${plant}-${protein}-${tissue}_${samplenum}_L001_R1_001.fastq
                                r1trimmedfile=trimmed/${plant}-${protein}-${tissue}_${samplenum}_L001_R1_001_trimmed.fq
                                mkdir trimmed >/dev/null 2>&1
                                echo "Running trim_galore on ${plant_map[$plant]}, $protein, $tissue, $samplenum"
                                trim_galore --three_prime_clip_R1 1 --fastqc --output_dir trimmed $r1file
                                echo "Running bowtie2 on <stuff> with ${plant_map[$plant]} index"
                                echo "Genome used is ${plantIndexPath}"
                                bowtie2 --threads $threads --very-sensitive -x ${plantIndexPath} -q -U $r1trimmedfile -S ${plant}-${protein}-${tissue}_${samplenum}_aligned.sam
                                echo "Creating Filtered copy of aligned SAM without XS tag"
                                grep -v "XS:" ${plant}-${protein}-${tissue}_${samplenum}_aligned.sam > ${plant}-${protein}-${tissue}_${samplenum}_aligned_no_multi.sam
                                echo "Sorting the aligned and XS filtered SAM"
                                samtools sort --threads $threads -O BAM -o ${plant}-${protein}-${tissue}_${samplenum}_aligned_sorted_no_multi.bam ${plant}-${protein}-${tissue}_${samplenum}_aligned_no_multi.sam
                                #echo "Remove duplicates with PICARD"
                                #java -jar $PICARD MarkDuplicates I=${plant}-${protein}-${tissue}_${samplenum}_aligned_sorted_no_multi.bam O=${plant}-${protein}-${tissue}_${samplenum}_aligned_sorted_no_multi_picard.bam METRICS_FILE=PICARD_${plant}-${protein}-${tissue}_${samplenum}.txt
                                echo "Running gem on the aligned, sorted, and filtered sample data"
								echo "This is the chromosome text file ${chromosome_file[$plant]}"
                                echo "This is the genome $plantGenomePath"
                                echo "This is the size of the chromosome ${genome_size[$plant]}"
                            	java -Xmx10G -jar $scriptsPath/$gemJar --d $scriptsPath/Read_Distribution_default.txt --s ${genome_size[$plant]} --g ${readdistributionPath}.txt --genome $plantGenomePath --expt ${plant}-${protein}-${tissue}_${samplenum}_aligned_no_multi.sam --f SAM --k_min 6 --k_max 20 --k_seqs 600 --k_neg_dinu_shuffle --t $threads --outNP
                                samtools index ${plant}-${protein}-${tissue}_${samplenum}_aligned_sorted_no_multi.bam
                                echo "copy .bam to bam_file"
                                #change bam file directory
                                mkdir -p $plantbamfilePath
                                cp *.bam $plantbamfilePath
                                cp *.bai $plantbamfilePath
                                echo "re-name .narrowpeak to include plant/protein/tissue etc."
                                #Change narrow peak directory
                                cp ./out/out_GEM_events.narrowPeak ./out/${plant}-${protein}-${tissue}_${samplenum}_out_GEM_events.narrowPeak
                                cp ./out/out_GPS_events.narrowPeak ./out/${plant}-${protein}-${tissue}_${samplenum}_out_GPS_events.narrowPeak
                                echo "copy .narrowpeak to narrow_file"
                                mkdir -p $plantnarrowPath
                                cp ./out/${plant}-${protein}-${tissue}_${samplenum}_out_GEM_events.narrowPeak $plantnarrowPath
                                echo "sorting .narrowpeak by chromosome position"
                                sort -k 1,1 -k 2,2n ./out/${plant}-${protein}-${tissue}_${samplenum}_out_GEM_events.narrowPeak > ./out/sorted_${plant}-${protein}-${tissue}_${samplenum}_out_GEM_events.narrowPeak 
                                echo "copy sorted to sorted_narrow_file"
                                mkdir -p $plantsortednarrowPath
                                cp ./out/sorted_${plant}-${protein}-${tissue}_${samplenum}_out_GEM_events.narrowPeak $plantsortednarrowPath                         
                                cd ../..
                        done
                        cd ..
                done
                cd ..
        done
        echo "Combine all sam files for all proteins in this plant"
        echo "Run gem, get peaks on combined plant data"
        cd ..
done

popd >/dev/null 2>$1


echo "All done Ying, good luck with the rest of the analysis!"
