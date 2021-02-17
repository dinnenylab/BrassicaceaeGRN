#/bin/bash

#######################
# Handle Command Line #
#######################
threads=$1
if [ "$1" == "" ]
then
        threads=5
fi


#######################
# Modified 04.10.2018 #
#######################

# Running RNA_Seq analysis

#####################
# Setup Environment #
#####################
genomesPath=/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/GENOMES
#genomesPath=/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/GENOMES
rawdataPath=$1
processingPath=$2
processedPath=$3
scriptsPath=/Carnegie/DPB/Data/Shared/Labs/Dinneny/Private/Ysun/2020/scripts
#scriptsPath=/home/ysun/scratch/scripts
#gemJar=gem.v2.5.jar

echo "Hello Ying! Today is $(date). I'm excited to analyze your data for you!"
echo "My name is RNA_Seq_analyzer and I was made to process RNA-Seq data"
echo "Your raw data came from $rawdataPath"
echo "Your genomes came from $genomesPath"
echo "Your output will go to $processingPath"


################
# Load Modules #
################
module load trim_galore/0.4.2
module load Bowtie/2.3.0
module load SAMtools/1.3.1
module load Java/8
module load Python/2.7.11
module load R/3.2.5
module load HISAT2/2.2.0
module load StringTie/1.3.3b


# Plant to genome map
declare -A plant_map
plant_map['At']="Arabidopsis_TAIRv10"
plant_map['Br']="Brapa_BRADv1.5"
plant_map['Cr']="Crubella_JGIv1.0"
plant_map['EsJ']="Esalsugineum_JGIv1.0"
plant_map['EsC']="Esalsugineum_CASv1.0"
plant_map['Si']="Sirio_CoGev2.51"
plant_map['Sp']="Sparvulav_2.1"

# Plant to chromosome file
declare -A chromosome_file
chromosome_file['At']="At_genome.chrom.sizes"
chromosome_file['Br']="Br_genome.chrom.sizes"
chromosome_file['Cr']="Cr_genome.chrom.sizes"
chromosome_file['EsJ']="EsJGI_genome.chrom.sizes"
chromosome_file['EsC']="EsCAS_genome.chrom.sizes"
chromosome_file['Si']="SiYm2sgenome.chrom.sizes.txt"
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
madebam_file['Br']="bam_files_Br"
madebam_file['Cr']="bam_files_Cr"
madebam_file['EsJ']="bam_files_EsJ"
madebam_file['EsC']="bam_files_EsC"
madebam_file['Si']="bam_files_Si"
madebam_file['Sp']="bam_files_Sp"

# Plant to madenarrow_file
# This should give me all the bam_files generated
declare -A madestring_file
madestring_file['At']="stringtie_At"
madestring_file['Br']="stringtie_Br"
madestring_file['Cr']="stringtie_Cr"
madestring_file['EsJ']="stringtie_EsJ"
madestring_file['EsC']="stringtie_EsC"
madestring_file['Si']="stringtie_Si"
madestring_file['Sp']="stringtie_Sp"

# Plant to madesorted_file
# This should give me all the bam_files generated
declare -A madesorted_file
madesorted_file['At']="sorted_narrowpeak_GEM_At"
madesorted_file['Br']="sorted_narrowpeak_GEM_Br"
madesorted_file['Cr']="sorted_narrowpeak_GEM_Cr"
madesorted_file['EsJ']="sorted_narrowpeak_GEM_EsJ"
madesorted_file['EsC']="sorted_narrowpeak_GEM_EsC"
madesorted_file['Si']="sorted_narrowpeak_GEM_Si"
madesorted_file['Sp']="sorted_narrowpeak_GEM_Sp"


# Plant to x file
declare -A plant_gff

plant_gff['At']="At.chr27445.gtf"
plant_gff['Cr']="Cr.chr26521.gtf"
plant_gff['EsJ']="Es.chr26351.gtf"
plant_gff['Si']="Si.chr27367.gtf"
plant_gff['Sp']="Sp.chr26847.gtf"


#plant_gff['At']="Araport11_chr_onlyminusCM.27445geneModels.gtf"
#plant_gff['Cr']="Crubella_183.primary.26521geneModels.gtf"
#plant_gff['EsJ']="Esalsugineum_173_v1.primary.26351geneModels.gtf"
#plant_gff['Si']="SiYm2s.gtf"
#plant_gff['Sp']="Sparvula_genome_annotation_v2.0.26847geneModels.fixed.gtf"

# Move to the processing folder
pushd $processingPath >/dev/null 2>&1

# Get Plants
for plant in $(ls)
do
        plantIndexPath="$genomesPath/${plant_map[$plant]}/index/${plant_map[$plant]}"
        readdistributionPath="$genomesPath/Read_distribution_txt/${chromosome_file[$plant]}"
        plantGenomePath="$genomesPath/${plant_map[$plant]}/"
        plantbamfilePath="$processedPath/${madebam_file[$plant]}/"
        plantstringPath="$processedPath/${madestring_file[$plant]}/"
        plantsortednarrowPath="$processedPath/${madesorted_file[$plant]}/"

        cd $plant
        #echo $plant
        for tissue in $(ls)
        do
        		cd $tissue
                #echo $tissue
                for treatment in $(ls)
                do
                		cd $treatment
                        #echo $treatment
                        for timepoint in $(ls)
                        do
                    			cd $timepoint
                        		#echo $timepoint
                        		for replicate in $(ls)
                       			do
                       					cd $replicate
                        				#echo $timepoint
                        				for dateN in $(ls)
                       					do
                       							cd $dateN
                        						#echo $timepoint
                        						for barcode in $(ls)
                       							do
                       									cd $barcode
                        								#echo $timepoint
                        								for samplenum in $(ls)
                       									do
                                								cd $samplenum
                                								echo "now trimming your files" 
                                								pwd
                                								r1file=${rawdataPath}/${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_L001_R1_001.fastq
                                								r1trimmedfile=trimmed/${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_L001_R1_001_trimmed.fq
                                								mkdir trimmed >/dev/null 2>&1
                                								echo "Running trim_galore on ${plant_map[$plant]}, $tissue, $treatment, $timepoint, $replicate, $dateN, $barcode, $samplenum"
                                								trim_galore --three_prime_clip_R1 1 --fastqc --output_dir trimmed $r1file
                                								echo "Running HISAT2 on <stuff> with ${plant_map[$plant]} index"
                                								echo "Genome used is ${plantIndexPath}"
                                								hisat2 --dta -x ${plantIndexPath} -q -U $r1trimmedfile -S ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_aligned.sam
                                								echo "convert SAM to BAM"
                                								samtools view -bS ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_aligned.sam > ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_aligned.bam
                                								echo "sort BAM files"
                                								samtools sort ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_aligned.bam -o ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_sorted.bam
                                								echo "make a .bai file"
                                								samtools index ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_sorted.bam
                                								echo "running stringtie to calculate FPKM"
                                								#stringtie -e ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_sorted.bam -o ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_stringtie.gtf
                                								##stringtie -B -G ${plantGenomePath}${plant_gff[$plant]} ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_sorted.bam -o ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_stringtie.gff
                                								stringtie -p 4 -e -o ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}.gtf -G ${plantGenomePath}${plant_gff[$plant]} ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}_sorted.bam &>> ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}.log
                                								echo "moving bam files"
                                								mkdir -p $plantbamfilePath
                                								cp *.bam $plantbamfilePath
                                								cp *.bai $plantbamfilePath
                                								echo "moving stringtie files"
                                								mkdir -p $plantstringPath
                                								cp ${plant}-${tissue}-${treatment}-${timepoint}-${replicate}-${dateN}-${barcode}_${samplenum}.gtf $plantstringPath
                                								echo "done! on to the next sample"
                                								cd ..
                        								done
                                						cd ..
                                				done
                        						cd ..
                        				done
                        				cd ..
                        		done
                        		cd ..
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
