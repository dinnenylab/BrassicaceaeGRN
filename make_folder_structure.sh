#!/bin/bash

################
# Make Folders #
################

#######################
# Handle Command Line #
#######################
rawdataPath=$1
folderPath=$2
if [[ "$1" == "" || "$2" == "" ]]
then
	echo "Usage: $0 <fastq_raw_data_folder > <path_to_create_folders>"
	exit 1
fi

# Loop through each .fastq file in the rawdata path
for fastq_file_path in `ls "$rawdataPath"/*.fastq`
do
	# Get file name
	fastq_file=`basename $fastq_file_path`
	plant=`echo $fastq_file | cut -d '_' -f 1 | cut -d '-' -f 1`
	protein=`echo $fastq_file | cut -d '_' -f 1 | cut -d '-' -f 2`
	tissue=`echo $fastq_file | cut -d '_' -f 1 | cut -d '-' -f 3-`
	samplenum=`echo $fastq_file | cut -d '_' -f 2`
	readnum=`echo $fastq_file | cut -d '_' -f 4`
	
	mkdir -p "${folderPath}/${plant}/${protein}/${tissue}/${samplenum}/${readnum}" >/dev/null 2>&1
done
