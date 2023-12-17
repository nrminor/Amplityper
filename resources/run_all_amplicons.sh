#!/bin/bash

# Check if a text file and fastq path are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <amplicon_bed_file> <path_to_data_dir>"
    exit 1
fi

# construct amplicon file
cut -f 4 $1 | sed 's/_LEFT//g' | sed 's/_RIGHT//g' | sort -u | uniq -u  > amplicons.txt

# check if the script has been run before
if [ -f "remaining_amplicons.txt" ]; then
    echo "Picking up from previous run."
    text_file="remaining_amplicons.txt"
else
    echo "No previous runs detected."
    echo "Starting from first amplicon."
    text_file=amplicons.txt
fi

# set fastq directory
fastq_dir="$2"

# Function to handle the interrupt signal
handle_interrupt() {
    echo "Interrupt received. Exiting the loop..."
    exit 1
}

# Trap the SIGINT signal (Ctrl+C)
trap handle_interrupt SIGINT

# double check that the provided directory actually exists
if [ ! -d "$fastq_dir" ]; then
    echo "Error: FASTQ directory '$fastq_dir' does not exist."
    exit 1
fi

# lock the remaining amplicons at the start of the run to avoid
# confusing the while loop
cat "$text_file" > remaining.lock

# Loop through each line in the text file
while IFS= read -r amplicon
do
    echo "Running on amplicon: $amplicon"

    # Run a command with the line as an argument
    # Replace 'your_command' with the actual command you want to run
    nextflow run . \
    --fastq_dir "$fastq_dir" \
    --desired_amplicon "$amplicon" \
    -profile local \
    -resume

    # Check if the interrupt signal has been received
    if [[ "$?" -ne 0 ]]; then
        echo "Command interrupted."
        break
    fi

    # log completed amplicons
    touch "complete_log.txt"
    echo $amplicon >> complete_log.txt

    # figure out which amplicons remain for resumability
    touch "remaining_amplicons.txt"
    grep -vxFf "complete_log.txt" "amplicons.txt" > "remaining_amplicons.txt"

    # Example: echo the line
    # echo "$line"
done < remaining.lock

# remove now unnecessary lock
rm remaining.lock
