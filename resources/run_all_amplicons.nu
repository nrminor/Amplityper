#!/usr/bin/env nu

# construct amplicon file
def define_amplicons [bed: string] -> string {
    [$bed]

    open $bed \
    | from tsv -n | get column4 \
    | str replace "_LEFT" "" \| str replace "_RIGHT" "" \
    | uniq | save amplicons.txt

    return "amplicons.txt"

}

# script entry point
def main [
    --bed: string,
    --data: string,
    --min_reads: int
    ] {
        [$bed, $data, $min_reads]

        # check if the script has been run before
        if ["./remaining_amplicons.txt" | path exists] {

            echo "Picking up from previous run."
            let text_file = "./remaining_amplicons.txt"

            # lock the remaining amplicons at the start of the run to avoid confusing the while loop
            open $text_file | save "remaining.lock"

        } else {

            echo "No previous runs detected."
            echo "Starting from first amplicon."

            let text_file = define_amplicons $bed
            open $text_file | save "remaining.lock"
        }

        # double check that the provided directory actually exists
        if ![$data | path exists] {
            echo "Error: FASTQ directory '$data' does not exist." ; exit 1
        }

        # if it does exist, set fastq directory
        let fastq_dir = $data

        # Loop through each line in the text file
        open "remaining.lock" | each { |amplicon|

            echo "Running on amplicon: $amplicon"

            # Run a command with the line as an argument
            (
                nextflow run .
                --fastq_dir $fastq_dir
                --desired_amplicon $amplicon
                --min_reads $min_reads
                --cleanup true
                -profile "local"
                -resume
            )

            # log completed amplicons
            touch "complete_log.txt"
            $amplicon | save --append complete_log.txt

        }

        # remove now unnecessary lock
        rm remaining.lock
}
