#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
from argparse import ArgumentParser, ArgumentTypeError
from subprocess import Popen, PIPE, run as subprocess_run
import gzip



### METHOD DEFINITIONS ###
### ------------------------------------------------------------------------------ ###

# define function that parses boolean command line argument strings into true booleans
def str2bool(v):
    """
    Standard function used across many D.A. Baker scripts for parsing command
    line argument strings into true boolean objects
    """
    
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')

# define a function that parses all command line arguments for this module
def parse_process_arrays_args():
    """
    Parses the python script arguments from bash and makes sure files/inputs 
    are valid
    """
    parser = ArgumentParser()
    
    parser.add_argument('--primer_fasta',
                        type=str,
                        help='where fasta primer, headers must be nomenclature: \n\
                        <primer_name>_<primer_group>_<primer_direction>_<optional_additional_primer>\n\
                        or\n\
                        <primer_name>_<primer_group>_<primer_direction>',
                        required=True)
    parser.add_argument('--in_path',
                        type=str,
                        help='input path to forward reads.',
                        required=True)
    parser.add_argument('--in2_path',
                        type=str,
                        help='input path to reverse-complement paired reads.',
                        required=False)
    parser.add_argument('--outm',
                        type=str,
                        help='output path to output_path forwared reads',
                        required=True)
    parser.add_argument('--outm2',
                        type=str,
                        help='output path to reverse-complement paired reads',
                        required=False)
    parser.add_argument('--temp_dir',
                        type=str,
                        help='path to intermediate files',
                        default=None,
                        required=False)
    parser.add_argument('--mem',
                        type=int,
                        help='ram to allocate in MB',
                        default=8000,
                        required=False)
    parser.add_argument('--bbmap_dir',
                        type=str,
                        help='dir where bbmap tools are located',
                        default='/bbmap',
                        required=False)
    parser.add_argument('--ktrim',
                        type=str,
                        help='set to f if you want to leave the primers, t to trim the primers',
                        default='f',
                        required=False)
    
    return parser.parse_args()

# function that concatenates fastq.gz in script. NOTE: this block will be removed in future implementations of this script.
def concat_fastq_gz(input_fastq_filepath, out_file):
    """
    This method loops through each FASTQ.gz file and unzips it into a new FASTQ. This will likely
    be removed in future versions to keep file I/O more lightweight.
    """
    
    fastq_sequences = FastqGeneralIterator(gzip.open(input_fastq_filepath, "rt"))
    for (f_id, f_seq, f_q) in fastq_sequences:
        out_file.write('@{0}\n{1}\n+\n{2}\n'.format(f_id, f_seq, f_q))

# function that uses the primer file to define which primers are part of which amplicons, and outputs a dictionary 
def fasta_groups_to_dict(fasta_path):
    """
    This method assigns a group number to each primer pair
    """

    primer_fasta_dict = {}
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    for fasta in fasta_sequences:
        header_name, sequence = fasta.id, str(fasta.seq)
        header_split = header_name.split('_')
        # print(header_split)
        # primer_name = header_split[0]
        primer_group = header_split[1]
        primer_direction = header_split[2]
        if primer_group not in primer_fasta_dict.keys():
            primer_fasta_dict[primer_group] = {}
        #print(primer_fasta_dict)
        if primer_direction not in primer_fasta_dict[primer_group]:
            primer_fasta_dict[primer_group][primer_direction] = []
        fastq_list = primer_fasta_dict[primer_group][primer_direction]
        fastq_list.append((header_name, sequence))
        primer_fasta_dict[primer_group][primer_direction] = fastq_list
    return primer_fasta_dict

# function that splits primer sequences into separate, per-amplicon FASTA files
def dict_to_grouped_fastas(grouped_fasta_dir, grouped_fasta_dict):
    """
    This method sorts the primers by amplicon group
    """
    
    primer_group_path_dict ={}
    for primer_group, primer_direction_dict in grouped_fasta_dict.items():
        for primer_direction, fasta_data_list in primer_direction_dict.items():
            grouped_ref_path = os.path.join(grouped_fasta_dir, '{0}_{1}.fasta'.format(primer_group, primer_direction))
            if primer_group not in primer_group_path_dict.keys():
                primer_group_path_dict[primer_group] = {}
            primer_group_path_dict[primer_group][primer_direction] = grouped_ref_path
            with open(grouped_ref_path, 'w') as out_file:
                for fasta_data in fasta_data_list:
                    out_file.write('>{0}\n'.format(fasta_data[0]))
                    out_file.write('{0}\n'.format(fasta_data[1]))

    return primer_group_path_dict

### ------------------------------------------------------------------------------ ###



### RUNNING THE WORKFLOW ###
### ------------------------------------------------------------------------------ ###

def main():
    # get the arguments
    args = parse_process_arrays_args()

    # save arguments to global variables by same names
    primer_fasta = args.primer_fasta
    in_path = args.in_path
    in2_path = args.in2_path
    outm = args.outm
    outm2 = args.outm2
    temp_dir = args.temp_dir
    mem = args.mem
    bbmap_dir = args.bbmap_dir
    ktrim = args.ktrim
    if temp_dir is None:
        temp_dir = os.path.join(os.path.dirname(outm), 'TMP')

    # create the required directories for outputting files
    os.makedirs(temp_dir, exist_ok=True)
    grouped_fasta_dir = os.path.join(temp_dir,'single_fasta')
    os.makedirs(grouped_fasta_dir, exist_ok=True)

    # Create grouped fasta files for the forward and reverse primer pairs.
    # create an iterable dictionary of paths.
    grouped_fasta_dict = fasta_groups_to_dict(primer_fasta)
    primer_group_path_dict = dict_to_grouped_fastas(grouped_fasta_dir, grouped_fasta_dict)
    print('Primer Group Count:', len(primer_group_path_dict.keys()))

    left_outpath_list = []
    right_outpath_list = []
    primer_group_count = len(primer_group_path_dict.keys())
    primer_group_count = len(primer_group_path_dict.keys())
    i = 1
    for primer_group, primer_dict in primer_group_path_dict.items():
        print('Processing Primer Group {0} of {1} -- name:{2}'.format(i, primer_group_count, primer_group))
        # create intermediate/temporary file names to be used by bbduk and repair sh to conduct the primer matching
        i += 1
        in_fasta_left_path = primer_dict['LEFT']
        base_left_fasta = os.path.basename(in_fasta_left_path)[:-6]
        in_fasta_right_path = primer_dict['RIGHT']
        base_right_fasta = os.path.basename(in_fasta_right_path)[:-6]
        fasta_left_temp_1 = os.path.join(temp_dir, '{}_temp_1.fastq.gz'.format(base_left_fasta))
        fasta_left_temp_2 = os.path.join(temp_dir, '{}_temp_2.fastq.gz'.format(base_left_fasta))
        fasta_left_temp_3 = os.path.join(temp_dir, '{}_temp3.fastq.gz'.format(base_left_fasta))
        left_outpath_list.append(fasta_left_temp_3)
        fasta_right_temp_1 = os.path.join(temp_dir, '{}_temp_1.fastq.gz'.format(base_right_fasta))
        fasta_right_temp_2 = os.path.join(temp_dir, '{}_temp_2.fastq.gz'.format(base_right_fasta))
        fasta_right_temp_3 = os.path.join(temp_dir, '{}_temp3.fastq.gz'.format(base_right_fasta))
        right_outpath_list.append(fasta_right_temp_3)
        # run bbduk on the forward primer set of the group
        print('Running bbduk left: ', primer_group)
        ssh_cmd_out = subprocess_run(['java', '-ea',
                                    '-Xmx{0}m'.format(mem),
                                    '-Xms{0}m'.format(mem),
                                    '-cp', '{0}/current/'.format(bbmap_dir),
                                    'jgi.BBDuk',
                                    'in={}'.format(in_path),
                                    'in2={}'.format(in2_path),
                                    'outm={}'.format(fasta_left_temp_1),
                                    'outm2={}'.format(fasta_right_temp_1),
                                    'ref={0}'.format(in_fasta_left_path),
                                    'ktrim={}'.format(ktrim),
                                    'k=21',
                                    'qtrim=f',
                                    'trimq=30',
                                    'hdist=3',
                                    'rcomp=f',
                                    'minlength=75',
                                    'restrictleft=32',
                                    'requireBothBad=f',
                                    'overwrite=t'],
                                    shell=False,
                                    stdout=PIPE,
                                    stderr=PIPE)
        # import os
        os.write(1, ssh_cmd_out.stderr)
        # run repair again as the threads can put the reads out of order and will be missing in bbduk
        # print('Running repair left: ', primer_group)
        # ssh_cmd_out = subprocess_run(['{0}/repair.sh'.format(bbmap_dir),
        #                               'in={}'.format(fasta_left_temp_1),
        #                               'in2={}'.format(in2_path),
        #                               'out={}'.format(fasta_left_temp_2),
        #                               'out2={}'.format(fasta_right_temp_1)],
        #                              shell=False,
        #                              stdout=PIPE,
        #                              stderr=PIPE)
        # print(ssh_cmd_out)
        # Run bbduk on the reversecomplement matched pair primer of the group
        print('Running bbduk right: ', primer_group)
        ssh_cmd_out = subprocess_run(['java', '-ea',
                                    '-Xmx{0}m'.format(mem),
                                    '-Xms{0}m'.format(mem),
                                    '-cp', '{0}/current/'.format(bbmap_dir),
                                    'jgi.BBDuk',
                                    'in={}'.format(fasta_left_temp_1),
                                    'in2={}'.format(fasta_right_temp_1),
                                    'outm={}'.format(fasta_left_temp_3),
                                    'outm2={}'.format(fasta_right_temp_3),
                                    'ref={0}'.format(in_fasta_right_path),
                                    'ktrim={}'.format(ktrim),
                                    'k=21',
                                    'qtrim=f',
                                    'trimq=30',
                                    'hdist=3',
                                    'rcomp=f',
                                    'minlength=75',
                                    'restrictleft=32',
                                    'requireBothBad=f',
                                    'overwrite=t'],
                                    shell=False,
                                    stdout=PIPE,
                                    stderr=PIPE)
        os.write(1, ssh_cmd_out.stderr)
        # run repair again as the threads can put the reads out of order and will be missing in bbduk
        # print('Running repair right: ', primer_group)
        # ssh_cmd_out = subprocess_run(['{0}/repair.sh'.format(bbmap_dir),
        #                               'in={}'.format(fasta_left_temp_2),
        #                               'in2={}'.format(fasta_right_temp_2),
        #                               'out={}'.format(fasta_left_temp_3),
        #                               'out2={}'.format(fasta_right_temp_3)],
        #                              shell=False,
        #                              stdout=PIPE,
        #                              stderr=PIPE)
        # print(ssh_cmd_out)

    # concatenate the  primer grouped fastq files into a single file for each direction
    print('Concatenating files to create a single file')
    with gzip.open(outm, 'wt') as out_file:
        for left_outpath in left_outpath_list:
            concat_fastq_gz(left_outpath, out_file)
            print(left_outpath)

    with gzip.open(outm2, 'wt') as out_file:
        for right_outpath in right_outpath_list:
            concat_fastq_gz(right_outpath, out_file)
            print(right_outpath)

if __name__ == "__main__":
    main()

### ------------------------------------------------------------------------------ ###
