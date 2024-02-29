### Compare HCMV Transcriptomes 2- and 6- days post-infection ###

# Step 1 in READMe.md

import argparse
import sys
import os
import subprocess
import glob
from Bio import SeqIO

# Input arguments

def check_args(args=None):
    parser = argparse.ArgumentParser(description='Homework #4')
    parser.add_argument('-i', '--input',help='path to input files', required='True') #input file
    return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_args(sys.argv[1:])
inpath = args.input

## Step 2: Build Index ##
def build_index(inpath):
    cwd = os.getcwd() # save the current working directory
    os.system("mkdir index_dir") # make directory to store HCMV index
    os.chdir("index_dir") # change to directory to store HCMV index
    # Create index for HCMV
    # Retrieve the complete HCMV reference genome fasta
    index_command = "datasets download genome accession GCA_027926585.1 --include gff3,rna,cds,protein,genome,seq-report"
    os.system(index_command) # run on the command line
    os.system("unzip ncbi_dataset.zip") # unzipping HCMV reference genome zip file
    # use unzipped data to build the HCMV index
    build_command = "bowtie2-build ncbi_dataset/data/GCA_027926585.1/GCA_027926585.1_ASM2792658v1_genomic.fna HCMV"
    os.system(build_command) # build HCMV index
    os.chdir("{0}/{1}".format(cwd,inpath)) # change to directory with fastq files
    download_list = glob.glob("SRR*") # adding elements to list
    transcriptomes = [ x for x in download_list if "fastq" in x ] # including only fastq files
    transcriptomes.sort() # put fastq files in order for mapping step
    id_list = [] #initialize list of file ids
    for i in range(0, len(transcriptomes)):
        if i % 2 == 0: # makes sure there are no repeat transcriptomes, only correct pairs of transcriptomes will be mapped to the index
            id = transcriptomes[i].split("_") # retrieving the SRR id
            id = id[0] # save SRR id of fastq file
            id_list.append(id) # append id to id_list
            os.chdir(cwd) # change back to Pipeline Project directory
            # mapping command
            map_command = "bowtie2 --quiet -x {0}/index_dir/HCMV -1 {1}/{2} -2 {1}/{3} -S {4}.sam --al-conc-gz {4}_mapped_%.fq.gz".format(cwd,inpath,transcriptomes[i],transcriptomes[i+1],id)
            os.system(map_command) # Calculate the number of read pairs before and after bowtie filtering
            os.chdir("{0}/{1}".format(cwd,inpath)) # change to directory with transcriptome fastq files
            before_command = "wc -l {}".format(transcriptomes[i]) #counting the lines before bowtie
            before = subprocess.check_output(before_command, shell = True, universal_newlines= True) # running the counting command and saving it to a string variable
            before = before.split(" ")[0]
            os.chdir(cwd) # changing directory back to Pipeline_Project
            os.system("gunzip {}_mapped_1.fq.gz".format(id))
            os.system("gunzip {}_mapped_2.fq.gz".format(id))
            after_command = "wc -l {}_mapped_1.fq".format(id) # counting the lines of the fastq files after bowtie
            after = subprocess.check_output(after_command, shell = True, universal_newlines= True) # running the counting command and saving output to string variable
            after = after.split(" ")[0] # split on empty space

            with open("PipelineProject.log","a") as f: # based on the SRR id of the file, write to the log file
                if "SRR5660030" in transcriptomes[i]:
                    f.write("Donor 1 (2dpi) had {0} read pairs before Bowtie2 filtering and {1} read pairs after.\n".format(str(int(before) // 4), str(int(after) // 4)))
                elif "SRR5660033" in transcriptomes[i]:
                    f.write("Donor 1 (6dpi) had {0} read pairs before Bowtie2 filtering and {1} read pairs after.\n".format(str(int(before) // 4), str(int(after) // 4)))
                elif "SRR5660044" in transcriptomes[i]:
                    f.write("Donor 3 (2dpi) had {0} read pairs before Bowtie2 filtering and {1} read pairs after.\n".format(str(int(before) // 4), str(int(after) // 4)))
                elif "SRR5660045" in transcriptomes[i]:
                    f.write("Donor 3 (6dpi) had {0} read pairs before Bowtie2 filtering and {1} read pairs after.\n".format(str(int(before) // 4), str(int(after) // 4)))
            os.chdir("{0}/{1}".format(cwd,inpath)) # change directory back to directory with transcriptomes fastq files
    assemble_spades(id_list,cwd) # go to step 3

## Step 3: Assemble Transcriptomes ##
def assemble_spades(id_list,cwd):
    os.chdir(cwd) # change to Pipeline_project
    # Run SPAdes with the Bowtie2 output files
    # This step assembles the inputted transcriptomes to produce one assembly
    spades_command = 'spades.py -k 127 -t 2 --only-assembler --pe-1 1 {0}_mapped_1.fq --pe-2 1 {0}_mapped_2.fq --pe-1 2 {1}_mapped_1.fq --pe-2 2 {1}_mapped_2.fq --pe-1 3 {2}_mapped_1.fq --pe-2 3 {2}_mapped_2.fq --pe-1 4 {3}_mapped_1.fq --pe-2 4 {3}_mapped_2.fq -o all4_assembly/'.format(id_list[0],id_list[1],id_list[2],id_list[3])
    os.system(spades_command)

    #write command to log file

    with open("PipelineProject.log","a") as f:
        f.write("{}\n".format(spades_command))
    count_contigs(cwd)

# Step 4: calculate length
# Calculate the # of contigs with a length > 1000
    
def count_contigs(cwd):
    contigs = list(SeqIO.parse(open("all4_assembly/contigs.fasta"), "fasta")) # parse SPAdes output contigs file
    count = 0 # initialize count of sequences over 1000 bp in length
    total_count = 0 # initialize total count of all lengths
    for contig in contigs: # iterate over contigs in file
        sequence = contig.seq # store sequence
        if len(sequence) > 1000: # check length of sequence
            count += 1 # add 1 if over 1000 bp 
            total_count += len(sequence) # add length to total length count

    # Write the numbers to the log file 
    with open("PipelineProject.log", "a") as f:
        f.write("There are {} contigs > 1000 bp in the assembly.\n".format(count))
        f.write("There are {} bp in the assembly.\n".format(total_count))
    step5(cwd) # Go to step 5

# Step 5: Running Blast on the longest contig
def step5(cwd):
    # Retrieve the longest contig from the SPAdes assembly
    contigs = list(SeqIO.parse(open("all4_assembly/contigs.fasta"), "fasta")) # open contigs files 
    seq_dict = {} # initialize sequence dict

    for contig in contigs: # iterate over contigs in file 
        seq_dict[len(contig.seq)] = contig.seq # make dictionary of lengths and sequences

    # write sequence to fasta file
    identifier_line = ">" + str(max(seq_dict)) + "\n" # create identifier line for fasta file
    sequence_line = str(seq_dict[max(seq_dict)]) # use sequence dict to find the max length and use that sequence in the fasta file
    
    with open("max_contig.fasta","w") as f: # write to fasta file
        f.write(identifier_line) # first identifier line 
        f.write(sequence_line + "\n") # sequence line

    # Make a local database of just sequences from the Betaherpesvirinae subfamily
    os.system("mkdir Betaherpesvirinae") # make new directory for Betaherpesvirinae
    os.chdir("Betaherpesvirinae") # change to that directory

    # command to download the sequences necessary for building the local database
    download_command = "datasets download virus genome taxon Betaherpesvirinae --include genome"
    os.system(download_command)
    os.system("unzip ncbi_dataset") # unzip data

    # command to make the local database with Betaherpesvirinae files
    makedb_command = "makeblastdb -in ncbi_dataset/data/genomic.fna -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl"
    os.system(makedb_command)
    # Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily  
    os.chdir(cwd) 

    # Blast command includes -max_hsps 1 flag to only include the best alignment (HSP) for any single query-subject pair of sequences
    # Includes -max_target_seqs 10 flag to only include top 10 hits in the output    
    blast_command = "blastn -query max_contig.fasta -db Betaherpesvirinae/Betaherpesvirinae -out longest_contig_results.csv -max_hsps 1 -max_target_seqs 10 -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    os.system(blast_command)
    with open("PipelineProject.log","a") as w: # write to .log file
        open_csv_cmd = "less longest_contig_results.csv" # uses output of less command to write to log file
        csv = subprocess.check_output(open_csv_cmd, shell = True, universal_newlines= True) # running the counting command and saving it to a string variable
        headers= "sacc   pident  length  qstart  qend    sstart  send    bitscore    evalue  stitle\n" # headers for the file, tab-separated
        w.write(headers) # headers line
        w.write(csv) # blast hits
    
# Start the wrapper script, running all of the steps
build_index(inpath) 