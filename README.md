# Pipeline Project
This is a python-based pipeline that is structured to compare HCMV transcriptomes 2 and 6 days post infection (dpi) and assemble the genome based on the given transcriptomes.

## Getting Started

- To clone this directory, run the following on the command line:

    git clone https://github.com/gchilton1102/Pipeline_Project_Grace_Chilton.git

- If you do not already have Python installed, follow the instructions on this page: https://www.python.org/downloads/

### You will need the following dependencies:

- The following libraries are part of the standard library of python version 2.7 or later:

    argparse

    sys

    os

    subprocess

    glob

- Separate installations, run the following on the command line:

    pip install biopython

### Retrieve transcriptomes from patients infected with HCMV 2 and 6 days post infection from SRA

- Run the following commands to download the transcriptome SRA files from each Donor to your machine/server:

    Donor 1 (2dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

    Donor 1 (6dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

    Donor 3 (2dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

    Donor 3 (6dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

- Run the following commands to uncompress the retrieved data:

    fastq-dump -I --split_files SRR5660030

    fastq-dump -I --split_files SRR5660033

    fastq-dump -I --split_files SRR566004

    fastq-dump -I --split_files SRR5660045

- Move all fastq files to the Pipeline_Project_Grace_Chilton directory, created with the git clone command

## Running the wrapper script

- Run the following on the command line:

    python3 wrapper.py -i name_of_transcriptome_folder/

- Substitue "name_of_transcriptome_folder" with the full filepath to the transcriptome fastq files downloaded in the previous step

- Make sure your fastq files are within the Pipeline_Project_Grace_Chilton directory

## Testing Pipeline with Sample Data

Sample data for this pipeline is included in the sample_data directory included in this GitHub repository

When testing with the sample data, include the sample data after the -i flag

The following is the command to run with the sample data on the command line, make sure you are working in the Pipeline_Project_Grace_Chilton directory:

python3 wrapper.py -i sample_data/

