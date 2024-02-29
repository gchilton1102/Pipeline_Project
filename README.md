# Pipeline Project
This pipeline is structured to compare HCMV transcriptomes 2 and 6 days post infection (dpi)

## Getting Started
### You will need the following dependencies:
argparse

sys

os

subprocess

glob

- Separate installations, run the following on the command line:
pip install biopython

- Retrieve transcriptomes from patients infected with HCMV 2 and 6 days post infection from SRA

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

## Running the wrapper script

- Run the following on the command line:

python3 wrapper.py -i ~/path_to_transcriptome_fastqs

- Substitue "path_to_transcriptome_fastqs" with the full filepath to the transcriptome fastq files downloaded in the previous step

## Testing Pipeline with Sample Data

Sample data for this pipeline is included in the sample_data directory included in this GitHub repository

When testing with the sample data, include the full path to the sample data after the -i flag

The following is an example command to run with the sample data:

python3 wrapper.py ~/sample_data/
