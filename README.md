# Pipeline Project

This pipeline is structured to compare HCMV transcriptomes 2 and 6 days post infection (dpi)

### You will need the following dependencies:
List all of the dependencies

## Step One: Retrieve transcriptomes from patients infected with HCMV 2 and 6 days post infection from SRA

Run the following commands to download the transcriptome SRA files from each Donor to your machine/server:

Donor 1 (2dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030

Donor 1 (6dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033

Donor 3 (2dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044

Donor 3 (6dpi): wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045

Run the following commands to uncompress the retrieved data:

fastq-dump -I --split_files SRR5660030

fastq-dump -I --split_files SRR5660033

fastq-dump -I --split_files SRR566004

fastq-dump -I --split_files SRR5660045