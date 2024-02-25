## Step 5: Does this assembly align with other virus strains? ##
print("Step 5: Does this assembly align with other virus strains?")
from Bio import SeqIO
import os

# Retrieve the longest contig from the SPAdes assembly
contigs = list(SeqIO.parse(open("all4_assembly/contigs.fasta"), "fasta"))
seq_dict = {}

for contig in contigs:
    seq_dict[len(contig.seq)] = contig.seq
print(max(seq_dict))

# write sequence to fasta file
identifier_line = ">" + str(max(seq_dict)) + "\n"
sequence_line = str(seq_dict[max(seq_dict)])
with open("max_contig.fasta","w") as f:
    f.write(identifier_line)
    f.write(sequence_line + "\n")
# Make a local database of just sequences from the Betaherpesvirinae subfamily
#download_command = "datasets download virus genome taxon Betaherpesvirinae --refseq --include genome"
#os.system(download_command)
#makedb_command = "makeblastdb -in ncbi_dataset/data/genomic.fna -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl"
#os.system(makedb_command)
# Use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily 
# Blast+ run should only keep the best alignment (HSP) for any single query-subject pair of sequences 

# For the top 10 hits, write the following to the log file:
    # Subject accession
    # Percent identity 
    # Alignment length 
    # Start of alignment in query
    # End of alignment in query
    # Start of alignment in subject
    # End of alignment in subject
    # Bit score
    # E-value
    # Subject Title
    
blast_command = "blastn -query max_contig.fasta -db Betaherpesvirinae -out longest_contig_results.csv -outfmt '6 sacc pident length qstart qend sstart send bitscore evalue stitle' -max_hsps 1"
os.system(blast_command)

