## Step 04: Calculate the number of contigs with a length > 1000 ##
from Bio import SeqIO
from Bio.Seq import Seq

print("### Step 04: Calculate the number of contigs > 1000 & Write the total number of base pairs of these contigs ###")
# Calculate the # of contigs with a length > 1000
contigs = list(SeqIO.parse(open("all4_assembly/contigs.fasta"), "fasta"))
count = 0
total_count = 0
for contig in contigs:
    sequence = contig.seq
    if len(sequence) > 1000:
        count += 1
        total_count += len(sequence)

# Write the # to the log file 
with open("PipelineProject.log", "a") as f:
    f.write("There are {} contigs > 1000 bp in the assembly.\n".format(count))
    f.write("There are {} bp in the assembly.\n".format(total_count))
# Calculate the length of the assembly(the total # of bp in all of the contigs > 1000 in length)

# Write the length to the log file