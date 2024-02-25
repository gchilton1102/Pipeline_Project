## Step 03 ##
import os
print("### Step 03: Assembling the Transcriptome")

# Using the output from Bowtie2, assemble all four transcriptomes together to produce one assembly
spades_command = 'nohup spades.py -k 77,99,127 -t 2 --only-assembler --pe-1 1 SRR5660030_mapped_1.fq --pe-2 1 SRR5660030_mapped_2.fq --pe-1 2 SRR5660033_mapped_1.fq --pe-2 2 SRR5660033_mapped_2.fq --pe-1 4 SRR5660045_mapped_1.fq --pe-2 4 SRR5660045_mapped_2.fq -o all4_assembly/ &'
os.system(spades_command)

#write command to log file
with open("PipelineProject.log","a") as f:
    f.write("{}\n".format(spades_command))