## Step 2 ##
import os
import subprocess

print("### Step 2: Building the HCMV Index and Mapping Transcriptomes ###")
# Create an index for HCMV
#index_command = "datasets download genome accession GCA_027926585.1 --include gff3,rna,cds,protein,genome,seq-report"
#os.system(index_command) 
# Save only the reads that map to the HCMV index for use in assembly
os.system("rm PipelineProject.log")
transcriptomes = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]
k = 1
for transcriptome in transcriptomes:
    #map_command = "bowtie2 --quiet -x HCMV -1 ~/transcriptomes/{0}_1.fastq -2 ~/transcriptomes/{0}_2.fastq -S {0}.sam --al-conc-gz {0}_mapped_%.fq.gz".format(transcriptome)
    #os.system(map_command)
    # Calculate the number of read pairs before and after bowtie filtering
    before_command = "wc -l transcriptomes/{}_1.fastq".format(transcriptome)
    before = subprocess.check_output(before_command, shell = True)
    after_command = "wc -l {}.sam".format(transcriptome)
    after = subprocess.check_output(after_command, shell = True)
    before_num = str(before).split("transcriptomes")[0]
    before_num = before_num.split("'")[1]
    after_num = str(after).split("transcriptomes")[0]
    after_num = after_num.split("'")[1]
    after_num = after_num.split(" ")[0]
    with open("PipelineProject.log","a") as f:
        f.write("Donor {} had {}read pairs before Bowtie2 filtering and {} read pairs after\n".format(k, before_num, after_num))
    k += 1