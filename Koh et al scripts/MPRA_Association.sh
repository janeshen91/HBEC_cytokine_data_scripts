#!/bin/bash
source /home/jane/miniconda3/etc/profile.d/conda.sh

#cd /home/jane/MPRAflow/
#conda activate MPRAflow
#nextflow -bg run association.nf \
#-w /data/jane//MPRAflow/Assoc_Basic/workKD2 \
#--fastq-insert "/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R1_001.fastq.gz"  \
#--fastq-insertPE "/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R3_001.fastq.gz" \
#--fastq-bc "/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R2_001.fastq.gz" \
#--design "//data/jane/MPRAflow/Feb2021/KD2expression/KD2_MPRA_RefFile.fa" \
#--name KD2.Feb2021.AssocBasicCIGAR  \
#--outdir /data/jane//MPRAflow/outs/KD2Cigar \
#--cigar 200M


#cd /home/jane/MPRAflow/
#conda activate MPRAflow
#nextflow -bg run association.nf \
#-w /data/jane//MPRAflow/Assoc_Basic/workKD2 \
#--fastq-insert "/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R1_001.fastq.gz"  \
#--fastq-insertPE "/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R3_001.fastq.gz" \
#--fastq-bc "/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R2_001.fastq.gz" \
#--design "//data/jane/MPRAflow/Feb2021/ID614.fa" \
#--name KD2.Feb2021.AssocBasicID614 \
#--outdir /data/jane//MPRAflow/outs/KD2ID614 \
#--cigar 200M


cd /home/jane/MPRAflow/
conda activate MPRAflow
nextflow -bg run association.nf \
-w /data/jane//MPRAflow/Assoc_Basic/workOYB3 \
--fastq-insert "/data/jane/MPRAflow/Sept2021/fastq/2021_0017_S1_R1_001.fastq.gz"  \
--fastq-insertPE "/data/jane/MPRAflow/Sept2021/fastq/2021_0017_S1_R3_001.fastq.gz" \
--fastq-bc "/data/jane/MPRAflow/Sept2021/fastq/2021_0017_S1_R2_001.fastq.gz" \
--design "/data/jane/MPRAflow/Sept2021/OYB3_LEAN_DESIGN_REF.fa" \
--name OYB3.AssocBasic \
--outdir /data/jane//MPRAflow/outs/KD2SatMut 


baseDir="/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71"
fastq_bc="/data/jane/MPRAflow/Feb2021/FASTQ2021_003/2021_003_S1_R2_001.fastq.gz"
bam="/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/toRunAndTest/s_merged.bam"
count_fastq="/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/count_fastq.txt"
count_bam="/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/toRunAndTest/count_merged.txt"
name="KD2baseq30NM0200M"
mapq=0
baseq=30
cigar='200M'
#samtools view -S -b sample.sam > sample.bam

cd /data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/MAPQANDNM


python /data/jane/MPRAflow/nf_ori_map_barcodes3.py ${baseDir} ${fastq_bc} ${count_fastq} $bam ${count_bam} ${name} ${mapq} ${baseq} ${cigar}

out="KD2baseq35200MandNM0"
map="/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/MAPQANDNM/KD2baseq30NM0200M_coords_to_barcodes.pickle"
table="/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/MAPQANDNM/KD2baseq30NM0200M_barcodes_per_candidate-no_repeats-no_jackpots.feather"
min_cov=3
min_frac=1
label="/data/jane/MPRAflow/Feb2021/KD2labels.txt"

python /data/jane/MPRAflow/nf_ori_map_barcodes3.py ${baseDir} ${fastq_bc} ${count_fastq} $bam ${count_bam} ${name} ${mapq} ${baseq} ${cigar}
