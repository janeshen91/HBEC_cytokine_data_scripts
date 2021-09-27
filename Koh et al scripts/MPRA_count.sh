#!/bin/bash
source /home/jane/miniconda3/etc/profile.d/conda.sh
cd ~//MPRAflow
conda activate MPRAflow

nextflow -bg run count.nf -w "/data/jane/MPRAflow/Count_Basic/work/kd2NewAssoc" \
--experiment-file "/data/jane/MPRAflow/Feb2021/KD2expression/KD2designCSV.csv" \
--dir "/data/jane/MPRAflow/Feb2021/KD2expression/" \
--outdir "/data/jane/MPRAflow/Count_Basic/KD2NewAssoc" \
--design "/data/jane/MPRAflow/Feb2021/KD2expression/KD2_MPRA_RefFile.fa" \
--umi-length 16 \
--association "/data/jane/MPRAflow/Assoc_Basic/workKD2MAPQ0/c6/1eb0adfd7d45c7f96040eb9f6b6b71/MAPQANDNM/KD2baseq35200MandNM0_filtered_coords_to_barcodes.pickle"