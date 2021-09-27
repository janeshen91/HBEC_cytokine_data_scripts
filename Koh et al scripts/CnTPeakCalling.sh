##this is the script for calling peaks on KLF CnT data

#first clone the alignall script,which uses STAR to align the fastq reads to the reference genome, and then removes duplicates
git clone -b test_mmu_snp https://janeshen91@bitbucket.org/joshucsf/alignall.git

##peak calling 
dir="/data/jane/KDcutNtag/Mar2021run/KLF5samples/"
allKLF5files=("KLF5_1043_IL13_merged.bam" "KLF5_1043_UT_merged.bam" "KLF5_1075_IL13_merged.bam" "KLF5_1075_UT_merged.bam" "KLF5_1332_IL13_merged.bam" "KLF5_1332_UT_merged.bam")

cd $dir
for file in `ls $dir`; do
echo $file
file2="${file}_p0p01"
macs2 callpeak -t $file -f BAM --keep-dup 1 -n $file  -p 0.01 --nomodel --outdir ./$file2
done


##then we run IDR to find the reproducible set
idr --samples /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak  \
/data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_firstRun_lowP0p01/KLF5_first_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S23_acrossSeqRuns \
--plot \
--log-output-file KLF5_UT_S23_acrossSeqRuns.log


cd /data/jane/KDcutNtag/Mar2021run/KLF5samples
#Sort and then run IDR
sort -k8,8nr KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.narrowPeak \
> KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.narrowPeak \
> KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.narrowPeak \
> KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.narrowPeak \
> KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.sorted.narrowPeak



idr --samples KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.sorted.narrowPeak  \
KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_1043_1075_acrossSeqRuns \
--plot \
--peak-list KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_1043_1075_acrossSeqRuns.log

idr --samples KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.sorted.narrowPeak  \
KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_1043_1332_acrossSeqRuns \
--plot \
--peak-list KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_1043_1332_acrossSeqRuns.log

idr --samples KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.sorted.narrowPeak  \
KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_1075_1332_acrossSeqRuns \
--plot \
--peak-list KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_1075_1332_acrossSeqRuns.log


#change into bed files and then intersect
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_1043_1075_acrossSeqRunsSig.txt >KLF5_UT_1043_1075_acrossSeqRunsSig.bed
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_1043_1332_acrossSeqRunsSig.txt >KLF5_UT_1043_1332_acrossSeqRunsSig.bed
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_1075_1332_acrossSeqRunsSig.txt >KLF5_UT_1075_1332_acrossSeqRunsSig.bed

bedtools intersect -a KLF5_UT_1043_1075_acrossSeqRunsSig.bed \
-b KLF5_UT_1043_1332_acrossSeqRunsSig.bed -wo \
>KLF5_UT_1043_1075_1043_1332_acrossSeqRunsIntersect.bed

bedtools intersect -a KLF5_UT_1043_1075_1043_1332_acrossSeqRunsIntersect.bed \
-b KLF5_UT_1075_1332_acrossSeqRunsSig.bed -wo \
>KLF5_UT_all3_acrossSeqRunsIntersect.bed

#NOW WE want to do all of this for IL13

KLF5_1043_IL13_merged.bam_p0p01

sort -k8,8nr KLF5_1043_IL13_merged.bam_q0p1/KLF5_1043_IL13_merged.bam_peaks.narrowPeak \
> KLF5_1043_IL13_merged.bam_q0p1/KLF5_1043_IL13_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.narrowPeak \
> KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_1332_IL13_merged.bam_q0p1/KLF5_1332_IL13_merged.bam_peaks.narrowPeak \
> KLF5_1332_IL13_merged.bam_q0p1/KLF5_1332_IL13_merged.bam_peaks.sorted.narrowPeak


macs2 callpeak -t KLF5_IL13_merged.bam -f BAM \
--keep-dup 1 -n KLF5_IL13_mergedP0p01  -p 0.01 --nomodel --outdir KLF5_IL13_mergedP0p01

sort -k8,8nr KLF5_IL13_mergedP0p01/KLF5_IL13_mergedP0p01_peaks.narrowPeak \
> KLF5_IL13_mergedP0p01/KLF5_IL13_mergedP0p01_peaks.sorted.narrowPeak



idr --samples KLF5_1043_IL13_merged.bam_q0p1/KLF5_1043_IL13_merged.bam_peaks.sorted.narrowPeak  \
KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_IL13_1043_1075_acrossSeqRuns \
--plot \
--peak-list KLF5_IL13_mergedP0p01/KLF5_IL13_mergedP0p01_peaks.sorted.narrowPeak \
--log-output-file KLF5_IL13_1043_1075_acrossSeqRuns.log

idr --samples KLF5_1043_IL13_merged.bam_q0p1/KLF5_1043_IL13_merged.bam_peaks.sorted.narrowPeak  \
KLF5_1332_IL13_merged.bam_q0p1/KLF5_1332_IL13_merged.bam_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_IL13_1043_1332_acrossSeqRuns \
--plot \
--peak-list KLF5_IL13_mergedP0p01/KLF5_IL13_mergedP0p01_peaks.sorted.narrowPeak \
--log-output-file KLF5_IL13_1043_1332_acrossSeqRuns.log

idr --samples KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.sorted.narrowPeak  \
KLF5_1332_IL13_merged.bam_q0p1/KLF5_1332_IL13_merged.bam_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_IL13_1075_1332_acrossSeqRuns \
--plot \
--peak-list KLF5_IL13_mergedP0p01/KLF5_IL13_mergedP0p01_peaks.sorted.narrowPeak \
--log-output-file KLF5_IL13_1075_1332_acrossSeqRuns.log


awk '{if($5 >= 540) print $0}' KLF5_IL13_1043_1075_acrossSeqRuns >KLF5_IL13_1043_1075_acrossSeqRunsSig.txt
awk '{if($5 >= 540) print $0}' KLF5_IL13_1043_1332_acrossSeqRuns >KLF5_IL13_1043_1332_acrossSeqRunsSig.txt
awk '{if($5 >= 540) print $0}' KLF5_IL13_1075_1332_acrossSeqRuns >KLF5_IL13_1075_1332_acrossSeqRunsSig.txt

#change into bed files and then intersect
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <KLF5_IL13_1043_1075_acrossSeqRunsSig.txt >KLF5_IL13_1043_1075_acrossSeqRunsSig.bed
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <KLF5_IL13_1043_1332_acrossSeqRunsSig.txt >KLF5_IL13_1043_1332_acrossSeqRunsSig.bed
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <KLF5_IL13_1075_1332_acrossSeqRunsSig.txt >KLF5_IL13_1075_1332_acrossSeqRunsSig.bed

bedtools intersect -a KLF5_IL13_1043_1075_acrossSeqRunsSig.bed \
-b KLF5_IL13_1043_1332_acrossSeqRunsSig.bed -wo \
>KLF5_IL13_1043_1075_1043_1332_acrossSeqRunsIntersect.bed

bedtools intersect -a KLF5_IL13_1043_1075_1043_1332_acrossSeqRunsIntersect.bed \
-b KLF5_IL13_1075_1332_acrossSeqRunsSig.bed -wo \
>KLF5_IL13_all3_acrossSeqRunsIntersect.bed

bedtools intersect -a

##combining the 2 all 3 intersect beds ie concatenate thhen merge them
cat KLF5_IL13_all3_acrossSeqRunsIntersect.bed KLF5_UT_all3_acrossSeqRunsIntersect.bed > KLF5_both_all3_acrossSeqRunsIntersect.bed
bedtools sort -i KLF5_both_all3_acrossSeqRunsIntersect.bed > KLF5_both_all3_acrossSeqRunsIntersectSorted.bed
bedtools intersect -a KLF5_IL13_all3_acrossSeqRunsIntersect.bed -b KLF5_UT_all3_acrossSeqRunsIntersect.bed >KLF5_common_all3_acrossSeqRunsIntersect.bed

bedtools merge -i KLF5_both_all3_acrossSeqRunsIntersectSorted.bed >KLF5_both_all3_acrossSeqRunsIntersectMerged.bed #THIS IS THE ALL file
##where is the common file?

cd /data/jane/KDcutNtag/Mar2021run/KLF5samples/
samtools merge -u KLF5_UT_IL13_merged.bam KLF5_UT_merged.bam KLF5_IL13_merged.bam
coverageBed -abam KLF5_UT_IL13_merged.bam -b KLF5_both_all3_acrossSeqRunsIntersectMerged.bed -counts > KLF5_both_all3_acrossSeqRunsIntersectMerged.cnt.bed
#take chr out of bed file
perl -ane 'print ""'
bedtools coverage -a KLF5_both_all3_acrossSeqRunsIntersectMerged.bed -b KLF5_UT_IL13_merged.bam >KLF5_both_all3_acrossSeqRunsIntersectMergedii.txt
bedtools coverage -a KLF5_both_all3_acrossSeqRunsIntersectMerged.bed -b KLF5_UT_merged.bam >countsOfKLF5_UT.txt
bedtools coverage -a KLF5_both_all3_acrossSeqRunsIntersectMerged.bed -b KLF5_IL13_merged.bam >countsOfKLF5_IL13.txt


##TRYING TO DO DIFF PEAK ON UNION OF PEAKS
#do I just concatenate them and then merge them?
cat KLF5_UT_1043_1075_acrossSeqRunsSig.bed KLF5_UT_1043_1332_acrossSeqRunsSig.bed KLF5_UT_1075_1332_acrossSeqRunsSig.bed >KLF5_UT_all3_acrossSeqRunsCat.bed
bedtools sort -i KLF5_UT_all3_acrossSeqRunsCat.bed >KLF5_UT_all3_acrossSeqRunsCat.sorted.bed
bedtools merge -i KLF5_UT_all3_acrossSeqRunsCat.sorted.bed>KLF5_UT_all3_acrossSeqRunsCat.merged.bed


cat KLF5_IL13_1043_1075_acrossSeqRunsSig.bed KLF5_IL13_1043_1332_acrossSeqRunsSig.bed KLF5_IL13_1075_1332_acrossSeqRunsSig.bed >KLF5_IL13_all3_acrossSeqRunsCat.bed
bedtools sort -i KLF5_IL13_all3_acrossSeqRunsCat.bed >KLF5_IL13_all3_acrossSeqRunsCat.sorted.bed
bedtools merge -i KLF5_IL13_all3_acrossSeqRunsCat.sorted.bed>KLF5_IL13_all3_acrossSeqRunsCat.merged.bed




bedtools coverage -a KLF5_IL13_all3_acrossSeqRunsCat.noChr.bed -b KLF5_IL13_merged.bam >countsOfKLF5_IL13.txt
bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_1332_IL13_merged.bam >countsOfKLF5_1332_IL13.txt
#bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_IL13_merged.bam >countsOfKLF5_IL13_all.txt
bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_1043_IL13_merged.bam >countsOfKLF5_1043_IL13_all.txt
bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_1075_IL13_merged.bam >countsOfKLF5_1075_IL13_all.txt

bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_1332_UT_merged.bam >countsOfKLF5_1332_UT.txt
#bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_UT_merged.bam >countsOfKLF5_UT_all.txt
bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_1043_UT_merged.bam >countsOfKLF5_1043_UT_all.txt
bedtools coverage -a KLF5_acrossRuns.merged.bed -b KLF5_1075_UT_merged.bam >countsOfKLF5_1075_UT_all.txt



##all
##common
bedtools intersect -a KLF5_IL13_all3_acrossSeqRunsCat.merged.bed -b KLF5_UT_all3_acrossSeqRunsCat.merged.bed -wo >KLF5_intersect_all3_acrossSeqRunsCat.bed
cat KLF5_IL13_all3_acrossSeqRunsCat.noChr.bed KLF5_UT_all3_acrossSeqRunsCat.noChr.bed >KLF5_acrossRuns.noChr.bed
bedtools sort -i KLF5_acrossRuns.noChr.bed >KLF5_acrossRuns.sorted.bed
bedtools merge -i KLF5_acrossRuns.sorted.bed >KLF5_acrossRuns.merged.bed

bedtools intersect -a KLF5_IL13_all3_acrossSeqRunsCat.noChr.bed -b KLF5_UT_all3_acrossSeqRunsCat.noChr.bed -wo >KLF5_acrossRuns.Intersect.noChr.bed
##ok finally have all the files I think? but I think we also need in dividual coverage files




