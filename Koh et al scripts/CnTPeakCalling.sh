##merge bam files for KLF5 UT
#S9, S17, S23
cd /data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments
samtools merge -u KLF5_UT_merged.bam ./STAR210301_153239/KLF5_S9_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam \
./STAR210303_144341/KLF5_S17_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam \
./STAR210305_110023/KLF5_S23_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam 

##calling peaks with a lower threshold
cd /data/jane/KDcutNtag/Mar2021run/
macs2 callpeak -t ./AlignAllOut/Alignments/STAR210301_153239/KLF5_S9_UTout//Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5_S9_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_S9_lowP0p01



cd /data/jane/KDcutNtag/Mar2021run/
macs2 callpeak -t ./AlignAllOut/Alignments//STAR210303_144341/KLF5_S17_UTout///Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5_S17_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_S17_lowP0p01



cd /data/jane/KDcutNtag/Mar2021run/
macs2 callpeak -t ./AlignAllOut/Alignments//STAR210305_110023/KLF5_S23_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5_S23_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_S23_lowP0p01

macs2 callpeak -t ./AlignAllOut/Alignments/KLF5_UT_merged.bam -f BAM \
--keep-dup 1 -n KLF5_merge  -q 0.1 --nomodel --outdir ./Batch210301_153239/KLF5_S23_mergeQop1

macs2 callpeak -t ./AlignAllOut/Alignments/KLF5_UT_merged.bam -f BAM \
--keep-dup 1 -n KLF5_merge  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_S23_mergelowP

#so compareing 1 to 2, 2 to 3 and 1 to 3
#and then??
#compare it to KLF5f?

sort -k8,8nr /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergeQop1/KLF5_merge_peaks.narrowPeak \
>/data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergeQop1/KLF5_merge_peaks.sorted.narrowPeak

sort -k8,8nr /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergelowP/KLF5_merge_peaks.narrowPeak \
>/data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergelowP/KLF5_merge_peaks.sorted.narrowPeak

cd /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/
sort -k8,8nr ./KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.narrowPeak \
>./KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak

cd /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/
sort -k8,8nr ./KLF5_S17_lowP0p01/KLF5_S17_lowP0p01_peaks.narrowPeak \
>./KLF5_S17_lowP0p01/KLF5_S17_lowP0p01_peaks.sorted.narrowPeak

cd /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/
sort -k8,8nr ./KLF5_S9_lowP0p01/KLF5_S9_lowP0p01_peaks.narrowPeak \
>./KLF5_S9_lowP0p01/KLF5_S9_lowP0p01_peaks.sorted.narrowPeak


sort -k8,8nr /data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout_peaks.narrowPeak \
> /data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout_peaks.sorted.narrowPeak 

idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S9_lowP0p01//KLF5_S9_lowP0p01_peaks.sorted.narrowPeak \
//data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S17_lowP0p01/KLF5_S17_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S9_S17_looseP_refMerge_idr \
--plot \
--peak-list /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergeQop1/KLF5_merge_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_S9_S17_looseP_refMerge_idr.log

awk '{if($5 >= 540) print $0}' KLF5_UT_S9_S17_looseP_refMerge_idr >KLF5_UT_S9_S17_looseP_refMerge_idrSig.txt | wc -l


##cd Batch210303_144344
idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S9_lowP0p01//KLF5_S9_lowP0p01_peaks.sorted.narrowPeak \
//data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S9_S23_looseP_refMerge_idr \
--plot \
--peak-list /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergeQop1/KLF5_merge_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_S9_S23_looseP_refMerge_idr.log
awk '{if($5 >= 540) print $0}' KLF5_UT_S9_S23_looseP_refMerge_idr >KLF5_UT_S9_S23_looseP_refMerge_idrSig.txt 
wc KLF5_UT_S9_S23_looseP_refMerge_idrSig.txt 
#11768


idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S17_lowP0p01//KLF5_S17_lowP0p01_peaks.sorted.narrowPeak \
//data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S17_S23_looseP_refMerge_idr \
--plot \
--peak-list /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergeQop1/KLF5_merge_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_S17_S23_looseP_refMerge_idr.log

awk '{if($5 >= 540) print $0}' KLF5_UT_S17_S23_looseP_refMerge_idr >KLF5_UT_S17_S23_looseP_refMerge_idrSig.txt 


##are the ones that pass better?
#check motif analysis
#loose P merge

cd IDR
idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S9_lowP0p01//KLF5_S9_lowP0p01_peaks.sorted.narrowPeak \
//data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S17_lowP0p01/KLF5_S17_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S9_S17_lowPForMergeToo \
--plot \
--peak-list /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergelowP/KLF5_merge_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_S9_S17_lowPForMergeToo.log

awk '{if($5 >= 540) print $0}' KLF5_UT_S9_S17_lowPForMergeToo >KLF5_UT_S9_S17_lowPForMergeTooSig.txt

idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S9_lowP0p01//KLF5_S9_lowP0p01_peaks.sorted.narrowPeak \
//data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S9_S23_lowPForMergeToo \
--plot \
--peak-list /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergelowP/KLF5_merge_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_S9_S23_lowPForMergeToo.log

awk '{if($5 >= 540) print $0}' KLF5_UT_S9_S23_lowPForMergeToo >KLF5_UT_S9_S23_lowPForMergeTooSig.txt

idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//KLF5_S17_lowP0p01/KLF5_S17_lowP0p01_peaks.sorted.narrowPeak \
/data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S17_S23_lowPForMergeToo \
--plot \
--peak-list /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_mergelowP/KLF5_merge_peaks.sorted.narrowPeak \
--log-output-file KLF5_UT_S17_S23_lowPForMergeToo.log

awk '{if($5 >= 540) print $0}' KLF5_UT_S17_S23_lowPForMergeToo >KLF5_UT_S17_S23_lowPForMergeTooSig.txt

##how do we know what is the final region of voerlaps?\
#I guess we can intersect all 3 to get the final region


perl ./configureHomer.pl -install homer
PATH=$PATH:/data/jane/HOMER/.//bin/

perl -ane '$count++; print "$F[0]\t$F[1]\t$F[2]\tPeak$count\n";' \
< //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout_peaks.narrowPeak.bed \
> //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout.narrowPeak.ii.bed
nohup findMotifsGenome.pl //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout.narrowPeak.ii.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/HOMER/KLF5_S9 -size given >./nohupKLF5IL13HOMER.txt &

perl -ane '$count++; print "$F[0]\t$F[1]\t$F[2]\tPeak$count\n";' \
< //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout_peaks.narrowPeak.bed \
> //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout.narrowPeak.ii.bed
nohup findMotifsGenome.pl //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S9_UTout/KLF5_S9_UTout.narrowPeak.ii.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/HOMER/KLF5_S9 -size given >./nohupKLF5IL13HOMER.txt &

file="/data/jane/KDcutNtag/Mar2021run/BatchSTAR210305_110023/KLF5.S23.IL13/KLF5.S23.IL13_peaks.narrowPeak"
perl -ane 'print "chr$F[0]\t$F[1]\t$F[2]\n";' <$file >$file.bed

perl -ane '$count++; print "$F[0]\t$F[1]\t$F[2]\tPeak$count\n";' \
< $file.bed \
> $file.ii.bed

nohup findMotifsGenome.pl $file.ii.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/HOMER/KLF5_S23 -size given >./nohupKLF5IL13HOMERs23.txt &

##get the intersect of all 3 idrsig
cd /data/jane/KDcutNtag/Mar2021run/IDR
perl -ane  'print "$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_S17_S23_lowPForMergeTooSig.txt >KLF5_UT_S17_S23_lowPForMergeTooSig.bed
perl -ane  'print "$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_S9_S23_lowPForMergeTooSig.txt >KLF5_UT_S9_S23_lowPForMergeTooSig.bed
perl -ane  'print "$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_S9_S17_lowPForMergeTooSig.txt >KLF5_UT_S9_S17_lowPForMergeTooSig.bed

bedtools intersect -a KLF5_UT_S17_S23_lowPForMergeTooSig.bed \
-b KLF5_UT_S9_S23_lowPForMergeTooSig.bed \
-wo >KLF5_UT_S17S23_S9S23.intersect.bed

bedtools intersect -a KLF5_UT_S9_S17_lowPForMergeTooSig.bed \
-b KLF5_UT_S17S23_S9S23.intersect.bed \
-wo >KLF_all3.bed

perl -ane '$count++; print "chr$F[0]\t$F[1]\t$F[2]\t$count\n";' <KLF_all3.bed >KLF_all3.bed
nohup findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/IDR/KLF_all3.bed\
hg38 /data/jane/KDcutNtag/Mar2021run/HOMER/KLF5_all3 -size given >./nohupKLF5IL13HOMERall3UT.txt &

##
perl -ane  'print "$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_S17_S23_looseP_refMerge_idrSig.txt >KLF5_UT_S17_S23_looseP_refMerge_idrSig.bed
perl -ane  'print "$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_S9_S17_looseP_refMerge_idrSig.txt >KLF5_UT_S9_S17_looseP_refMerge_idrSig.bed
perl -ane  'print "$F[0]\t$F[1]\t$F[2]\n";' <KLF5_UT_S9_S23_looseP_refMerge_idrSig.txt >KLF5_UT_S9_S23_looseP_refMerge_idrSig.bed

bedtools intersect -a KLF5_UT_S17_S23_looseP_refMerge_idr.bed \
-b KLF5_UT_S9_S17_looseP_refMerge_idrSig.bed \
-wo >KLF5_UT_S17S23_S9S23.intersect.QRef.bed

bedtools intersect -a KLF5_UT_S17S23_S9S23.intersect.QRef.bed \
-b KLF5_UT_S9_S23_looseP_refMerge_idrSig.bed \
-wo >KLF5_UT_all3.intersect.QRef.bed


###let's see if the 2 KLF and KLF5f's are compatible, and see if we are 
#ie s13 and s14
##and also check like the concordance between this and the last time
##

perl -ane '$count++; print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n";' \
< //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S13_IL13out/KLF5_S13_IL13out_peaks.narrowPeak \
> //data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S13_IL13out/KLF5_S13_IL13out_peaks.narrowPeak.bed
##call the peaks with looser values

cd /data/jane/KDcutNtag/Mar2021run/
macs2 callpeak -t ./AlignAllOut/Alignments/STAR210301_153239/KLF5_S13_IL13out//Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5_S13_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_S13_lowP0p01

cd /data/jane/KDcutNtag/Mar2021run/
macs2 callpeak -t /data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments/STAR210303_144341/KLF5f_S14_IL13out/Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5f_S14_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5f_S14_lowP0p01

cd /data/jane/KDcutNtag/Mar2021run/IDR
idr --samples /data/jane/KDcutNtag/Mar2021run/./Batch210301_153239/KLF5f_S14_lowP0p01/KLF5f_S14_lowP0p01_peaks.narrowPeak \
//data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S13_lowP0p01/KLF5_S13_lowP0p01_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_IL13_S13_S14_lowPForMergeToo \
--plot \
--log-output-file KLF5_IL13_S13_S14_lowPForMergeToo.log

##merge things pairwise??
#KLF5 S23 with the ut klf5
cd /data/jane/KDcutNtag/Mar2021run/
macs2 callpeak -t /data/jane/KDcutNtag/AlignAllOut/Alignments/STAR201214_165848/KLF5_UT_rep1out//Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5_first_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_firstRun_lowP0p01

macs2 callpeak -t /data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments/STAR210305_110023/KLF5_S23_UTout//Aligned.sortedByCoord.duplicateRemoved.out.bam -f BAM \
--keep-dup 1 -n KLF5_S23_lowP0p01  -p 0.01 --nomodel --outdir ./Batch210301_153239/KLF5_S23_lowP0p01


sort -k8,8nr /data/jane/KDcutNtag/Mar2021run//Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.narrowPeak \
> /data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak 


sort -k8,8nr /data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_firstRun_lowP0p01//KLF5_firstRun_lowP0p01t_peaks.narrowPeak \
> /data/jane/KDcutNtag/Mar2021run/Batch210301_153239//Batch210301_153239/KLF5_firstRun_lowP0p01/KLF5_firstRun_lowP0p01_peaks.sorted.narrowPeak 

sort -k8,8nr  KLF5_first_lowP0p01_peaks.narrowPeak >KLF5_first_lowP0p01_peaks.sorted.narrowPeak

idr --samples /data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak  \
/data/jane/KDcutNtag/Mar2021run/Batch210301_153239/KLF5_firstRun_lowP0p01/KLF5_first_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S23_acrossSeqRuns \
--plot \
--log-output-file KLF5_UT_S23_acrossSeqRuns.log

##combing bam files across samples
cd /data/jane/KDcutNtag/Mar2021run/KLF5samples
cd /data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments/
##combing S9 and S10

samtools merge -u /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_1043_UT_merged.bam ./STAR210301_153239/KLF5_S9_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam \
./STAR210302_102841/KLF5f_S10_UT/Aligned.sortedByCoord.duplicateRemoved.out.bam 

samtools merge -u /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_1043_IL13_merged.bam ./STAR210301_153239/KLF5_S13_IL13out/Aligned.sortedByCoord.duplicateRemoved.out.bam \
./STAR210303_144341/KLF5f_S14_IL13out/Aligned.sortedByCoord.duplicateRemoved.out.bam 

samtools merge -u /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_1332_UT_merged.bam ./STAR210303_144341/KLF5_S17_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam \
./STAR210301_153239/KLF5f_S18_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam 

samtools merge -u /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_1332_IL13_merged.bam ./STAR210305_110023/KLF5_S21_IL13out/Aligned.sortedByCoord.duplicateRemoved.out.bam \
./STAR210301_153239/KLF5f_S22_IL13out/Aligned.sortedByCoord.duplicateRemoved.out.bam 

samtools merge -u /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_1075_UT_merged.bam ./STAR210305_110023/KLF5_S23_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam \
/data/jane/KDcutNtag/AlignAllOut/Alignments/STAR201214_165848/KLF5_UT_rep1out/Aligned.sortedByCoord.duplicateRemoved.out.bam 

samtools merge -u /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_1075_IL13_merged.bam ./STAR210303_144341/KLF5_S24_IL13out/Aligned.sortedByCoord.duplicateRemoved.out.bam \
/data/jane/KDcutNtag/AlignAllOut/Alignments/STAR201214_165848//KLF5_IL13_rep1out/Aligned.sortedByCoord.duplicateRemoved.out.bam 

cd /data/jane/KDcutNtag/Mar2021run/KLF5samples/
samtools merge -u KLF5_IL13_merged.bam KLF5_1075_IL13_merged.bam KLF5_1332_IL13_merged.bam KLF5_1043_IL13_merged.bam
samtools merge -u KLF5_UT_merged.bam KLF5_1075_UT_merged.bam KLF5_1332_UT_merged.bam KLF5_1043_UT_merged.bam

macs2 callpeak -t KLF5_IL13_merged.bam -f BAM --keep-dup 1 -n KLF5_IL13_merged  -p 0.01 --nomodel --outdir ./KLF5_IL13_merged0p01
macs2 callpeak -t KLF5_IL13_merged.bam -f BAM --keep-dup 1 -n KLF5_IL13_merged  -q 0.1 --nomodel --outdir ./KLF5_IL13_mergedQ0p1

#peak calling at loose p value
dir="/data/jane/KDcutNtag/Mar2021run/KLF5samples/"
allKLF5files=("KLF5_1043_IL13_merged.bam" "KLF5_1043_UT_merged.bam" "KLF5_1075_IL13_merged.bam" "KLF5_1075_UT_merged.bam" "KLF5_1332_IL13_merged.bam" "KLF5_1332_UT_merged.bam")

cd $dir
for file in `ls $dir`; do
echo $file
file2="${file}_p0p01"
file3="${file}_q0p1"
macs2 callpeak -t $file -f BAM --keep-dup 1 -n $file  -p 0.01 --nomodel --outdir ./$file2
macs2 callpeak -t $file -f BAM --keep-dup 1 -n $file  -q 0.1 --nomodel --outdir ./$file3
done

for file in ${allKLF5files[@]}; do
echo $file
file3="${file}_q0p1"
file4="${file}_peaks.narrowPeak"
scp jane@10.37.30.76:/data/jane/KDcutNtag/Mar2021run/KLF5samples/$file3/$file4 .
done
##how to get the bam.bai for the merged

##then we run IDR to find the reproducible set
idr --samples /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_S23_lowP0p01/KLF5_S23_lowP0p01_peaks.sorted.narrowPeak  \
/data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_firstRun_lowP0p01/KLF5_first_lowP0p01_peaks.sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file KLF5_UT_S23_acrossSeqRuns \
--plot \
--log-output-file KLF5_UT_S23_acrossSeqRuns.log


cd /data/jane/KDcutNtag/Mar2021run/IgGsamples/
samtools merge -u /data/jane/KDcutNtag/Mar2021run/IgGsamples/IgG_IL13_merged.bam \
/data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments//STAR210303_144341/IgG_S7_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam \
/data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments/STAR210301_153239/IgG_S15_UTout/Aligned.sortedByCoord.duplicateRemoved.out.bam 


samtools merge -u IgG_UT_merged.bam KLF5_1075_UT_merged.bam KLF5_1332_UT_merged.bam KLF5_1043_UT_merged.bam

##let's just...get the IDR sets
cd /data/jane/KDcutNtag/Mar2021run/KLF5samples
#Sort and then run KLF5_UT IDR
sort -k8,8nr KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.narrowPeak \
> KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.narrowPeak \
> KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.narrowPeak \
> KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.sorted.narrowPeak

sort -k8,8nr KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.narrowPeak \
> KLF5_UT_mergedP0p01/KLF5_UT_mergedP0p01_peaks.sorted.narrowPeak

macs2 callpeak -t KLF5_UT_merged.bam -f BAM \
--keep-dup 1 -n KLF5_UT_mergedP0p01  -p 0.01 --nomodel --outdir KLF5_UT_mergedP0p01


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

#then find the intersection of all 3 IDR sets
#first find the significnat IDRs
awk '{if($5 >= 540) print $0}' KLF5_UT_1043_1075_acrossSeqRuns >KLF5_UT_1043_1075_acrossSeqRunsSig.txt
awk '{if($5 >= 540) print $0}' KLF5_UT_1043_1332_acrossSeqRuns >KLF5_UT_1043_1332_acrossSeqRunsSig.txt
awk '{if($5 >= 540) print $0}' KLF5_UT_1075_1332_acrossSeqRuns >KLF5_UT_1075_1332_acrossSeqRunsSig.txt

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

#checking the enhancer region
KLF510431332<-read.csv(file="KLF5_IL13_1075_1332_acrossSeqRunsSig.bed",sep="\t",header=F)
> KLF510751332<-read.csv(file="KLF5_IL13_1075_1332_acrossSeqRunsSig.bed",sep="\t",header=F)
> chr67532<-subset(KLF510751332,((KLF510751332$V1=="chr6"))
+ )
> a<-subset(chr67532,chr67532$V3<35586932)
> a<-a[order(a$V2),]

##let's do intersection
bedtools -intersection -a chr6EnhancerKD.bed

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


namesOfFiles=("KLF5f_S22_IL13out")
bamFile="Aligned.sortedByCoord.duplicateRemoved.out.bam.bai"
for file in ${namesOfFiles[@]}; do
echo $file
scp jane@10.37.30.76:/data/jane/KDcutNtag/Mar2021run/AlignAllOut/Alignments/STAR210301_153239/$file/$bamFile $file.bam.bai
done

###trying to use MANorm2
##can I plot hte MA plto somehow?
cd /data/jane/KDcutNtag/Mar2021run/KLF5samples
samtools view -h KLF5_1332_IL13_merged.bam >KLF5_1332_IL13_merged.sam
samtools view -h KLF5_1332_UT_merged.bam >KLF5_1332_UT_merged.sam
samtools view -h KLF5_1043_IL13_merged.bam >KLF5_1043_IL13_merged.sam
samtools view -h KLF5_1043_UT_merged.bam >KLF5_1043_UT_merged.sam
samtools view -h KLF5_1075_IL13_merged.bam >KLF5_1075_IL13_merged.sam
samt_peaks.narrowPeakools view -h KLF5_1075_UT_merged.bam >KLF5_1075_UT_merged.sam

sam2bed -i KLF5_1332_IL13_merged.sam -o KLF5_1332_IL13_merged.bed
sam2bed -i KLF5_1332_UT_merged.sam -o KLF5_1332_UT_merged.bed
sam2bed -i KLF5_1043_IL13_merged.sam -o KLF5_1043_IL13_merged.bed
sam2bed -i KLF5_1043_UT_merged.sam -o KLF5_1043_UT_merged.bed
sam2bed -i KLF5_1075_IL13_merged.sam -o KLF5_1075_IL13_merged.bed
sam2bed -i KLF5_1075_UT_merged.sam -o KLF5_1075_UT_merged.bed

KLF5_1332_IL13_merged.bam_q0p1/KLF5_1332_IL13_merged.bam_peaks.sorted.narrowPeak

##change the summits.bed

 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1332_IL13_merged.bam_peaks.narrowPeak >KLF5_1332_IL13_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.narrowPeak >KLF5_1332_UT_merged.bam_q0p1/KLF5_1332_UT_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.narrowPeak >KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1043_IL13_merged.bam_q0p1/KLF5_1043_IL13_merged.bam_peaks.narrowPeak >KLF5_1043_IL13_merged.bam_q0p1/KLF5_1043_IL13_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.narrowPeak >KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.narrowPeak.bed
perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.narrowPeak >KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.narrowPeak.bed

perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1332_UT_merged.bam_p0p01/KLF5_1332_IL13_merged.bam_peaks.narrowPeak >KLF5_1332_UT_merged.bam_p0p01/KLF5_1332_IL13_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1332_UT_merged.bam_p0p01/KLF5_1332_UT_merged.bam_peaks.narrowPeak >KLF5_1332_UT_merged.bam_p0p01/KLF5_1332_UT_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1043_UT_merged.bam_p0p01/KLF5_1043_UT_merged.bam_peaks.narrowPeak >KLF5_1043_UT_merged.bam_p0p01/KLF5_1043_UT_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1043_IL13_merged.bam_p0p01/KLF5_1043_IL13_merged.bam_peaks.narrowPeak >KLF5_1043_IL13_merged.bam_p0p01/KLF5_1043_IL13_merged.bam_peaks.narrowPeak.bed
 perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1075_UT_merged.bam_p0p01/KLF5_1075_UT_merged.bam_peaks.narrowPeak >KLF5_1075_UT_merged.bam_p0p01/KLF5_1075_UT_merged.bam_peaks.narrowPeak.bed
perl -ane 'print "$F[0]\t$F[1]\t$F[2]\n";'  <KLF5_1075_IL13_merged.bam_p0p01/KLF5_1075_IL13_merged.bam_peaks.narrowPeak >KLF5_1075_IL13_merged.bam_p0p01/KLF5_1075_IL13_merged.bam_peaks.narrowPeak.bed


profile_bins --peaks=KLF5_1332_IL13_merged.bam_p0p01/KLF5_1332_IL13_merged.bam_peaks.narrowPeak.bed,KLF5_1332_UT_merged.bam_p0p01/KLF5_1332_UT_merged.bam_peaks.narrowPeak.bed,KLF5_1043_IL13_merged.bam_p0p01/KLF5_1043_IL13_merged.bam_peaks.narrowPeak.bed,KLF5_1043_UT_merged.bam_p0p01/KLF5_1043_UT_merged.bam_peaks.narrowPeak.bed,KLF5_1075_IL13_merged.bam_p0p01/KLF5_1075_IL13_merged.bam_peaks.narrowPeak.bed,KLF5_1075_UT_merged.bam_p0p01/KLF5_1075_UT_merged.bam_peaks.narrowPeak.bed\
             --reads=KLF5_1332_IL13_merged.bed,KLF5_1332_UT_merged.bed,KLF5_1043_IL13_merged.bed,KLF5_1043_UT_merged.bed,KLF5_1075_IL13_merged.bed,KLF5_1075_UT_merged.bed\
             --labs=1332_IL13,1332_UT,1043_IL13,1043_UT,1075_IL13,1075_UT  -n DonorAll3Pvalue
             
profile_bins --peaks=KLF5_acrossRuns.merged.bed,KLF5_acrossRuns.merged.bed,KLF5_acrossRuns.merged.bed,KLF5_acrossRuns.merged.bed,KLF5_acrossRuns.merged.bed,KLF5_acrossRuns.merged.bed \
             --reads=KLF5_1332_IL13_merged.bed,KLF5_1332_UT_merged.bed,KLF5_1043_IL13_merged.bed,KLF5_1043_UT_merged.bed,KLF5_1075_IL13_merged.bed,KLF5_1075_UT_merged.bed\
             --labs=1332_IL13,1332_UT,1043_IL13,1043_UT,1075_IL13,1075_UT  -n DonorAll3Merged
             
profile_bins --peaks=KLF5_1043_IL13_merged.bam_p0p01/KLF5_1043_IL13_merged.bam_peaks.narrowPeak.bed,KLF5_1043_UT_merged.bam_q0p1/KLF5_1043_UT_merged.bam_peaks.narrowPeak.bed \
             --reads=KLF5_1043_IL13_merged.bed,KLF5_1043_UT_merged.bed \
             --labs=1043_IL13,1043_UT -n Donor1043

profile_bins --peaks=KLF5_1075_IL13_merged.bam_q0p1/KLF5_1075_IL13_merged.bam_peaks.narrowPeak.bed,KLF5_1075_UT_merged.bam_q0p1/KLF5_1075_UT_merged.bam_peaks.narrowPeak.bed \
             --reads=KLF5_1075_IL13_merged.bed,KLF5_1075_UT_merged.bed \
             --labs=1075_IL13,1075_UT -n Donor1043

##Apr2 analyses
#homer analyses on these
perl /data/jane/HOMER/configureHomer.pl -install homer
 PATH=$PATH:/data/jane/HOMER/bin/
mkdir KLF5_all_merged_HOMER
cd /data/jane/KDcutNtag/Mar2021run/KLF5samples/


perl -ane '$count++; print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n";' </data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_acrossRuns.merged.bed >/data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_acrossRuns.merged.homer.bed

perl ./findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_acrossRuns.merged.homer.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_all_merged_HOMER -size given

perl ./findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_acrossRuns.merged.homer.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_all_merged_HOMER -size given -bg /data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/nonSig2.bed

perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/resDESEQ2SigIL13enriched.txt >/data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/resDESEQ2SigIL13enriched.bed

perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/nonSig.bed >/data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/nonSig2.bed

nohup perl ./bin/findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/resDESEQ2SigIL13enriched.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_UP_bgnonSig -size given -bg /data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/nonSig2.bed >nohupUPNonSig.txt &

nohup perl ./bin/findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/KLF5samples/STAT6_KLF5/resDESEQ2SigIL13enricheddown.sorted.Chr.txt \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/KLF5_DOWN_bgnonSig -size given -bg /data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/nonSig2.bed >nohupDownNonSig.txt &



perl ./bin//findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/KLF5samples/DESEQ2/resDESEQ2SigIL13enriched.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_IL13UPii -size given 

##runing homer on downregulated genes

perl ./bin//findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/KLF5samples/STAT6_KLF5/resDESEQ2SigIL13enricheddown.sorted.Chr.txt \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_IL13Down -size given


perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.txt >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.bed 
bedtools sort -i downregulatedDE.bed >downregulatedDE.sorted.bed
nohup perl ./bin//findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.sorted.bed hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_WalterPeaksDown -size given >walterPeaksDown.txt &


perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.txt >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.bed 
bedtools sort -i upregulatedDE.bed >upregulatedDE.sorted.bed
nohup perl ./bin//findMotifsGenome.pl /data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.sorted.bed hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_WalterPeaksUp -size given >walterPeaksUp.txt &

###
perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.txt >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.bed 
bedtools sort -i PeakOnlyGreaterp01.bed >PeakOnlyGreaterp01.sorted.bed

nohup perl ./bin//findMotifsGenome.pl \
/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.sorted.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_WalterPeaksUpBG -size given >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/walterPeaksUpBG.txt \
-bg /data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.sorted.bed &

nohup perl ./bin//findMotifsGenome.pl \
/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.sorted.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_WalterPeaksDownBG -size given >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/walterPeaksDownBG.txt \
-bg /data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.sorted.bed &
