##install homer and set path
perl ./configureHomer.pl -install homer
PATH=$PATH:/data/jane/HOMER/.//bin/



#make into the right bed file format and sort
perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.txt >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.bed 
bedtools sort -i downregulatedDE.bed >downregulatedDE.sorted.bed

perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.txt >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.bed 
bedtools sort -i upregulatedDE.bed >upregulatedDE.sorted.bed

###get the background
perl -ane '$count++;print "chr$F[0]\t$F[1]\t$F[2]\tPeak$count\n"' </data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.txt >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.bed 
bedtools sort -i PeakOnlyGreaterp01.bed >PeakOnlyGreaterp01.sorted.bed


#run homer with the peaks and the background from above
nohup perl ./bin//findMotifsGenome.pl \
/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/upregulatedDE.sorted.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_WalterPeaksUpBG -size given >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/walterPeaksUpBG.txt \
-bg /data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.sorted.bed &

nohup perl ./bin//findMotifsGenome.pl \
/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/downregulatedDE.sorted.bed \
hg38 /data/jane/KDcutNtag/Mar2021run/KLF5samples/HOMER_WalterPeaksDownBG -size given >/data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/walterPeaksDownBG.txt \
-bg /data/jane/KDcutNtag/Mar2021run/WalterDEpeaks/PeakOnlyGreaterp01.sorted.bed &
