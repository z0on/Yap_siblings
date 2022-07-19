cd /work2/01211/cmonstr/reyap


GENOME_REF=/work/03121/sbarfie/yap/Amil_v1.0/Amil_ref.fa
BAMS=/work/03121/sbarfie/yap/*.bam
mkdir yap_bams
cp /work/03121/sbarfie/yap/*.bam yap_bams/. &

ls $BAMS>bams
cat bams

# -----remapping to chromosome-level assembly

# extracting original reads with scores from recalibrated file

# batch job:
>unbam
for F in `cat bams`;do
OUT=`echo $F | perl -pe 's/\/.+\///' | perl -pe 's/\\..+/\.fq/'`
echo $OUT
echo "samtools view $F | awk '{print \"@\"\$1\"\n\"\$10\"\n+\n\"\$11}'  >$STOCKYARD/reyap/reads/$OUT" >>unbam;
done
ls5_launcher_creator.py -j unbam -n unbam -t 0:10:00 -e matz@utexas.edu -w 24 -a tagmap -q normal
sbatch unbam.slurm


# mapping, converting to bams, indexing
>maps
REF=$WORK/db/amilV2_chroms.fasta
for F in `ls reads/*.fq`; do
OUT=`echo $F | perl -pe 's/.+\///' | perl -pe 's/\\..+//'`
echo $OUT
echo "bowtie2 --no-unal --local -x $REF -U $F -S ${OUT}.sam && samtools sort -O bam -o ${OUT}.bam ${OUT}.sam && samtools index ${OUT}.bam">>maps
done
ls5_launcher_creator.py -j maps -n maps -a tagmap -e matz@utexas.edu -t 2:00:00 -w 24 -q normal
# mapsjob=$(sbatch --dependency=afterok:$mergejob  maps.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')
mapsjob=$(sbatch maps.slurm | grep "Submitted batch job" | perl -pe 's/\D//g')

grep overall maps.e*

# quality check
ls *[JA0-9].bam > bams.local
ls *nl.bam >bams
cat bams.local

export MinIndPerc=0.5
FILTERSQ='-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minInd $MI'
TODOQ="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"
echo 'export NIND=`cat bams.local | wc -l`; export MI=`echo "($NIND*$MinIndPerc+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams.local -r chr1:1-2000000 -GL 1 $FILTERSQ $TODOQ -P 12 -out dd && Rscript ~/bin/plotQC.R prefix=dd">a0
bash a0 &

cat quality.txt
ll -tr
wc -l bams.qc
# 291
wc -l bams
# 299

# ------ running ANGSD on each bam to get to calculate nsites per sample (and heterozygosity just because the code is this way)

cp bams.qc goodbams
>hets
>goodbams.allsites.het
REF=$STOCKYARD/db/amilV2_chroms.fasta
FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doMajorMinor 5 -setMaxDepth 50'
I=0
for F in `cat goodbams`; do
echo $I
echo "sleep $I && angsd -i $F -anc $REF $FILTERS -GL 1 -doSaf 1 -doCounts 1 -dumpCounts 3 -out ${F/.bam/} && realSFS ${F/.bam/}.saf.idx >${F/.bam/}.ml && awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml >>goodbams.allsites.het">>hets;
I=`echo "$I+0.1" |bc` # so lines are not overwritten
done
# check if it works (ctl-C the process if you see no errors immediately)
# head -1 hets | bash
# execute all lines in the file hets (code for TACC's stampede2)
ls6_launcher_creator.py -j hets -n hets -t 0:10:00 -e matz@utexas.edu -w 48 -a IBN21018 -q normal
sbatch hets.slurm


#---- which samples contain the least sites?
conda activate ReadProcessing
R
  zz=read.table("goodbams.allsites.het")
  zz=zz[order(zz$V2),]
  zz
  write(zz$V1[1:3],"3worst")
  q()

# ------- extracting sites shared between all three lowest-sites samples

REF=$STOCKYARD/db/amilV2_chroms.fasta
FILTERS='-minInd 3 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-3'
TODO="-doSaf 1 -doCounts 1 -anc $REF -ref $REF -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 1 -doGlf 2 -makeMatrix 1 -doIBS 1"
echo "angsd -b 3worst -GL 1 -P 6 $FILTERS $TODO -out 3w">3w
bash 3w
        # -> Total number of sites analyzed: 924850
        # -> Number of sites retained after filtering: 316484

# saving and indexing sites
zcat 3w.mafs.gz | cut -f 1,2 | tail -n +2 >goodsites
angsd sites index goodsites

#----------------------- for everybody's IBS

FILTERS='-minInd 236 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -minMaf 0.05 -snp_pval 1e-6 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO1='-doMajorMinor 1 -dosnpstat 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3 -makeMatrix 1 -doIBS 1'
echo "angsd -b bams.qc -sites goodsites -GL 1 $FILTERS $TODO1 -P 12 -out g3all">g3all
bash g3all

        # -> Total number of sites analyzed: 3317294
        # -> Number of sites retained after filtering: 14365
        # -> Number of sites retained after filtering: 3552 (with minmaf 0.05)
        
#----------------------- for only good bams

# edit goodba,s 

# glf3 for ngsRelate and IBS 
FILTERS='-minInd 236 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -minMaf 0.05 -snp_pval 1e-6 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO1='-doMajorMinor 1 -dosnpstat 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 3 -makeMatrix 1 -doIBS 1'
echo "angsd -b goodbams -sites goodsites -GL 1 $FILTERS $TODO1 -P 12 -out g3">g3
bash g3
        # -> Total number of sites analyzed: 3225791
        # -> Number of sites retained after filtering: 12400 

# ----- relatedness with NgsRelate
zcat g3.mafs.gz | cut -f5 |sed 1d >freq
ngsRelate -g g3.glf.gz -n 262 -f freq -O g3.relatedness

# ----- pcangsd

# rerunning with -dpGlf 2 (beagle output)
FILTERS='-minInd 236 -uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doHWE 1 -maxHetFreq 0.5 -snp_pval 1e-6 -sb_pval 1e-3 -hetbias_pval 1e-3'
TODO1='-doMajorMinor 1 -dosnpstat 1 -doMaf 1 -doCounts 1 -doPost 1 -doGlf 2 '
echo "angsd -b goodbams -sites goodsites -GL 1 $FILTERS $TODO1 -P 12 -out g2">g2
bash g2

pcangsd --beagle g2.beagle.gz --admix -o pcangsd.min --inbreedSamples -t 12 --minMaf 0.05 --selection --sites_save --maf_save

# ---- heterozygosity on min sites only

>hets
>goodbams.minsites.het
REF=$STOCKYARD/db/amilV2_chroms.fasta
FILTERS='-uniqueOnly 1 -skipTriallelic 1 -minMapQ 20 -minQ 20 -doMajorMinor 5 -setMaxDepth 50'
I=0
for F in `cat goodbams`; do
echo $I
echo "sleep $I && angsd -sites goodsites -i $F -anc $REF $FILTERS -GL 1 -doSaf 1 -doCounts 1 -dumpCounts 3 -out ${F/.bam/} && realSFS ${F/.bam/}.saf.idx >${F/.bam/}.ml && awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ${F/.bam/}.ml >>goodbams.minsites.het">>hets;
I=`echo "$I+0.1" |bc` # so lines are not overwritten
done
# check if it works (ctl-C the process if you see no errors immediately)
# head -1 hets | bash
# execute all lines in the file hets (code for TACC's stampede2)
ls6_launcher_creator.py -j hets -n hets -t 0:10:00 -e matz@utexas.edu -w 48 -a IBN21018 -q normal
sbatch hets.slurm

