# @Author: Lijia Yu <lijia>
# @Date:   2017-02-14 13:59 -06:00
# @Email:  yu@lijiayu.net
# @Last modified by:   lijia
# @Last modified time: 2017-02-16 15:05 -06:00



#!/bin/bash

#source /etc/profile.d/modules.sh
#module load bwa/0.7.10
#module load samtools/1.3


# test if softwares installed
command -v bwa >/dev/null 2>&1 || { echo >&2 "bcWES requires bwa but it's not installed.  Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "bcWES requires samtools but it's not installed.  Aborting."; exit 1; }
command -v fastqc  >/dev/null 2>&1 || { echo >&2 "bcWES requires fastqc but it's not installed.  Aborting."; exit 1; }

# pass paramters
usage(){
cat << EOF

Program: bcWES (bulk cell WES SNV/Indel analysis pipeline)
Version: 0.1.0.151228.Beta
Contact: Lijia Yu <yu@lijiayu.net>
StartDate: 2015.12.25

Map and processing Bulk WES/WGS reads, calling SNV/Indel/CNV
(1) FASTQ quality control and quality statistics
(2) trim sequences from the reads
(3) map FASTQ1/FASTQ2 using BWA, reads realign with GATK
(4) call SNV/Indel with GATK

usage: ${0} [-h] [-t THREADS] [-m MAX_MEM] [-a FASTQ1] [-b FASTQ2] [-p PREFIX] [-r EXON_REGION] [-d TMP_DIR] [-o OUT_DIR]

Example:
${0} -t 4 -m 4G -a data/demo_R1.fastq.gz -b data/demo_R2.fastq.gz -p ouput_test -r /home/xuelab05/files/ex_region.sort.new1. -d /home/lijia/project/zheng/tmp -o /home/xuelab05/output
${0} -t 4 -m 4G -a XZ01252017/TTP-001_S1_L003_R1_001.fastq.gz -b XZ01252017/TTP-001_S1_L003_R2_001.fastq.gz -p TEST -r SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed -d tmp -o out

Options:
        -h, --help                      show this help message and exit.
        -t  THREADS                     threads [4].
        -m  MAX_MEM                     max memory usage [4G].
        -a  FASTQ1                      first mate of pair-end sequencing data [.fq/.fastq/.gz].
        -b  FASTQ2                      second mate of pair-end sequencing data [.fq/.fastq/.gz].
        -p  PREFIX                      prefix of output files.
        -r  EXON_REGION                 exon target region fils.
        -d  TMP_DIR                     temporary folder, which is needed by GATK/fastqc.
        -o  OUT_DIR                     output folder.
EOF
}

THREADS=2
MAX_MEM="4G"

while getopts ":t:m:a:b:p:r:d:o:" opt;
do
        case "$opt" in
                t) THREADS=$OPTARG;;
                m) MAX_MEM=$OPTARG;;
                a) FASTQ1=$OPTARG;;
                b) FASTQ2=$OPTARG;;
                p) PREFIX=$OPTARG;;
                r) EXON_REGION=$OPTARG;;
                d) TMP_DIR=$OPTARG;;
                o) OUT_DIR=$OPTARG;;
                \?) usage
                exit 1
                ;;
        esac
done


if [ $# -lt 12 ] ; then
   usage
   echo "error: too few arguments"
   exit 1
fi

re='^[0-9]+$'
if ! [[ $THREADS =~ $re ]] ; then
   echo "error: '$THREADS' Not a number" >&2;
   exit 1
fi

# check if input file exists
if [ ! -f $FASTQ1 ]; then
        usage
    echo "error: '$FASTQ1' not exists.";
        exit 1
fi

if [ ! -f $FASTQ2 ]; then
        usage
    echo "error: '$FASTQ2' not exists.";
        exit 1
fi


if [ ! -f $EXON_REGION ]; then
        usage
    echo "error: '$EXON_REGION' not exists.";
        exit 1
fi

if [ -z "$TMP_DIR" ]
then
   usage
   exit
fi

if [[ ! -e $OUT_DIR ]]; then
    mkdir $OUT_DIR
elif [[ ! -d $OUT_DIR ]]; then
    echo "$OUT_DIR already exists, all files will output at here.";
fi


#OUT_DIR=${OUT_DIR%/}

HAPMAP=/home/lijia/projects/zheng/database/hapmap_3.3.hg19.sites.vcf
OMNI=/home/lijia/projects/zheng/database/1000G_omni2.5.hg19.sites.vcf
DBSNP=/home/lijia/projects/zheng/database/dbsnp_138.hg19.vcf
MILLS=/home/lijia/projects/zheng/database/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
G1000_INDEL=/home/lijia/projects/zheng/database/1000G_phase1.indels.hg19.sites.vcf
G1000_SNP=/home/lijia/projects/zheng/database/1000G_phase1.snps.high_confidence.hg19.sites.vcf

PICARD=/home/lijia/bin/picard.jar
GENOME=/home/lijia/projects/zheng/genome/ucsc.hg19.fasta
GATK=/home/lijia/bin/GenomeAnalysisTK.jar
#TRIM=/home/lijia/project/zheng/project/perl/bin/clean_malbac_reads/trim_malbac_adapter_lowqual
TRIM=/home/lijia/bin/trimmomatic-0.36.jar

#M_primer=/home/lijia/project/zheng/project/perl/bin/clean_malbac_reads/MALBAC_primer.fa
#I_primer=/home/lijia/project/zheng/project/perl/bin/clean_malbac_reads/illumina_adapter.fa


echo "fastQC start" > $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

## (1)
fastqc -t $THREADS -o $OUT_DIR --noextract -d $TMP_DIR -f fastq $FASTQ1 $FASTQ2

echo "fastQC end" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

echo "trim start" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

## (2)

FQ1=`echo $FASTQ1 | sed s/.gz$//g`
FQ2=`echo $FASTQ2 | sed s/.gz$//g `

FQ1=`basename $FQ1`
FQ2=`basename $FQ2`


java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR  -jar $TRIM PE -threads 20 -phred33 $FASTQ1 $FASTQ2 $OUT_DIR/$FQ1.trimed.gz $OUT_DIR/$FQ1.forward_unpaired.fq.gz \
        $OUT_DIR/$FQ2.trimed.gz  $OUT_DIR/$FQ2.reverse_unpaired.fq.gz SLIDINGWINDOW:4:15 MINLEN:36 LEADING:3 TRAILING:3

##java -Xmx20g -jar $TRIM/trimmomatic-0.33.jar PE -threads 20 -phred33 $FASTQ1 $FASTQ2 $FQ1.forward_paired.fq.gz $FQ1.forward_unpaired.fq.gz \
##        $FQ2.reverse_paired.fq.gz $FQ2.reverse_unpaired.fq.gz SLIDINGWINDOW:4:15 MINLEN:36 LEADING:3 TRAILING:3

#$TRIM -t $THREADS -m $M_primer -i $I_primer $FASTQ1 $OUT_DIR/$FQ1.trimed.gz $OUT_DIR/$FQ1.trimed.stat
#$TRIM -t $THREADS -m $M_primer -i $I_primer $FASTQ2 $OUT_DIR/$FQ2.trimed.gz $OUT_DIR/$FQ2.trimed.stat

echo "trim end" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

echo "mapping start" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

## (3)


bwa aln -n 2 -t $THREADS $GENOME $OUT_DIR/$FQ1.trimed.gz > $OUT_DIR/$FQ1.sai
bwa aln -n 2 -t $THREADS $GENOME $OUT_DIR/$FQ2.trimed.gz > $OUT_DIR/$FQ2.sai

bwa sampe -a 800 $GENOME $OUT_DIR/$FQ1.sai $OUT_DIR/$FQ2.sai \
        $OUT_DIR/$FQ1.trimed.gz $OUT_DIR/$FQ2.trimed.gz |\
        samtools view -q 10 -h - > $OUT_DIR/$PREFIX.sam

java -jar -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD SortSam I=$OUT_DIR/$PREFIX.sam \
        O=$OUT_DIR/$PREFIX.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

java -jar -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT \
        RGLB=hg19 RGPL=illumina RGPU=hg19PU RGSM=$PREFIX \
        I=$OUT_DIR/$PREFIX.sorted.bam O=$OUT_DIR/$PREFIX.addheader.bam

samtools rmdup $OUT_DIR/$PREFIX.addheader.bam $OUT_DIR/$PREFIX.sortedRMdup.bam
samtools index $OUT_DIR/$PREFIX.sortedRMdup.bam


java -jar -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME -T RealignerTargetCreator \
        -L $EXON_REGION -I $OUT_DIR/$PREFIX.sortedRMdup.bam \
        -o $OUT_DIR/$PREFIX.RealignerTargetCreator.intervals \
        --known $G1000_INDEL \
        --known $MILLS \
        --interval_padding 100
      #  --fix_misencoded_quality_scores \


java -jar -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME -T IndelRealigner \
        -L $EXON_REGION -I $OUT_DIR/$PREFIX.sortedRMdup.bam \
        -targetIntervals $OUT_DIR/$PREFIX.RealignerTargetCreator.intervals \
        -o $OUT_DIR/$PREFIX.Realigner.bam \
        -known $G1000_INDEL \
        -known $MILLS \
        --interval_padding 100
        #--fix_misencoded_quality_scores \


java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME -T BaseRecalibrator -nct $THREADS \
        -L $EXON_REGION -I $OUT_DIR/$PREFIX.Realigner.bam \
        -knownSites $DBSNP \
        -knownSites $G1000_INDEL \
        -knownSites $MILLS \
        -o $OUT_DIR/$PREFIX.recal.grp \
        --interval_padding 100

java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME -T BaseRecalibrator -nct $THREADS \
        -L $EXON_REGION -I $OUT_DIR/$PREFIX.Realigner.bam \
        -BQSR $OUT_DIR/$PREFIX.recal.grp \
        -o $OUT_DIR/$PREFIX.post_recal.grp \
        -knownSites $DBSNP \
        -knownSites $G1000_INDEL \
        -knownSites $MILLS \
        --interval_padding 100

java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME -T PrintReads -nct $THREADS \
        -I $OUT_DIR/$PREFIX.Realigner.bam -BQSR $OUT_DIR/$PREFIX.post_recal.grp \
        -o $OUT_DIR/$PREFIX.recal.bam

echo "mapping end" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log


echo "GATK call SNV/Indel start" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

## (4)


java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME \
        -T HaplotypeCaller -L $EXON_REGION -I $OUT_DIR/$PREFIX.recal.bam \
        -D $DBSNP \
        -mbq 20 -stand_call_conf 30.00 \
        -o $OUT_DIR/$PREFIX.raw.vcf
#-stand_emit_conf 10.00 \
#
#  -T HaplotypeCaller -I $OUT_DIR/$PREFIX.recal.bam \


java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME \
        -T VariantRecalibrator -mode SNP -input $OUT_DIR/$PREFIX.raw.vcf \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
        -resource:ommi,known=false,training=true,truth=false,prior=12.0 $OMNI \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000_SNP \
        -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $DBSNP \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP --maxGaussians 4\
        -recalFile $OUT_DIR/$PREFIX.snp.raw.recal -tranchesFile $OUT_DIR/$PREFIX.snp.raw.tranches \
        --TStranche 95.0 --TStranche 96.0 --TStranche 97.0 --TStranche 98.0 --TStranche 99.0 --TStranche 99.5

java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME \
        -T VariantRecalibrator -mode INDEL -input $OUT_DIR/$PREFIX.raw.vcf \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS \
        -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $DBSNP \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000_SNP \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP --maxGaussians 4 \
        -recalFile $OUT_DIR/$PREFIX.indel.raw.recal -tranchesFile $OUT_DIR/$PREFIX.indel.raw.tranches


java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME \
        -T ApplyRecalibration -mode SNP -input $OUT_DIR/$PREFIX.raw.vcf \
        -tranchesFile $OUT_DIR/$PREFIX.snp.raw.tranches \
        -recalFile $OUT_DIR/$PREFIX.snp.raw.recal \
        -o $OUT_DIR/$PREFIX.snp.filtered.vcf --ts_filter_level 98

java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME \
        -T ApplyRecalibration -mode INDEL -input $OUT_DIR/$PREFIX.raw.vcf \
        -tranchesFile $OUT_DIR/$PREFIX.indel.raw.tranches \
        -recalFile $OUT_DIR/$PREFIX.indel.raw.recal \
        -o $OUT_DIR/$PREFIX.indel.filtered.vcf --ts_filter_level 97.0

echo "GATK call SNV/Indel end" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

rm -rf $OUT_DIR/*sort*bam*
rm -rf $OUT_DIR/*sam
rm -rf $OUT_DIR/*Realigner.bam*
rm -rf $OUT_DIR/*raw.[rt]*
rm -rf $OUT_DIR/*.grp
rm -rf $OUT_DIR/*.intervals

