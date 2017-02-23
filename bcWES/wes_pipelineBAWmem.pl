#!/usr/bin/env perl

# @Author: Lijia Yu <lijia>
# @Date:   2017-02-15 12:39 -06:00
# @Email:  yu@lijiayu.net
# @Last modified by:   lijia
# @Last modified time: 2017-02-21 09:53 -06:00

# Purpose: Whole Exome sequencing calling SNV/Indel pipeline
# require tools: GATK, samtools, fastqc, picard, trimmomatic, bwa
#

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

#Usage part

my $scpl_VER = '0.1.0.170215.Beta';
my (
    $opt_help, $opt_man, $opt_versions,
    $IN_DIR,  $OUT_DIR,          $TMP_DIR,     $THREADS,
    $MAX_MEM, $GENOME_REFERENCE, $EXON_REGION, $PREFIX,
    $FASTQ1,  $FASTQ2
);

GetOptions(
    'help' => \$opt_help,
    'man'  => \$opt_man,
    'versions'      => \$opt_versions,
    'input_dir=s'   => \$IN_DIR,
    'output_dir=s'  => \$OUT_DIR,
    'tmp_dir=s'     => \$TMP_DIR,
    'thread_num=i'  => \$THREADS,
    'max_mem=s'     => \$MAX_MEM,
    'genome_ref=s'  => \$GENOME_REFERENCE,
    'exon_region=s' => \$EXON_REGION,
    'prefix=s'      => \$PREFIX,
    'forward_fq=s'  => \$FASTQ1,
    'reverse_fq=s'  => \$FASTQ2,

) or pod2usage( -verbose => 1 ) && exit;

if ( defined $opt_versions ) {
    print
      "\nModules, Perl, OS, Program info:\n",
      "  Pod::Usage            $Pod::Usage::VERSION\n",
      "  Getopt::Long          $Getopt::Long::VERSION\n",
      "  strict                $strict::VERSION\n",
      "  Perl                  $]\n",
      "  OS                    $^O\n", "  $0         $scpl_VER\n", "\n\n"
      and exit;
}

pod2usage( -verbose => 1 ) && exit if defined $opt_help;
pod2usage( -verbose => 2 ) && exit if defined $opt_man;


pod2usage(
    -verbose  => 99,
    -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
  )
  && exit
  if !defined $IN_DIR;
pod2usage(
    -verbose  => 99,
    -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
  )
  && exit
  if !defined $OUT_DIR;
pod2usage(
    -verbose  => 99,
    -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
  )
  && exit
  if !defined $EXON_REGION;
pod2usage(
    -verbose  => 99,
    -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
  )
  && exit
  if !defined $PREFIX;
pod2usage(
    -verbose  => 99,
    -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
  )
  && exit
  if !defined $FASTQ1;
pod2usage(
    -verbose  => 99,
    -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
  )
  && exit
  if !defined $FASTQ2;

#end Usage part

=head1 NAME

 wes_pipeline.pl Whole Exome Sequencing bulk data pipeline

=head1 SYNOPSIS

usage: wes_pipeline.pl [--help] [--man] [--opt_name] [--versions]
               [--input_dir INPUT_DIRECTROY] [--output_dir OUTPUT_DIRECTROY]
               [--tmp_dir TMP_DIRECTROY_FOR_JAVA] [--thread_num THREADS][--max_mem MAX_MEM]
               [--genome_ref GENOME_REFERENCE] [--exon_region EXON_REGION] [--prefix PREFIX]
               [--forward_fq FASTQ1] [--reverse_fq FASTQ2]



=head1 DESCRIPTION

 Bulk data Whole Exome sequencing pipeline

 You need to provide input file, output file name at least.

 For help write:
 perl wes_pipeline.pl --help
 perl wes_pipeline.pl --man
 perl wes_pipeline.pl --opt_name


=head1 ARGUMENTS

 --help      print Options and Arguments instead of splitting fasta
 --man       print complete man page instead of splitting fasta
 --opt_name  print name of the script and synopsis

=head1 OPTIONS

 --versions               print Modules, Perl, OS, Program info
 --input_dir              input fastq file directory
 --output_dir             output results directory
 --tmp_dir                tmp directory used by Java
 --thread_num             thread number (default is 4 threads)
 --max_mem                maximun memory used by Java (default is 4G )
 --genome_ref             Reference genome (default is /home/lijia/projects/zheng/genome/ucsc.hg19.fasta)
 --exon_region            exon sequencing target regions bed file (example: /home/lijia/projects/zheng/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed)
 --prefix                 prefix of output file name
 --forward_fq             input fastq data, forward strand
 --reverse_fq             input fastq data, reverse strand

=head1 TESTED

  Modules, Perl, OS, Program info:
  Pod::Usage            1.63
  Getopt::Long          2.4
  strict                1.07
  Perl                  5.016003
  OS                    linux
  wes_pipeline.pl               0.1.0.170215.Beta


=head1 BUGS

none


=head1 TODO

 I have no idea.


=head1 UPDATES

 2015.11.30
 Initial working code

 2017.02.15
 rewrite the pipeline, it will runs on choir cluster at UAB

=head1 EXAMPLE


wes_pipeline.pl  --input_dir=./input --output_dir=./output --tmp_dir=./tmp
         --thread_num=4 --max_mem=4G --exon_region=./exon_region.bed  --prefix=sample1
         --forward_fq=demo_R1.fastq.gz  --reverse_fq=demo_R2.fastq.gz

=cut

if ( !-d $IN_DIR ) {
    print "**ERROR: I cound't find the input directory!**\n\n";
    pod2usage(
        -verbose  => 99,
        -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
    ) and exit;
}

if ( !-d $OUT_DIR ) {
    print "Cound't find the output directory, I will creat it right now.\n";
    unless ( mkdir $OUT_DIR ) {
        die "Unable to create $OUT_DIR\n";
    }
}

if ( !defined $TMP_DIR ) {
    $TMP_DIR = "/tmp";
}
else {
    if ( !-d $TMP_DIR ) {
        print "Cound't find the temp directory, I will creat it right now.\n";
        unless ( mkdir $TMP_DIR ) {
            die "Unable to create $TMP_DIR\n";
        }
    }
}

if ( !defined $THREADS ) {
    $THREADS = 4;
}

if ( !defined $MAX_MEM ) {
    $MAX_MEM = "4G";
}

if ( !defined $GENOME_REFERENCE ) {
    $GENOME_REFERENCE = "/home/lijia/projects/zheng/genome/ucsc.hg19.fasta";
}
else {
    if ( !-f $GENOME_REFERENCE ) {
        print "I couldn't find the genome reference: $GENOME_REFERENCE!";
        pod2usage(
            -verbose  => 99,
            -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
        ) && exit;
    }
}

if ( !-f $EXON_REGION ) {
    print "I couldn't find the bed file: $EXON_REGION!";
    pod2usage(
        -verbose  => 99,
        -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
    ) && exit;
}

if ( !-f "$IN_DIR/$FASTQ1" ) {
    print "I couldn't find the fastq file: $IN_DIR/$FASTQ1!";
    pod2usage(
        -verbose  => 99,
        -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
    ) && exit;
}

if ( !-f "$IN_DIR/$FASTQ2" ) {
    print "I couldn't find the fastq file: $IN_DIR/$FASTQ2!";
    pod2usage(
        -verbose  => 99,
        -sections => [qw(NAME SYNOPSIS DESCRIPTION EXAMPLE)]
    ) && exit;
}

( my $FQ1 = $FASTQ1 ) =~ s/\.[^.]+$//;
( my $FQ2 = $FASTQ2 ) =~ s/\.[^.]+$//;

my $HAPMAP = "/home/lijia/projects/zheng/database/hapmap_3.3.hg19.sites.vcf";
my $OMNI   = "/home/lijia/projects/zheng/database/1000G_omni2.5.hg19.sites.vcf";
my $DBSNP  = "/home/lijia/projects/zheng/database/dbsnp_138.hg19.vcf";
my $MILLS =
"/home/lijia/projects/zheng/database/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
my $G1000_INDEL =
  "/home/lijia/projects/zheng/database/1000G_phase1.indels.hg19.sites.vcf";
my $G1000_SNP =
"/home/lijia/projects/zheng/database/1000G_phase1.snps.high_confidence.hg19.sites.vcf";

my $PICARD   = "/home/lijia/bin/picard.jar";
my $GATK     = "/home/lijia/bin/GenomeAnalysisTK.jar";
my $TRIM     = "/home/lijia/bin/trimmomatic-0.36.jar";
my $JAVA     = "/usr/bin/java";
my $FASTQC   = "/home/lijia/bin/fastqc";
my $BWA      = "/home/lijia/bin/bwa";
my $SAMTOOLS = "/home/lijia/bin/samtools";

## Step 1 FASTQC quality control

mkdir "$OUT_DIR/fastqc";

mkdir "$OUT_DIR/trimed";

!system(
"$FASTQC -t $THREADS -o $OUT_DIR/fastqc --noextract -d $TMP_DIR -f fastq $IN_DIR/$FASTQ1 $IN_DIR/$FASTQ2"
) or die "fastqc is wrong \n";

## Step 2 Trim  sequences

!system(
"$JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $TRIM PE -threads 20 -phred33 $IN_DIR/$FASTQ1 $IN_DIR/$FASTQ2 $OUT_DIR/trimed/$FQ1.trimed.gz $OUT_DIR/trimed/$FQ1.forward_unpaired.fq.gz $OUT_DIR/trimed/$FQ2.trimed.gz  $OUT_DIR/trimed/$FQ2.reverse_unpaired.fq.gz SLIDINGWINDOW:4:15 MINLEN:36 LEADING:3 TRAILING:3"
) or die "trime R1.fq.gz is wrong\n";

## Step 3 BWA pair-end reads alignment (MEM method)

=begin comment
BWA-MEM with Picard markDuplicates
The BWA-MEM algorithm performs local alignment.
It may produce multiple primary alignments for different part of a query sequence.
This is a crucial feature for long sequences.
However, some tools such as Picard’s markDuplicates does not work with split alignments.
One may consider to use option -M to flag shorter split hits as secondary.
=end comment

=cut

!system("$BWA mem -M -t $THREADS  $GENOME_REFERENCE $OUT_DIR/trimed/$FQ1.trimed.gz  $OUT_DIR/trimed/$FQ2.trimed.gz  | $SAMTOOLS view -q 10 -h - > $OUT_DIR/$PREFIX.sam"
) or die "system bwa mem failed: $?";

## Step 4 sort file


!system("$JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD SortSam I=$OUT_DIR/$PREFIX.sam O=/dev/stdout SO=coordinate VALIDATION_STRINGENCY=SILENT | $JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT RGLB=hg19 RGPL=illumina RGPU=hg19PU RGSM=$PREFIX I=/dev/stdin O=$OUT_DIR/$PREFIX.addheader.bam") or die "picard sort and add header has problem\n";


!system("/usr/bin/rm -rf $OUT_DIR/*sam") or die "couldn't remove sam";
!system("/usr/bin/rm -rf $OUT_DIR/trimed") or die "couldn't remove trimed fastq";

=begin comment
!system(
"$JAVA  -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD SortSam I=$OUT_DIR/$PREFIX.sam O=$OUT_DIR/$PREFIX.sorted.bam SO=coordinate VALIDATION_STRINGENCY=SILENT"
) or die "picard sort has problem\n";

## Step 5 add header

!system(
    "$JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR  -jar $PICARD AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT RGLB=hg19 RGPL=illumina RGPU=hg19PU RGSM=$PREFIX I=$OUT_DIR/$PREFIX.sorted.bam O=$OUT_DIR/$PREFIX.addheader.bam"
) or die "picard add header has problem\n";





Q: Can Picard tools read from stdin and write to stdout?

A: Most Picard tools can do this.
To read from stdin, specify /dev/stdin as the input file.
To write to stdout, specify /dev/stdout as the output file,
and also add the option QUIET=true to suppress other messages that otherwise would be written to stdout.
Not that Picard determines whether to write in SAM or BAM format by examining the file extension of the output file.
Since /dev/stdout ends in neither .sam nor .bam, Picard defaults to writing in BAM format.
Some Picard tools, e.g. MarkDuplicates, cannot read their input from stdin,
because it makes multiple passes over the input file.
When writing a BAM to stdout so that it can be read by another program,
passing the argument COMPRESSION_LEVEL=0 to the program writing the BAM to stdout can reduce unnecessary computation.
=end comment
=cut


## Step 6 remove duplication by Picard

!system("$JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD MarkDuplicates I=$OUT_DIR/$PREFIX.addheader.bam O=$OUT_DIR/$PREFIX.sortedRMdup.bam M=$OUT_DIR/$PREFIX.marked_dup_metrics.txt"
) or die "Picard: remove duplication has problem \n";
#"$SAMTOOLS rmdup $OUT_DIR/$PREFIX.addheader.bam $OUT_DIR/$PREFIX.sortedRMdup.bam"

## Step 7 build bam index by Picard

!system("$JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $PICARD BuildBamIndex I=$OUT_DIR/$PREFIX.sortedRMdup.bam")
  or die "picard: index bam file has problem \n";

## Step 8 Define intervals to target for local realignment

!system(
    "$JAVA -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME_REFERENCE -T RealignerTargetCreator -L $EXON_REGION -I $OUT_DIR/$PREFIX.sortedRMdup.bam -o $OUT_DIR/$PREFIX.RealignerTargetCreator.intervals  --known $G1000_INDEL --known $MILLS --interval_padding 100"
) or die "system GATK RealignerTargetCreator failed: $?";

## Step 9 Perform local realignment of reads around indels

!system(
    "java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME_REFERENCE -T IndelRealigner -L $EXON_REGION -I $OUT_DIR/$PREFIX.sortedRMdup.bam -targetIntervals $OUT_DIR/$PREFIX.RealignerTargetCreator.intervals -o $OUT_DIR/$PREFIX.Realigner.bam -known $G1000_INDEL -known $MILLS --interval_padding 100"
) or die "system GATK IndelRealigner failed: $?";

## Step 10 Detect systematic errors in base quality scores

!system(
"java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME_REFERENCE -T BaseRecalibrator -nct $THREADS -L $EXON_REGION -I $OUT_DIR/$PREFIX.Realigner.bam -knownSites $DBSNP -knownSites $G1000_INDEL -knownSites $MILLS -o $OUT_DIR/$PREFIX.recal.grp --interval_padding 100"
) or die "system GATK BaseRecalibrator step 1 failed: $?";

## Step 11  Base quality score recalibration

!system(
"java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME_REFERENCE -T BaseRecalibrator -nct $THREADS -L $EXON_REGION -I $OUT_DIR/$PREFIX.Realigner.bam -BQSR $OUT_DIR/$PREFIX.recal.grp -o $OUT_DIR/$PREFIX.post_recal.grp -knownSites $DBSNP -knownSites $G1000_INDEL -knownSites $MILLS --interval_padding 100"
) or die "system GATK BaseRecalibrator step 2 failed: $?";

## Step 12 After Base quality score recalibration, write out sequence read data

!system(
"java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME_REFERENCE -T PrintReads -nct $THREADS -I $OUT_DIR/$PREFIX.Realigner.bam -BQSR $OUT_DIR/$PREFIX.post_recal.grp -o $OUT_DIR/$PREFIX.recal.bam"
) or die "system GATK PrintReads failed: $?";

## Step 13 SNP/INDEL calling

!system(
"java -Xmx$MAX_MEM -Djava.io.tmpdir=$TMP_DIR -jar $GATK -R $GENOME_REFERENCE -T HaplotypeCaller -I $OUT_DIR/$PREFIX.recal.bam -D $DBSNP -mbq 20 -stand_call_conf 30.00 -o $OUT_DIR/$PREFIX.raw.vcf"
) or die "system GATK HaplotypeCaller failed: $?";

## Step 14 remove files

!system("/usr/bin/rm -rf $OUT_DIR/*sort*.ba*") or die "couldn't remove sort.bam";

!system("/usr/bin/rm -rf $OUT_DIR/*add*.ba*") or die "couldn't remove addheader.bam";

!system("/usr/bin/rm -rf $OUT_DIR/*Realigner.ba*") or die "couldn't remove Realigner.bam";

!system("/usr/bin/rm -rf $OUT_DIR/*raw.[rt]*  $OUT_DIR/*.grp $OUT_DIR/*.intervals $OUT_DIR/*.txt") or die "couldn't remove intermediate file";
