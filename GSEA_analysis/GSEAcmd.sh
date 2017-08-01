# @Author: Lijia Yu <lijia>
# @Date:   2017-07-14 13:50 -05:00
# @Email:  yu@lijiayu.net
# @Last modified by:   lijia
# @Last modified time: 2017-07-14 15:33 -05:00



#!/bin/bash



command -v java  >/dev/null 2>&1 || { echo >&2 "GSEAcmd requires java but it's not installed.  Aborting."; exit 1; }


usage(){
cat << EOF

Program: GSEAcmd (GSEA preranked gene enrichment analysis)
Version: 0.1.0.170714.Beta
Contact: Lijia Yu <yu@lijiayu.net>
StartDate: 2017.07.14


Before use this tool, please make sure you have the correct input files:
(1) Download MSigDB database from http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C3
(2) Generate Ranked list file http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#RNK:_Ranked_list_file_format_.28.2A.rnk.29

GSEA preranked gene enrichment analysis:
(3) preranked enrichment analysis.
(4) extract the Term with ES and NES value (STDOUT).


usage: ${0} [-h] [-r rnk file] [-m MAX_MEM] [-d MSigDB database] [-o OUT_DIR] [-p prefix] [-g gsea_software]

Example:
${0} -m 4G -r data/test.rnk -d database/c2.cp.reactome.v6.0.symbols.gmt -p ouput_test  -o /home/lijia/gsea_out -g /home/lijia/bin/gsea2-2.2.4.jar

Options:
        -h, --help                      show this help message and exit.
        -m  MAX_MEM                     max memory usage [default 1G].
        -r  RNK                         Ranked list file [.rnk].
        -d  MSigDB                      MSigDB database [.gmt].
        -o  OUT_DIR                     output folder.
        -p  PREFIX                      prefix of output files.
        -g  GSEA                        GSEA software path.
EOF
}


MAX_MEM="1G"

while getopts ":m:p:r:d:o:g:" opt;
do
        case "$opt" in
                m) MAX_MEM=$OPTARG;;
                r) RNK=$OPTARG;;
                d) MSigDB=$OPTARG;;
                p) PREFIX=$OPTARG;;
                o) OUT_DIR=$OPTARG;;
                g) GSEA=$OPTARG;;
                \?) usage
                exit 1
                ;;
        esac
done

if [ $# -lt 5 ] ; then
   usage
   echo "error: too few arguments"
   exit 1
fi

if [ ! -f $GSEA ]; then
        usage
    echo "error: '$GSEA' not exists. Please give me the GSEA software with path";
        exit 1
fi

if [ ! -f $RNK ]; then
        usage
    echo "error: '$RNK' not exists.";
        exit 1
fi


if [ ! -f $MSigDB ]; then
        usage
    echo "error: '$MSigDB' not exists.";
        exit 1
fi


if [[ ! -e $OUT_DIR ]]; then
    mkdir $OUT_DIR
elif [[ ! -d $OUT_DIR ]]; then
    echo "$OUT_DIR already exists, all files will output at here.";
fi


echo "GASE analysis start" > $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log


java  -Xmx$MAX_MEM  -cp $GSEA xtools.gsea.GseaPreranked \
      -gmx $MSigDB -collapse false -mode Max_probe -norm meandiv -nperm 1000 \
      -rnk $RNK -scoring_scheme weighted -rpt_label $PREFIX \
      -include_only_symbols true -make_sets true -plot_top_x 20 \
      -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
      -out $OUT_DIR -gui false


echo "GASE analysis start" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

echo "Generate enrichment term table start" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log

for i in `ls $OUT_DIR/$PREFIX*/gsea_report_for*.xls`;
do
  tail -n +2 $i | cut -f 1,4,5,6 >> $OUT_DIR/$PREFIX.tsv.tmp 
done

echo -e "NAME\tSIZE\tES\tNES" > $OUT_DIR/$PREFIX.tsv
cat $OUT_DIR/$PREFIX.tsv.tmp | sort -nrk 4 >> $OUT_DIR/$PREFIX.tsv

rm -rf $OUT_DIR/$PREFIX.tsv.tmp

echo "Generate enrichment term table end" >> $OUT_DIR/$PREFIX.log
date >> $OUT_DIR/$PREFIX.log
