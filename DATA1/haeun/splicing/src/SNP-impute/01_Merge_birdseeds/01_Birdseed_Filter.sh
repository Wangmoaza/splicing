
source 00_source

NAME=$1  # "TCGA-D8-A27L"
INPUT_BIRD=${INPUT_DIR}/$1


echo -e "probeset_id\tCall" > $TMP_DIR/$NAME.input
awk '{if ($3<0.1) print $1"\t"$2}' $INPUT_BIRD >> $TMP_DIR/$NAME.input
$APT_DIR/apt-format-result --calls-file $TMP_DIR/$NAME.input --annotation-file $AFFY_DB --export-plinkt-file $TMP_DIR/$NAME.out

