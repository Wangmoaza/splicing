#!/bin/bash

source 00_source
NAME=$1
genderlist=$INPUT_DIR/list.gender.txt

declare -A dic
while read -r -a array
do
        dic["${array[0]}"]="${array[1]}"
done < $genderlist

echo -e $NAME"\t"$NAME"\t0\t0\t"${dic[$NAME]}"\t1" > $TMP_DIR/$NAME.out_filt.tfam  #FID IID fatherID motherID Gender(1=male/2=female) Phenotype(1=control/2=case)

cat <(cut -f2,5 $TMP_DIR/$NAME.out.tped) ../Affy6.0.map | awk 'BEGIN{FS=OFS="\t"} (NF==2){dic[$1]=$2} (NF>2){if ($2 in dic) print $1,$2,$3,$4,dic[$2]; else print $1,$2,$3,$4,"0 0"}' > $TMP_DIR/$NAME.out_filt.tped
ln -s ../Affy6.0.map $TMP_DIR/$NAME.out_filt.map

$PLINK_DIR/plink --noweb --tfile $TMP_DIR/$NAME.out_filt --recode --tab --out $OUTPUT_DIR/$NAME
