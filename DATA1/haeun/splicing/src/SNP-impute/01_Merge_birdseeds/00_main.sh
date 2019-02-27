#!/bin/bash

source 00_source


mkdir $TMP_DIR
mkdir $OUTPUT_DIR

cat $LIST | \
#parallel --eta -j 24 --load 80% --noswap 'echo {}; bash 01_Birdseed_Filter.sh {}; bash 02_tped_recode.sh {}'
#bash 03_merge_tped.sh

parallel --eta -j 24 --load 80% --noswap 'echo {}; bash 01_Birdseed_Filter.sh {}'

#while read name
#do
#	echo $name
#	bash 01_Birdseed_Filter.sh $name
#	bash 02_tped_recode.sh $name
#done < $LIST

mv *.log logs/
mv *.errors logs/
