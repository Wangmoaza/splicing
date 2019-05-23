#!/bin/bash
# https://josephcckuo.wordpress.com/2016/11/17/modify-chromosome-notation-in-bam-file/

IN_DIR="/home/haeun/DATA2/splicing/CCLE_BRCA_bam/bamfiles.old"

dothis() {
	B=$1
	Bout=$2
	samtools view -h $B | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' \
		| samtools view -bS - > /home/haeun/DATA2/splicing/CCLE_BRCA_bam/bamfiles/$Bout
	samtools index /home/haeun/DATA2/splicing/CCLE_BRCA_bam/bamfiles/$Bout
}

export -f dothis
find $IN_DIR -type f | egrep "bam$" | \
	parallel --eta -k --noswap -j3 --load 80% dothis {} {/}


#for B in $DIR/*bam; 
#do
#    echo $B
#    samtools view -h $B | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > temp.bam
#    mv -f temp.bam $DIR/$B
#    samtools index $DIR/$B
#done
