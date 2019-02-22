#!/bin/bash
while read SRR filename
do
	#move file
	echo ln -s $SRR.fastq $filename.fastq
done < NewSRAs.txt

