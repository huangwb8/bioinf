#!/bin/bash
# Program
#  sra2fq is used to convert .sra into .fastq.gz
# History:
#  2019/09/21 Weibin First release
# Parameters:
#  1: full path of a .txt file with first col as sample id, and the second col as case name
#  2: path of .sra space
#  3: path of saved space of .fastq.gz results 
#  4: numeric;the number of treads 
PATH=$PATH:~/bin
export PATH

## help information
echo "1: full path of a .txt file with first col as sample id, and the second col as case name"
echo "2: path of .sra space"
echo "3: path of saved space of .fastq.gz results"
echo "4: numeric;the number of treads"

## create space 
mkdir $3
touch $3/sra2fq.log

## read id
cat $1 | while read id # cat ./sra2case.txt | while read id
do
  arr=(${id})
  sample=${arr[0]}.sra
  case=${arr[1]}
			    
  ## -p show progress; -x print details; -3 writes single reads in special file;-N use row-id as name; -P print read-numbers; -f force to overwrite existing file(s);-e how many thread; --skip-technical skip technical reads
  time fasterq-dump -p -x -3 -N -P -f --skip-technical -e $4 $2/${sample} -O $3 >> $3/sra2fq.log 2>&1 
			          
  sed s/${sample}/${case}/ $3/${sample}_1.fastq > $3/${case}_1.fastq
  pigz -p $4 $3/${case}_1.fastq
				        
  sed s/${sample}/${case}/ $3/${sample}_2.fastq > $3/${case}_2.fastq
  pigz -p $4 -f $3/${case}_2.fastq
done
