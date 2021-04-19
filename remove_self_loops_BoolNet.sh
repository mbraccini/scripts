#!/bin/bash

filename="${1%.*}"
line=$(head -n 1 $1)
OUTPUT_FILENAME=$filename"_no_selfloop"
echo $line > $OUTPUT_FILENAME
nodes=$2
cat $1 | awk -v var=$nodes '{FS=","; OFS = FS} /Gene*/ {
	text=$2;
	gsub(/(!|\(|\)|&|\|)/," ",text); 
	max=split(text,array," "); 
   
		
    do {
		condition=0
		guess="Gene"int(rand()*var) + 1; 
	    for (i in array){
         	if (guess==array[i]){ 
				condition=1
			}
		}
	         
    } while (condition==1);
	
	gsub($1,guess,$2); 
	print $0
} ' >> $OUTPUT_FILENAME

