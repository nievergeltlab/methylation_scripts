#Split EWAS file into blocks for parallel processing
filename=cleaned_v3d #get a line count on the data to know how many chunks we need to make... 
probecount=$(wc -l "$filename".methylation | cut -d " " -f1)
#431032 probes
head -n1 "$filename".methylation > temporary_files/"$filename".header
#chunk defines number of probes in each file
CHUNK=2500
i=2
j=1
while [ $i -lt $probecount ]; do
echo $i $j
awk -v start=$i -v chunk=$CHUNK '{if (NR >= start  && NR <= start + chunk) print }' "$filename".methylation | cat temporary_files/"$filename".header -  > qced_blocks_v3d/"$filename"_chunk"$j".methylation

i=$(( $i + $CHUNK + 1))
j=$(($j +1))
done
