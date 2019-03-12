#!/usr/bin/env bash

FILENAME=$1
TIMES=$2

awk '{if ($1 ~ ">") {printf "\n%s\n", $1} else {print $0}}' $FILENAME | sed '/^$/d' > tmp
grep ">" tmp > tmp.grep
awk '{if ($1 ~ ">") {printf "%s %s\n", NR-1, $1}}' tmp.grep > tmp.awk
FILE_LEN=$(cat tmp.awk | wc -l)

> headers.grep
for i in $(seq 1 $TIMES);
do
	RAND_NUM=$RANDOM
	SELECTED_NUM=$(echo "$RAND_NUM % $FILE_LEN" | bc)
	grep "^$SELECTED_NUM " tmp.awk >> headers.grep
done

cut -d" " -f2 headers.grep > out; mv out headers.grep

> ${FILENAME%.*}.rand.fa
while read -r LINE;
do
	grep -A 1 "$LINE" tmp >> ${FILENAME%.*}.rand.${TIMES}.fa
done < headers.grep

rm tmp tmp.grep tmp.awk headers.grep

echo "Success!"
