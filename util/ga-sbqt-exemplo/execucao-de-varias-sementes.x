#!/bin/bash

STARTSEED=1
ENDSEED=5
PREVIOUSI=0


for i in $(seq $STARTSEED $ENDSEED); do
	SEEDOPTION="s/seed = "$PREVIOUSI"/seed = "$i"/g"
	sed -i "$SEEDOPTION" GaInput.txt
	./zerg.exe
	cp GaInput.txt GaInput.txt$i 
	mv ga-output.txt ga-output.txt$i
	mv operators-histogram.txt operators-histogram.txt$i

	PREVIOUSI=$i
done

#sed -i 's/original-string/new-string/g' file.txt




