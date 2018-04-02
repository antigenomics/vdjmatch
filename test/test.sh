#!/usr/bin/env bash

VDJDB="java -Xmx4G -jar `ls ../build/libs/vdjmatch-*.jar`"
RES="../src/test/resources"

$VDJDB -S human -R TRB --search-scope 3,1,3 -v-match -j-match $RES/sergey_anatolyevich.gz test
if [[ ! -s test.sergey_anatolyevich.txt ]]; then echo "No results file"; exit 1; fi; rm test.sergey_anatolyevich.annot.txt
if [[ ! -s test.annot.summary.txt ]]; then echo "No summary file"; exit 1; fi; cat test.annot.summary.txt | head; rm test.annot.summary.txt

$VDJDB -S human -R TRB --search-scope 3,0,3 -v-match -j-match --scoring-vdjmatch $RES/sergey_anatolyevich.gz test
if [[ ! -s test.sergey_anatolyevich.txt ]]; then echo "No results file"; exit 1; fi
if [[ ! -s test.annot.summary.txt ]]; then echo "No summary file"; exit 1; fi; cat test.annot.summary.txt | head

$VDJDB -S human -R TRB $RES/empty.txt test_empty
if [[ ! -s test_empty.empty.txt ]]; then echo "No results file"; exit 1; fi;
if [[ ! -s test_empty.annot.summary.txt ]]; then echo "No summary file"; exit 1; fi; cat test_empty.annot.summary.txt | head

$VDJDB -S human -R TRB $RES/no_found.txt test_empty2
if [[ ! -s test_empty2.no_found.txt ]]; then echo "No results file"; exit 1; fi;
if [[ ! -s test_empty2.annot.summary.txt ]]; then echo "No summary file"; exit 1; fi; cat test_empty2.annot.summary.txt | head

# Test compatibility with VDJtools
java -Xmx4G -jar ../vdjtools/target/vdjtools-*.jar Annotate test.sergey_anatolyevich.annot.txt annot/
if [[ ! -s annot/metadata.txt ]]; then echo "No metadata file"; exit 1; fi