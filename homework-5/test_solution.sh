#!/usr/bin/env bash

AB=$(date +%s%3N)
ACTUAL=$(echo 123 | ./student_submission)
AA=$(date +%s%3N)
DIFF1=$((AA - AB))

AB=$(date +%s%3N)
EXPECTED=$(echo 123 | ./sequential_implementation)
AA=$(date +%s%3N)
DIFF2=$((AA - AB))

if [ "$ACTUAL" = "$EXPECTED" ]; then echo "PASSED"; else echo  "FAILED"; fi

echo $DIFF1
echo $DIFF2