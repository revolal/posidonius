#!/bin/bash
echo "Into Launch-create "

for i in {0..2}
#let "i=1"
do
    python create_PLANETi.py target/test$i/trappist1_multi.json target/test$i/trappist1_multi_history.bin
    echo "sleep for 1s into launch "
    sleep 1.0
    echo ""
done