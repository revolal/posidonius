#!/bin/sh
python setup.py install --user
cargo install --path . --force
python cases/trappist1_multi.py target/trappist1_multi.json
posidonius start target/trappist1_multi.json target/trappist1_multi.bin target/trappist1_multi_history.bin 
