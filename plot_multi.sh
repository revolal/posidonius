#!/bin/sh
python scripts/explore_history.py target/trappist1_multi.json   target/trappist1_multi_history.bin
python scripts/explore_history_multi.py target/trappist1_multi.json   target/trappist1_multi_history.bin
python scripts/explore_history_multi_orbital_freq.py target/trappist1_multi.json   target/trappist1_multi_history.bin
python scripts/explore_history_multi_semi_maj.py target/trappist1_multi.json   target/trappist1_multi_history.bin