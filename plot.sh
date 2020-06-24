#!/bin/sh
python scripts/explore_history.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_history_1e.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_history_frequ.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_history_K2.py target/trappist_1e.json   target/trappist_1e_history.bin
