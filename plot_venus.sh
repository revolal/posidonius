#!/bin/sh
python scripts/explore_history.py target/venus_hot.json   target/venus_hot_history.bin
python scripts/explore_history_1e.py target/venus_hot.json   target/venus_hot_history.bin
python scripts/explore_history_orbital_freq.py target/venus_hot.json   target/venus_hot_history.bin
python scripts/explore_history_frequ.py target/venus_hot.json   target/venus_hot_history.bin
python scripts/explore_history_K2.py target/venus_hot.json   target/venus_hot_history.bin
