#!/bin/sh
python scripts/explore_history_1e.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_orbital_parameters_history.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_imaginaryk2_evolution.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_realk2_evolution.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_frequency_evolution.py target/trappist_1e.json   target/trappist_1e_history.bin
python scripts/explore_history_K2.py target/trappist_1e.json   target/trappist_1e_history.bin