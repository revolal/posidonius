#!/bin/sh
python setup.py install --user
cargo install --path . --force
python cases/trappist_1e.py target/trappist_1e.json
posidonius start target/trappist_1e.json target/trappist_1e.bin target/trappist_1e_history.bin