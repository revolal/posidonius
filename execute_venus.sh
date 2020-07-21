#!/bin/sh
python setup.py install --user
cargo install --path . --force
python cases/venus_hot.py target/venus_hot.json
posidonius start target/venus_hot.json target/venus_hot.bin target/venus_hot_history.bin 
