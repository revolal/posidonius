#!/bin/sh
python setup.py install --user
cargo install --path . --force
if [ -e target/trappist_1e.json ]
then
    rm ./target/trappist_1e.json
    rm ./target/trappist_1e.bin
    rm ./target/trappist_1e_history.bin
else
    echo "Let's go"
fi

python cases/trappist_1e.py target/trappist_1e.json
posidonius start target/trappist_1e.json target/trappist_1e.bin target/trappist_1e_history.bin 
