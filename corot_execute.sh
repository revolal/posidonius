#!/bin/sh
python setup.py install --user
cargo install --path . --force
if [ -e target/corot-7-bc.json ]
then
    rm ./target/corot-7-bc.json
    rm ./target/corot-7-bc.bin
    rm ./target/corot-7-bc_history.bin
else
    echo "Let's go"
fi

python cases/corot-7-bc.py target/corot-7-bc.json
posidonius start target/corot-7-bc.json target/corot-7-bc.bin target/corot-7-bc_history.bin 
