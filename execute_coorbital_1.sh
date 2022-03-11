#!/bin/sh
python setup.py install --user
cargo install --path . --force
if [ -e target/co-orbital.json ]
then
    rm ./target/co-orbital.json
else 
    echo "no jason let's go "
fi
if [ -e target/co-orbital.bin ]
then
    rm ./target/co-orbital.bin
else
    echo "no bin let's go"
fi
if [ -e target/co-orbital_history.bin ]
then
    rm ./target/co-orbital_history.bin
else
    echo "no history let's go"
fi

# python cases/example_kaula.py target/kaula_test.json
# posidonius start target/kaula_test.json target/kaula_test.bin target/kaula_test_history.bin 
python cases/co-orbital/cas5_-5.46_helio_bary.py target/co-orbital.json
posidonius start target/co-orbital.json target/co-orbital.bin target/co-orbital_history.bin 

