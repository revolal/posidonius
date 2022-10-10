#!/bin/sh
python setup.py install --user
cargo install --path . --force
# if [ -e target/kaula_test.json ]
# then
#     rm ./target/kaula_test.json
# else 
#     echo "no jason let's go "
# fi
# if [ -e target/kaula_test.bin ]
# then
#     rm ./target/kaula_test.bin
# else
#     echo "no bin let's go"
# fi
# if [ -e target/kaula_test_history.bin ]
# then
#     rm ./target/kaula_test_history.bin
# else
#     echo "no history let's go"
# fi

python cases/example_kaula.py target/kaula_test.json
posidonius start target/kaula_test.json target/kaula_test.bin target/kaula_test_history.bin 
# python cases/example_kaula.py target/ctl_test.json
# posidonius start target/ctl_test.json target/ctl_test.bin target/ctl_test_history.bin 

