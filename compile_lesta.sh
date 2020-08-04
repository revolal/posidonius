#!/bin/sh
python setup.py install --user
cargo install --path . --force
data_dir=/hpcstorage/revolal/Posidonius/T1e
python cases/trappist_1e.py $data_dir/trappist_1e.json
