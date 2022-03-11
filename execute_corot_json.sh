#!/bin/sh
python setup.py install --user
cargo install --path . --force
python cases/corot-7-bc.py corot_7bc_tide_7b_tides_05.json
