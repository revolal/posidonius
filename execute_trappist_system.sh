#!/bin/sh
python setup.py install --user
cargo install --path . --force
python cases/trappist_system.py trappist_system.json

