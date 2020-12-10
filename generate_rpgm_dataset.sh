#!/bin/bash

python3 data_generate/generate_random_pgm.py --type hops --name "./synthetic_data/hop_train.dat " --hop True --size 900000
python3 data_generate/generate_random_pgm.py --type hops --name "./synthetic_data/hop_test.dat " --hop True --size 100000


python3 data_generate/generate_random_pgm.py --type raw --name "./synthetic_data/raw_train.dat " --hop True --size 90000
python3 data_generate/generate_random_pgm.py --type raw --name "./synthetic_data/raw_test.dat " --hop True --size 10000
