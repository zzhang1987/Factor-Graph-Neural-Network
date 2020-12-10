#!/bin/bash

if [ ! -d "./dataset" ]; then
    mkdir dataset
fi 

python3 data_generate/ldpc.py --save_path dataset/ldpc_valid.pt --num 1000


