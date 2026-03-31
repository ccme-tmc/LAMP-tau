#!/usr/bin/bash

num_freq=216
num_pt=10

python3 Extract_CFPs.py --name CoL2 -el perb/eq/stdout2 -nd $num_freq -np $num_pt -e
