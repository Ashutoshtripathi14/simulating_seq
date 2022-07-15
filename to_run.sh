#!/usr/bin/env bash

#BSUB -J simulating_seq
#BSUB -q standard
#BSUB -W 6:10
#BSUB -e ~/LSFexample/simseq.err
#BSUB -o ~/LSFexample/simseq.out
#BSUB -n 1
#BSUB -M 1GB

module load python-3.9.0-gcc-9.3.0-5t75egs 
module load py-numpy-1.19.4-gcc-9.3.0-x2neh6p 
module load py-pandas-1.1.5-gcc-9.3.0-rjsya74

python ~/LSFexample/simulating_seq/script_to_sim.py 
