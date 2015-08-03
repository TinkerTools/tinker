#!/bin/bash

#$ -S /bin/bash
#$ -N iEL-SCF
#$ -pe threaded 16
#$ -q all.q
#$ -cwd
#$ -o info.out
#$ -e error.out

printenv > env.txt

/home/aalbaugh/tinker/tinker_iEL-SCF/source/dynamic.x water512.xyz 1000000 1.0 10.0 1  -k water.key > test.log

