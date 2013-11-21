#!/bin/bash
mkdir scratch
mkdir scratch/Left
mkdir scratch/Right

for i in `seq 0 12`;
do
mkdir scratch/Left/site_$i
done

for i in `seq 1 13`;
do
mkdir scratch/Right/site_$i
done
