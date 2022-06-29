#!/bin/bash


# Using bash script to run SEDR. This is the only way around multiple
# CPU usage by pytorch. All other trials failed.


cd /home/pcnmartin/SEDR

files=$(ls /home/pcnmartin/Vesalius/Simulation/*.csv)

for i in $files
do
  echo $i
  python3 /home/pcnmartin/Vesalius/Simulation/SEDR_bash.py --file $i
done
