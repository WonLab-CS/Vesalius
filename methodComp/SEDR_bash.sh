#!/bin/bash

cd /home/pcnmartin/SEDR

files=$(ls /home/pcnmartin/Vesalius/Simulation/*.csv)

for i in $files
do
  echo $i
  python3 /home/pcnmartin/Vesalius/Simulation/SEDR_bash.py --file $i
done
