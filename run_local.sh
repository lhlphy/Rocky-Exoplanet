#!/bin/bash
# test script for wsl, test project on local machine
# this script is for wsl, not for windows
echo "Maxwell" 
export lavatype=zero
# Read parameters from params.txt
PARAMS=$(cat params.txt)
# Extract the value of --id parameter and remove any trailing whitespace
ID=$(grep -oP '(?<=--id )\d+' params.txt | tr -d '[:space:]')
echo "ID: $ID"
echo "PARAMS: $PARAMS"
# Create the directory temp/Rn
mkdir -p temp/R$ID  
# Copy params.txt to temp/Rn
cp params.txt temp/R$ID
# Run the Python script with the parameters
# python lib/Full_spectrum_main.py $PARAMS
echo "DONE"