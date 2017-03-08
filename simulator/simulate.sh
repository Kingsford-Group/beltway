#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "illegal number of parameters"
    echo "USAGE: ./uniqWeights_simulate.sh <expNumber>"
    exit 1
fi
#int_simul1 --> integer weights 
outDirInt=data/int_simul$1
outDirFloat=data/float_simul$1

mkdir -p $outDirInt
java -jar build/GenerateCharacterSet.jar 50 50 400 charset.txt $outDirInt N Y
seq 5 1 20 | xargs -I{} java -jar build/Simulator.jar {} 1 1.0 spectrum_{}.txt charset.txt $outDirInt N

#float_simul1 --> floating-point weights
mkdir -p $outDirFloat
java -jar build/GenerateCharacterSet.jar 50 50 400 charset.txt $outDirFloat N N
seq 5 1 20 | xargs -I{} java -jar build/Simulator.jar {} 1 1.0 spectrum_{}.txt charset.txt $outDirFloat N
