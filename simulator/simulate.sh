#!/bin/bash

#simul1 
mkdirhier data/int_simul1
java -jar build/GenerateCharacterSet.jar 50 50 400 charset.txt data/int_simul1 Y
seq 5 1 20 | xargs -I{} java -jar build/Simulator.jar {} 10 1.0 spectrum_{}.txt charset.txt data/int_simul1

#simul2
mkdirhier data/float_simul1
java -jar build/GenerateCharacterSet.jar 50 50 400 charset.txt data/float_simul1 Y
seq 5 1 20 | xargs -I{} java -jar build/Simulator.jar {} 10 1.0 spectrum_{}.txt charset.txt data/float_simul1
