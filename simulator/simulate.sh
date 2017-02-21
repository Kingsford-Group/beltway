#!/bin/bash

#int_simul1 --> integer weights 
mkdirhier data/int_simul2
java -jar build/GenerateCharacterSet.jar 50 50 400 charset.txt data/int_simul2 Y
seq 5 1 20 | xargs -I{} java -jar build/Simulator.jar {} 1 1.0 spectrum_{}.txt charset.txt data/int_simul2

#float_simul1 --> floating-point weights
mkdirhier data/float_simul2
java -jar build/GenerateCharacterSet.jar 50 50 400 charset.txt data/float_simul2 Y
seq 5 1 20 | xargs -I{} java -jar build/Simulator.jar {} 1 1.0 spectrum_{}.txt charset.txt data/float_simul2
