#!/bin/bash
rm *LOG;
source /home/user/Sim/geant4.10.04.p02-install/bin/geant4.sh;
make;
./exec_ONE run.mac >> inLOG;
