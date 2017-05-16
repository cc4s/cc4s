#!/bin/bash

H0_I=$(paste HoleEigenEnergies.dat DrccdHoleOccupancies.dat  | tail -n +3 | awk "{h+=\$1*(\$2-2)} END {printf(\"%21.17e\\n\", h)}")
H0=$(paste ParticleEigenEnergies.dat DrccdParticleOccupancies.dat  | tail -n +3 | awk "{h+=\$1*\$2} END {printf(\"%21.17e\\n\", h+$H0_I)}")
H1=$(tail DrccdCoulombExpectationValue.dat -n +3)
H=$(tail DrccdEnergy.dat -n +3)
D=$(awk "BEGIN {print ($H0)+($H1) - ($H)}")
echo "<H0>_cc + <H1>_cc = E_cc + error."
echo "$H0 + $H1 = $H + $D"
