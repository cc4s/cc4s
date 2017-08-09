#!/bin/bash
# specify number of imaginary grid points
N=${1-7}
echo "N = $N"
# turn into 2 digit number
N0=$(awk "BEGIN {printf(\"%02d\", $N)}")

FIRST_HOLE=$(head -n 3 HoleEigenEnergies.dat | tail -n 1)
LAST_HOLE=$(tail -n 1 HoleEigenEnergies.dat)
FIRST_PARTICLE=$(head -n 3 ParticleEigenEnergies.dat | tail -n 1)
LAST_PARTICLE=$(tail -n 1 ParticleEigenEnergies.dat)
MU=$(awk "BEGIN {print ($FIRST_PARTICLEi+($LAST_HOLE))/2}")
MIN_EXCITATION=$(awk "BEGIN {print 2*($FIRST_PARTICLE-($LAST_HOLE))}")
MAX_EXCITATION=$(awk "BEGIN {print 2*($LAST_PARTICLE-($FIRST_HOLE))}")
R=$(awk "BEGIN {print $MAX_EXCITATION/($MIN_EXCITATION)}")
echo "min(eps_a-eps_i) = $MIN_EXCITATION"
echo "max(eps_a-eps_i) = $MAX_EXCITATION"
echo "R = max/min = $R"

# get base and exponent in scientif notation
R_SCIENTIFIC=$(awk "BEGIN {printf(\"%1.0e\n\", $R)}")
R_EXPONENT=$(echo $R_SCIENTIFIC | awk 'BEGIN {FS="e"} {print $2}')
# make sure FACTOR*10^EXPONENT is >= R
R=$(awk "BEGIN {print $R+0.5e$R_EXPONENT}")
R_SCIENTIFIC=$(awk "BEGIN {printf(\"%1.0e\n\", $R)}")
R_FACTOR=$(echo $R_SCIENTIFIC | awk 'BEGIN {FS="e"} {print $1+0}')
R_EXPONENT=$(echo $R_SCIENTIFIC | awk 'BEGIN {FS="e"} {print $2+0}')

# list suitable files
FILE=${0%/*}/reciprocalApproximation/1_xk${N0}.${R_FACTOR}_${R_EXPONENT}
if test ! -e $FILE
then
	echo "No samples available for the system, try larger N"
	exit 1
fi

# write taus
echo "ImaginaryTimePoints 1 $N" >ImaginaryTimePoints.dat
echo "i " >>ImaginaryTimePoints.dat
awk "/alpha/{print \$1/$MIN_EXCITATION}" $FILE >>ImaginaryTimePoints.dat

# write weights
echo "ImaginaryTimeWeights 1 $N" >ImaginaryTimeWeights.dat
echo "i " >>ImaginaryTimeWeights.dat
awk "/omega/{print \$1/$MIN_EXCITATION}" $FILE >>ImaginaryTimeWeights.dat

# write mu
echo "ChemicalPotential 0" >ChemicalPotential.dat
echo "  " >>ChemicalPotential.dat
echo "$MU" >>ChemicalPotential.dat
