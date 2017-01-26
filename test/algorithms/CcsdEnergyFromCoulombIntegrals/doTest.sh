# Here the class of this test case is defined
# @CLASS=essential,ccd,trd
TEST_DESCRIPTION="CCSD energy: check if the energy is within tolerance"


# RUN_COMMAND and CC4S_PATH path are globally
# defined in the main test.sh script
${RUN_COMMAND} ${CC4S_PATH} -file ccsd.cc4s


ENERGY=$(readScalar CcsdEnergy.dat)
ENERGY_REEVALUATED=$(readScalar CcsdEnergyReevaluated.dat)
ENERGYSLICE=$(readScalar CcsdEnergySlice.dat)
ENERGYFACTORS=$(readScalar CcsdEnergyFactors.dat)
COMPARE_ENERGY_VERTEX=$COULOMBVERTEX_CCSD
COMPARE_ENERGY_FACTORS=$COULOMBFACTORS_CCSD
TOLERANCE=1e-11

# If the test succeeds then TEST_RESULT=0
# If the test fails then    TEST_RESULT=1

TEST_RESULT=$(
python <<EOF 
if abs(${ENERGY} - ${COMPARE_ENERGY_VERTEX})<${TOLERANCE} and abs(${ENERGYSLICE} - ${COMPARE_ENERGY_VERTEX})<${TOLERANCE} and abs(${ENERGYFACTORS} - ${COMPARE_ENERGY_FACTORS})<${TOLERANCE}:
    print(0)
else:
    print(1)
EOF
)

echoDebug Energy from Integrals: $ENERGY
echoDebug Energy reevaluated from given amplitudes: $ENERGY_REEVALUATED
echoDebug Energy from sliced Integrals: $ENERGYSLICE
echoDebug Energy from Coulomb Factors: $ENERGYFACTORS
