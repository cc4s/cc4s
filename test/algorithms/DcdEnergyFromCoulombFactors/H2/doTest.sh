 # Here the class of this test case is defined
# @CLASS=H2,ccd,dcd,trd
TEST_DESCRIPTION="Testing system: H2 \n      10x10x10 Angstrom supercell \n      Ecut=300eV, EcutGW=150eV \n      AVTZ: NG=2070 No=1 Nv=46 \n      G-vector reduction: Ng=0.5*NG=1035 \n      TRD: NR=2*Ng=2070 \n      DCD correlation energy from reducted G-vector and decomposed Coulomb vertex: check if the energy is within tolerance"


# RUN_COMMAND and CC4S_PATH path are globally
# defined in the main test.sh script
${RUN_COMMAND} ${CC4S_PATH} -file dcd.cc4s


ENERGY=$(readScalar DcdEnergy.dat) 
ENERGYSLICE=$(readScalar DcdEnergySlice.dat) 
ENERGYREDUCE=$(readScalar DcdEnergyReduce.dat)
ENERGYSLICEREDUCE=$(readScalar DcdEnergySliceReduce.dat)
ENERGYREDUCETRD=$(readScalar DcdEnergyReduceTRD.dat)
ENERGYSLICEREDUCETRD=$(readScalar DcdEnergySliceReduceTRD.dat)
ENERGYFACTORSREDUCETRD=$(readScalar DcdEnergyFactorsReduceTRD.dat)
TOLERANCE=1e-11

# If the test succeeds then TEST_RESULT=0
# If the test fails then    TEST_RESULT=1

TEST_RESULT=$(
python <<EOF 
if abs(${ENERGY} - ${DCDENERGY})<${TOLERANCE} and abs(${ENERGYSLICE} - ${DCDENERGY})<${TOLERANCE} and abs(${ENERGYREDUCE} - ${DCDENERGYREDUCE})<${TOLERANCE} and abs(${ENERGYSLICEREDUCE} - ${DCDENERGYREDUCE})<${TOLERANCE} and abs(${ENERGYREDUCETRD} - ${DCDENERGYREDUCETRD})<${TOLERANCE} and abs(${ENERGYSLICEREDUCETRD} - ${DCDENERGYREDUCETRD})<${TOLERANCE} and abs(${ENERGYFACTORSREDUCETRD} - ${DCDENERGYREDUCETRD})<${TOLERANCE}:
    print(0)
else:
    print(1)
EOF
)

echoDebug Energy from Integrals: $ENERGY
echoDebug Energy from sliced Integrals: $ENERGYSLICE
echoDebug Energy from Integrals: Reduced G-vector: $ENERGYREDUCE
echoDebug Energy from sliced Integrals: Reduced G-vector: $ENERGYSLICEREDUCE
echoDebug Energy from Integrals: Reduced G-vector - TRD: $ENERGYREDUCETRD
echoDebug Energy from sliced Integrals: Reduced G-vector - TRD: $ENERGYSLICEREDUCETRD
echoDebug Energy from Coulomb Factors: Reduced G-vector - TRD: $ENERGYFACTORSREDUCETRD
