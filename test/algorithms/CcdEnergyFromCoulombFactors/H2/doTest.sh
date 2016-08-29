 # Here the class of this test case is defined
# @CLASS=H2,ccd,trd
TEST_DESCRIPTION=" Testing system: H2"

cat <<LONG_TEST_DESCRIPTION
      10x10x10 Angstrom supercell 
      Ecut=300eV, EcutGW=150eV 
      AVTZ: NG=2070 No=1 Nv=46 
      G-vector reduction: Ng=0.5*NG=1035 
      TRD: NR=2*Ng=2070 

      CCD correlation energy from reducted G-vector and decomposed Coulomb
      vertex: check if the energy is within tolerance
LONG_TEST_DESCRIPTION

# RUN_COMMAND and CC4S_PATH path are globally
# defined in the main test.sh script
${RUN_COMMAND} ${CC4S_PATH} -file ccd.cc4s


ENERGY=$(readScalar CcdEnergy.dat) 
ENERGYSLICE=$(readScalar CcdEnergySlice.dat) 
ENERGYREDUCE=$(readScalar CcdEnergyReduce.dat)
ENERGYSLICEREDUCE=$(readScalar CcdEnergySliceReduce.dat)
ENERGYREDUCETRD=$(readScalar CcdEnergyReduceTRD.dat)
ENERGYSLICEREDUCETRD=$(readScalar CcdEnergySliceReduceTRD.dat)
ENERGYFACTORSREDUCETRD=$(readScalar CcdEnergyFactorsReduceTRD.dat)
TOLERANCE=1e-11

# If the test succeeds then TEST_RESULT=0
# If the test fails then    TEST_RESULT=1

TEST_RESULT=$(
python <<EOF 
if abs(${ENERGY} - ${CCDENERGY})<${TOLERANCE} and abs(${ENERGYSLICE} - ${CCDENERGY})<${TOLERANCE} and abs(${ENERGYREDUCE} - ${CCDENERGYREDUCE})<${TOLERANCE} and abs(${ENERGYSLICEREDUCE} - ${CCDENERGYREDUCE})<${TOLERANCE} and abs(${ENERGYREDUCETRD} - ${CCDENERGYREDUCETRD})<${TOLERANCE} and abs(${ENERGYSLICEREDUCETRD} - ${CCDENERGYREDUCETRD})<${TOLERANCE} and abs(${ENERGYFACTORSREDUCETRD} - ${CCDENERGYREDUCETRD})<${TOLERANCE}:
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
