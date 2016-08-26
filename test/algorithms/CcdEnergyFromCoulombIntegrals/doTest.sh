# Here the class of this test case is defined
# @CLASS=essential,ccd


# RUN_COMMAND and CC4S_PATH path are globally
# defined in the main test.sh script
${RUN_COMMAND} ${CC4S_PATH} -file ccd.cc4s


ENERGY=$(readScalar CcdEnergy.dat)
COMPARE_ENERGY=$COULOMBVERTEX_CCD
TOLERANCE=1e-11

#echo ${ENERGY} >&2
#echo ${COMPARE_ENERGY} >&2

# If the test succeeds then TEST_RESULT=0
# If the test fails then    TEST_RESULT=1

TEST_RESULT=$(
python -c "print(0 if (abs(${ENERGY} - ${COMPARE_ENERGY})<${TOLERANCE}) else 1)"
)




