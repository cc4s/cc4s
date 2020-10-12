OBJECTS= \
build/${CONFIG}/obj/main/Options.o \
build/${CONFIG}/obj/main/Setting.o \
build/${CONFIG}/obj/main/util/Log.o \
build/${CONFIG}/obj/main/util/Timer.o \
build/${CONFIG}/obj/main/util/TensorIo.o \
build/${CONFIG}/obj/main/util/BlacsWorld.o \
build/${CONFIG}/obj/main/tcc/Tcc.o \
build/${CONFIG}/obj/main/engines/DryTensor.o \
build/${CONFIG}/obj/main/algorithms/Algorithm.o \
build/${CONFIG}/obj/main/algorithms/TensorReader.o \
build/${CONFIG}/obj/main/algorithms/TensorWriter.o \
build/${CONFIG}/obj/main/algorithms/TensorNetwork.o \
build/${CONFIG}/obj/main/algorithms/DefineHolesAndParticles.o \
build/${CONFIG}/obj/main/algorithms/SliceOperator.o \
build/${CONFIG}/obj/main/algorithms/CoulombIntegralsFromVertex.o \
build/${CONFIG}/obj/main/algorithms/Mp2EnergyFromCoulombIntegrals.o \
build/${CONFIG}/obj/main/algorithms/ClusterSinglesDoublesAlgorithm.o \
build/${CONFIG}/obj/main/algorithms/DrccdEnergyFromCoulombIntegrals.o \
build/${CONFIG}/obj/main/algorithms/CcsdEnergyFromCoulombIntegralsReference.o \
build/${CONFIG}/obj/main/algorithms/CcsdEnergyFromCoulombIntegrals.o \
build/${CONFIG}/obj/main/mixers/Mixer.o \
build/${CONFIG}/obj/main/mixers/LinearMixer.o \
build/${CONFIG}/obj/main/mixers/DiisMixer.o \

TESTS_OBJECTS = \
build/${CONFIG}/obj/test/math/Antisymmetrize.o \
build/${CONFIG}/obj/test/math/Symmetrize.o \
build/${CONFIG}/obj/test/math/PermutationSign.o \
build/${CONFIG}/obj/test/math/FockVector.o \
build/${CONFIG}/obj/test/util/LapackGeneralEigenSystem.o \
build/${CONFIG}/obj/test/Test.o  \
build/${CONFIG}/obj/test/algorithms/UccsdAmplitudesFromCoulombIntegrals.o  \

# build/${CONFIG}/obj/test/math/EigenSystemDavidson.o \