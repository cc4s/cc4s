SRC_FILES += \
main/Cc4s.cxx \
main/Options.cxx \
main/Writer.cxx \
main/Reader.cxx \
main/Setting.cxx \
main/util/Log.cxx \
main/util/Timer.cxx \
main/util/TensorIo.cxx \
main/util/BlacsWorld.cxx \
main/tcc/Tcc.cxx \
main/engines/DryTensor.cxx \
main/algorithms/Algorithm.cxx \
main/algorithms/Write.cxx \
main/algorithms/Read.cxx \
main/algorithms/TensorReader.cxx \
main/algorithms/TensorWriter.cxx \
main/algorithms/DefineHolesAndParticles.cxx \
main/algorithms/SliceOperator.cxx \
main/algorithms/CoulombIntegralsFromVertex.cxx \
main/algorithms/Mp2EnergyFromCoulombIntegrals.cxx \
main/algorithms/ClusterSinglesDoublesAlgorithm.cxx \
main/algorithms/DrccdEnergyFromCoulombIntegrals.cxx \
main/algorithms/CcsdEnergyFromCoulombIntegralsReference.cxx \
main/algorithms/CcsdEnergyFromCoulombIntegrals.cxx \
main/algorithms/PerturbativeTriplesReference.cxx \
main/algorithms/StructureFactor.cxx \
main/algorithms/FiniteSizeCorrection.cxx \
main/algorithms/BasisSetCorrection.cxx \
main/mixers/Mixer.cxx \
main/mixers/LinearMixer.cxx \
main/mixers/DiisMixer.cxx \
main/math/TensorUnion.cxx \

ifeq ($(ATRIP), yes)
SRC_FILES += main/algorithms/Atrip.cxx
endif

TEST_SRC_FILES = \
test/math/Antisymmetrize.cxx \
test/math/Symmetrize.cxx \
test/math/PermutaticxxnSign.cxx \
test/math/TensorUnion.cxx \
test/util/LapackGeneralEigenSystem.cxx \
test/Test.cxx  \
test/algorithms/UccsdAmplitudesFromCoulombIntegrals.cxx  \
