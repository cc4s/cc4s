SRC_FILES += \
main/Cc4s.cxx \
main/Options.cxx \
main/Writer.cxx \
main/Reader.cxx \
main/Setting.cxx \
main/Log.cxx \
main/Timer.cxx \
main/TensorIo.cxx \
main/TensorSet.cxx \
main/tcc/Tcc.cxx \
main/engines/DryTensor.cxx \
main/mixers/Mixer.cxx \
main/mixers/LinearMixer.cxx \
main/mixers/DiisMixer.cxx \
main/algorithms/Algorithm.cxx \
main/algorithms/Write.cxx \
main/algorithms/Read.cxx \
main/algorithms/TensorReader.cxx \
main/algorithms/TensorWriter.cxx \
main/algorithms/DefineHolesAndParticles.cxx \
main/algorithms/SliceOperator.cxx \
main/algorithms/DimensionProperty.cxx \
main/algorithms/NonZeroCondition.cxx \
main/algorithms/VertexCoulombIntegrals.cxx \
main/algorithms/SecondOrderPerturbationTheory.cxx \
main/algorithms/CoupledCluster.cxx \
main/algorithms/coupledcluster/CoupledClusterMethod.cxx \
main/algorithms/coupledcluster/Ccsd.cxx \
main/algorithms/coupledcluster/CcsdReference.cxx \
main/algorithms/coupledcluster/Drccd.cxx \
main/algorithms/PerturbativeTriplesReference.cxx \
main/algorithms/FiniteSizeCorrection.cxx \
main/algorithms/BasisSetCorrection.cxx \
main/algorithms/UegVertexGenerator.cxx \

ifeq ($(ATRIP), yes)
SRC_FILES += main/algorithms/PerturbativeTriples.cxx
endif

TEST_SRC_FILES = \
test/math/Antisymmetrize.cxx \
test/math/Symmetrize.cxx \
test/math/PermutaticxxnSign.cxx \
test/math/TensorUnion.cxx \
test/util/LapackGeneralEigenSystem.cxx \
test/Test.cxx  \
test/algorithms/UccsdAmplitudesFromCoulombIntegrals.cxx  \
