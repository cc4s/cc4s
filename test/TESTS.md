
# hf molecule dz
- path: ./test/tests/hf
- tags:
	- essential
	- nwchem
	- molecule
	- gaussian
	- ccsd
	- ccsd(t)
	- cT
- Run and check command: `./test/run 4 --run --check --rx "hf molecule dz"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- PerturbativeCT
	- PerturbativeTriples
	- PerturbativeTriplesReference

# h2o2 molecule mrcc
- path: ./test/tests/h2o2
- tags:
	- essential
	- molecule
	- mrcc
	- ccsd
	- ccsd(t)
	- fno
- Run and check command: `./test/run 4 --run --check --rx "h2o2 molecule mrcc"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- MrccReader
	- PerturbativeCT
	- PerturbativeTriples

# h2o molecule nwchem cc-pvdz
- path: ./test/tests/h2o_cc-pvdz
- tags:
	- essential
	- molecule
	- nwchem
	- gto
	- ccsd
	- lambda-ccsd
	- (t-2)
	- (t-3)
	- (t-4)
	- (3cc)_nonit
- Run and check command: `./test/run 4 --run --check --rx "h2o molecule nwchem cc-pvdz"`
- Algorithms tested:
	- CcsdtNonit 
	- CcsdtNonit 3CC
	- CoupledCluster Ccsd 
	- CoupledCluster LambdaCcsd 

# h2o molecule aug-cc-pvdz
- path: ./test/tests/h2o_dz
- tags:
	- essential
	- nwchem
	- molecule
	- gaussian
	- ccsd
	- ccsd(t)
- Run and check command: `./test/run 4 --run --check --rx "h2o molecule aug-cc-pvdz"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- PerturbativeTriples

# h2o mp2 fno aug-cc-pvqz
- path: ./test/tests/h2o_fno-qz
- tags:
	- essential
	- nwchem
	- molecule
	- gaussian
	- ccsd
	- ccsd(t)
	- basis-set-correction
- Run and check command: `./test/run 4 --run --check --rx "h2o mp2 fno aug-cc-pvqz"`
- Algorithms tested:
	- BasisSetCorrection
	- CoupledCluster Ccsd 
	- PerturbativeTriples

# n2 molecule nwchem cc-pvdz chan-kallay-gauss
- path: ./test/tests/n2_cc-pvdz-chan
- tags:
	- essential
	- molecule
	- nwchem
	- gto
	- ccsd
	- dcsd
	- ccsdt
	- dcsdt
	- 3cc
- Run and check command: `./test/run 4 --run --check --rx "n2 molecule nwchem cc-pvdz chan-kallay-gauss"`
- Algorithms tested:
	- CoupledCluster 3CC 
	- CoupledCluster Ccsd 
	- CoupledCluster Ccsdt 
	- CoupledCluster Ccsdt-1b 
	- CoupledCluster Ccsdt-3 
	- CoupledCluster Dcsd 
	- CoupledCluster Dcsdt 

# n2 molecule nwchem cc-pvdz HEAT benchmark
- path: ./test/tests/n2_cc-pvdz-heat
- tags:
	- essential
	- molecule
	- nwchem
	- gto
	- ccsd
	- ccsdt
	- ccsdt-1b
	- ccsdt-3
	- 3cc
	- dcsdt
- Run and check command: `./test/run 4 --run --check --rx "n2 molecule nwchem cc-pvdz HEAT benchmark"`
- Algorithms tested:
	- CoupledCluster 3CC 
	- CoupledCluster Ccsd 
	- CoupledCluster Ccsdt 
	- CoupledCluster Ccsdt-1b 
	- CoupledCluster Ccsdt-3 
	- CoupledCluster Dcsdt 

# n2 molecule vasp
- path: ./test/tests/n2_vasp
- tags:
	- essential
	- molecule
	- vasp
	- ccsd
	- ccsd(t)
	- cT
- Run and check command: `./test/run 4 --run --check --rx "n2 molecule vasp"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- PerturbativeCT
	- PerturbativeTriples

# r-BN Γ-point
- path: ./test/tests/r-bn_gamma
- tags:
	- essential
	- vasp
	- basis-set-correction
	- finite-size-correction
	- real
- Run and check command: `./test/run 4 --run --check --rx "r-BN Γ-point"`
- Algorithms tested:
	- BasisSetCorrection
	- CoupledCluster Ccsd 
	- FiniteSizeCorrection
	- PerturbativeTriples

# r-BN shifted single k-point
- path: ./test/tests/r-bn_shifted
- tags:
	- essential
	- vasp
	- finite-size-correction
	- basis-set-correction
	- complex
- Run and check command: `./test/run 4 --run --check --rx "r-BN shifted single k-point"`
- Algorithms tested:
	- BasisSetCorrection
	- CoupledCluster Ccsd 
	- FiniteSizeCorrection
	- PerturbativeTriples

# ueg model rs1_no7_nv26 complex
- path: ./test/tests/ueg_complex-7occ-26virt
- tags:
	- essential
	- ueg
	- drccd
	- ccsd
	- ccsd(t)
	- (cT)
	- complex
- Run and check command: `./test/run 4 --run --check --rx "ueg model rs1_no7_nv26 complex"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- CoupledCluster Drccd 
	- PerturbativeTriples
	- PerturbativeTriplesReference

# ueg model rs1_no7_nv26 real
- path: ./test/tests/ueg_real-rs1.0-7occ-26virt
- tags:
	- essential
	- ueg
	- drccd
	- ccsd
	- ccsd(t)
	- (cT)
	- real
- Run and check command: `./test/run 4 --run --check --rx "ueg model rs1_no7_nv26 real"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- CoupledCluster Drccd 
	- PerturbativeTriples
	- PerturbativeTriplesReference

# ueg model rs1_no123_nv1114
- path: ./test/tests/ueg_rs1.0-123occ-1114virt
- tags:
	- ueg
	- ccsd(t)-performance
- Run and check command: `./test/run 4 --run --check --rx "ueg model rs1_no123_nv1114"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- UegVertexGenerator

# ueg model rs1_no19-nv232
- path: ./test/tests/ueg_rs1.0-19occ-232virt
- tags:
	- ueg
	- ct-atrip
- Run and check command: `./test/run 4 --run --check --rx "ueg model rs1_no19-nv232"`
- Algorithms tested:
	- CoupledCluster Ccsd 
	- PerturbativeCT
	- TensorSetToTensor
	- UegVertexGenerator

# h2o molecule vasp
- path: ./test/tests/h2o_vasp
- tags:
	- essential
	- molecule
	- ccsd
	- ccsd(t)
- Run and check command: `./test/run 4 --run --check --rx "h2o molecule vasp"`
- Algorithms tested:
	- CcsdtNonit 
	- CcsdtNonit 3CC
	- CoupledCluster 3CC 
	- CoupledCluster ApxRpa 
	- CoupledCluster Ccsd 
	- CoupledCluster Ccsdt-1b 
	- CoupledCluster Ccsdt-3 
	- CoupledCluster Dcsd 
	- CoupledCluster Dcsdt 
	- CoupledCluster LambdaCcsd LambdaCcsd
	- CoupledCluster LambdaCcsd LambdaCcsdReference
	- CoupledCluster Ldrccd 
	- PerturbativeCT
	- PerturbativeTriples
	- PerturbativeTriplesReference
