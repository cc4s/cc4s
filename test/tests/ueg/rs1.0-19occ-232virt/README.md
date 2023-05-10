This is a hands-on for running (cT) calculations with atrip. \
Make sure that your atrip version is pretty new:
https://github.com/alejandrogallo/atrip commit: 282ef27 


It is assumend that the calculation is that large that a reference-(T) calculation is no more feasible.

1.) Run a normal ccsd calculation \
2.) run the algorithm PerturbativeTriplesReference using the flag:   ```intermediateOnly: 1``` and put ```cTIntermediate: CtIntermediate``` as the algoirthm's output.\
3.) write: CtIntermediate, EigenEnergies, Amplitudes, CoulombIntegrals to disk. ( binary: 1 !!!!)


!!!IMPORTANT!!!\
It is necessary to write the EigenEnergies in the following way (otherwise you will enter the land of pain)
```


- name: Write
  in:
    fileName: "EigenEnergiesSliced.yaml"
    source: EigenEnergies
    binary: 1
  out: {}
```

The run will most probably crash because he tries to write the PPPP integral to disk.\
If not one should remove CoulombIntegrals which are not needed.
You just need hhhp pphh ppph for atrip.

rm CoulombIntegrals.components.hhhh.elements CoulombIntegrals.components.hhpp.elements \
   CoulombIntegrals.components.hhph.elements CoulombIntegrals.components.hphp.elements \
   CoulombIntegrals.components.hppp.elements CoulombIntegrals.components.phpp.elements \
   CoulombIntegrals.components.pphp.elements CoulombIntegrals.components.phhp.elements \
   CoulombIntegrals.components.phph.elements CoulombIntegrals.components.phhh.elements \
   CoulombIntegrals.components.hpph.elements CoulombIntegrals.components.hphh.elements

 


follow the instructions on
https://github.com/airmler/atrip_resources


```
irmler@andoria01:~/Programs/cc4s-develop/test/tests/ueg/rs1.0-27occ-560virt$ mpirun -np 64 ~/Programs/cc4s-develop/build/gcc-oblas-ompi/bin/Cc4s

                __ __
     __________/ // / _____
    / ___/ ___/ // /_/ ___/
   / /__/ /__/__  __(__  )
   \___/\___/  /_/ /____/
  Coupled Cluster for Solids

version: heads/develop-0-g983439bb-dirty, date: Wed Mar 1 14:03:13 2023 +0100
build date: May  4 2023 12:06:44
compiler: g++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
mpi: mpirun (Open MPI) 4.0.3
total processes: 64
calculation started on: Sat May  6 22:02:56 2023

DRY RUN ONLY - nothing will be calculated
Dry run finished. Estimates provided for 64 ranks.
Memory estimate (per Rank/Total): 1.07601 / 68.8645 GB
Operations estimate (per Rank/Total): 12787.2 / 818383 GFLOPS
Time estimate with assumed performance of 10 GFLOPS/core/s: 1278.72 s (0.355201 h)
--
step: 1, UegVertexGenerator
System Information:
  rs 1, No 19, Nv232
  Volume 159.1740278, madelung 0.5235371088
  HOMO 0.4177839456, LUMO 1.612902262
  Reference Energy per Electron/total 0.5666218456/21.53163013
realtime 1.033732682 s
--
step: 2, DefineHolesAndParticles
number of holes     No: 19
number of particles Nv: 232
number of states    Np: 251
realtime 0.000104767 s
--
step: 3, SliceOperator
Slicing CoulombVertex into holes and particles.
realtime 0.350618145 s
--
step: 4, VertexCoulombIntegrals
Using complex Coulomb integrals
number of field variables NF: 1665
realtime 0.626437181 s
--
step: 5, CoupledCluster
Using method Ccsd. integralsSliceSize: 56
Using mixer DiisMixer. maxResidua: 4
Maximum number of iterations: 20
Unless reaching energy convergence dE: 1e-08
and amplitudes convergence dR: 1e-08
Iter         Energy         dE           dR         time   GF/s/core
   1  -1.52132677e+00  -1.5213e+00   4.6696e-01      1.3   16.1
   2  -1.42634181e+00   9.4985e-02   5.6879e-02     48.4   11.8
   3  -1.45015637e+00  -2.3815e-02   8.1449e-03     47.4   11.7
   4  -1.45065845e+00  -5.0207e-04   9.6968e-04     47.2   11.8
   5  -1.45075071e+00  -9.2265e-05   1.5774e-04     48.2   11.6
   6  -1.45075259e+00  -1.8759e-06   2.2932e-05     47.8   11.7
   7  -1.45075253e+00   5.4670e-08   3.6412e-06     48.1   11.6
   8  -1.45075292e+00  -3.9087e-07   7.1757e-07     48.2   11.6
   9  -1.45075295e+00  -2.4598e-08   9.7528e-08     47.9   11.6
  10  -1.45075296e+00  -1.3979e-08   2.7457e-08     48.4   11.5
  11  -1.45075296e+00   1.7020e-09   4.1276e-09     49.7   11.2

Ccsd correlation energy:          -1.4507529611
2nd-order correlation energy:     -1.5213267733
realtime 483.033570705 s
--
step: 6, PerturbativeTriplesReference
realtime 86.253544494 s
--
step: 7, PerturbativeTriples
Progress(%)  time(s)   GFLOP/s
1            2         6.857
10           18        6.813
20           20        6.795
30           20        6.800
40           20        6.805
50           20        6.798
60           20        6.794
70           20        6.793
80           20        6.792
90           20        6.791
:q100          20        6.792
(T) correlation energy:      -0.064410258562274
```



```
irmler@cqc05:~/Projects/renormTriples/nikos_ueg$ ./twistUeg -e 38 -T -g 11
Number of threads: 1
- Estimate for the memory required: 0.00811201 GB

- General information of the system
  - Lattice type                  : sc
  - Lattice constant in Bohr units: 5.41948
  - Box-volume in Bohr units      : 159.174
  - WS radius in Bohr units       : 1
  - Reference interacting/free    : interacting
  - Madelung constant             : -0.523537
  - Occupied orbitals             : 19
  - Virtuals orbitals             : 232
  - Shift                         : 0 0 0


[Energy] Hartree Fock Energy: 0.566622
[Energy] HOMO/LUMO/BandGap: 0.417784 / 1.6129 / 1.19512

[Energy]             MP2: ( -0.05472909 | +0.01469417 | -0.04003491497677 | 232 )
It.      Energy       Difference    Time
  1 -0.037535310698 -0.037535310698    0.1 s
  2 -0.038162009812 -0.000626699114    0.1 s
  3 -0.038175222299 -0.000013212487    0.1 s
  4 -0.038177740714 -0.000002518415    0.1 s
  5 -0.038177690203 +0.000000050511    0.1 s
  6 -0.038177699543 -0.000000009340    0.1 s
[Energy]             CCD: ( -0.05433114 | +0.01615344 | -0.03817769954300 | 232 )
++++++++++
[Energy]       CRTriples: (      -      |      -      | -0.00152173790975 | 232 )
[Energy]      CRTriplesD: (      -      |      -      | -0.00138641357742 | 232 )
[Energy]         Triples: (      -      |      -      | -0.00169500541235 | 232 )
```
