#!/bin/bash


GVECfile="../GVECDUMP"
tail -n +3 "${GVECfile}" > "${GVECfile}".tmp
lines=`cat ${GVECfile}.tmp | wc -l`
#lines=$(($lines-2))
echo "Momenta 2 3 $lines" > Momenta.dat
echo "ij " >> Momenta.dat
awk -v OFS='\n' '{print $1, $2, $3}' ${GVECfile}.tmp >> Momenta.dat

###make CoulombKernel.dat file############
echo "CoulombKernel 1 $lines" > CoulombKernel.dat
echo "i " >> CoulombKernel.dat
awk '{print $4}' ${GVECfile}.tmp >> CoulombKernel.dat

###parameters needed in cc4s file############
K1=$(grep -A1 Gamma ../KPOINTS | tail -1 | awk '{print $1}')
K2=$(grep -A1 Gamma ../KPOINTS | tail -1 | awk '{print $2}')
K3=$(grep -A1 Gamma ../KPOINTS | tail -1 | awk '{print $3}')
kpoints=$(($K1*$K2*K3))
echo "kpoints=$kpoints"
volume=`grep "volume" ../OUTCAR | tail -1 | awk '{print $5}'`
echo "volume=$volume"
constantFactor=`awk '{print sprintf("%.11f", ($1^2+$2^2+$3^2)*$4)}' $GVECfile.tmp | tail -n1`
echo "constantFactor=$constantFactor"
