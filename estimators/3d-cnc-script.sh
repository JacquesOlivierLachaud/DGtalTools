#!/bin/bash

# Path and filename of estimator program
EXEC=${HOME}/GITHUB/DGtalTools/build-release/estimators/3dFastCorrectedNormalCurrent

# H, G, HII, GII
QUANTITY=GII
H_FACTOR=0.75
H_START=1
H_END=0.02
# H_START=.05631351470947265625
# H_END=.05631351470947265625
# crixxi cylinder diabolo distel heart
# leopold rcube ellipsoid goursat sphere1 sphere9 torus
# or a polynomial
SHAPE=goursat

# Normal estimation
# True II VCM CTrivial Trivial
N_ESTIMATOR=II
N_ALPHA=0.33
N_BR=10
N_R=3
N_T=3

# Measure estimation
M_COEF=3
M_POW=0.5
# nothing or --crisp
CRISP=

# NOISE
NOISE=0

###############################################################################
case ${QUANTITY} in
    H|G|HII|GII) ;;
    *) echo "Invalid quantity: ${QUANTITY}"; exit 1;;
esac

# Computing filename
ERRFILE="${SHAPE}-${QUANTITY}-${N_ESTIMATOR}-a${N_ALPHA}-br${N_BR}-r${N_R}-mc${M_COEF}-mp${M_POW}${CRISP}-N${NOISE}.txt"
H=${H_START}
while test `echo "$H >= ${H_END}" | bc -l` -eq 1; do
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "------------------ ${SHAPE} at ${H}, ${QUANTITY} -------------------- "
    ${EXEC} -p ${SHAPE} -Q ${QUANTITY} -e ${N_ESTIMATOR} -R ${N_BR} -r ${N_R} --alpha ${N_ALPHA} -t ${N_T} -g ${H} --m-coef ${M_COEF} --m-pow ${M_POW} -N ${NOISE} ${CRISP} -V None --error ${ERRFILE}
    H=`echo "$H * ${H_FACTOR}" | bc -l`
done
    
