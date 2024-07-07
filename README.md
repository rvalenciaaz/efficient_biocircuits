# efficient_biocircuits

This repository contains reimplementations in C/SUNDIALS of gene circuits ODE models found in the biocircuits course by the Elowitz lab at Caltech. The motivation behind this project is to first introduce the typical format of a SUNDIALS ODE solver and how to integrate the different solver components in a C script, and second to have a set of very efficient examples for biological ODE solving, that can be useful on high-throughput testing of big circuit architectures.

`gcc sundials_code.c -o sundials_code -lsundials_cvode -lsundials_nvecserial -lm`

`gcc -shared -o sundials_code_ctypes.so -fPIC sundials_code_ctypes.c -lsundials_cvode -lsundials_nvecserial -lm`
