# efficient_biocircuits
**gcc sundials_code.c -o sundials_code -lsundials_cvode -lsundials_nvecserial -lm**
**gcc -shared -o sundials_code_ctypes.so -fPIC sundials_code_ctypes.c -lsundials_cvode -lsundials_nvecserial -lm**
