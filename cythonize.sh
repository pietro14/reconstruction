echo "cythonizing the noise part"
cython cython_cygno.pyx
cythonize -a -i cython_cygno.pyx
