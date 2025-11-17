# Pywinter
Python lib for Reading/Creating MPAS/WRF-WPS intermediate file

Pywinter is a Python3 library designed for handling files in MPAS/WRF-WPS intermediate file format. Usually, users don't need to manipulate intermediate files directly because that is the role of ungrib.exe, but sometimes the meteorological data is not avalaible in GRIB format. Pywinter offers a simple way to read and create intermediate files.

Please check pywinter documentation in /docs and /examples.

# Dependencies
pywinter requires:

Fortran compiler (gfortran)

Python (>= 3.9 )

Numpy

Pywinter has been tested in linux distributions.

# User installation

pip install pywinter

or

pip install git+https://github.com/dniloash/Pywinter
