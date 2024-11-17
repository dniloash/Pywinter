# Pywinter
Python lib for Reading/Creating MPAS/WRF-WPS intermediate file

Pywinter is a Python3 library designed for handling files in MPAS/WRF-WPS intermediate file format. Usually you don't need to deal with the intermediate files by your own because that is the function of ungrib.exe, but sometimes you don't have your meteorological data in GRIB format. Pywinter allows to read and create intermediate files by a simple way.

Please check pywinter documentation in /docs and /examples.

# Installation

# Dependencies
pywinter requires:

Fortran compiler (gfortran)

Python (>= 3.9 )

Numpy

Pywinter has been tested in linux distributions.

# User installation

pip install pywinter

or

pip3 install pywinter
