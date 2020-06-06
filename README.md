# Pywinter
WRF-WPS intermediate file python writer

Pywinter is a Python3 library designed to create your own files in WPS intermediate file format (one step back to execute metgrid.exe).  Usually you dont need to create the intermediate files by your ownbecause that is the function of ungrib.exe, but sometimes you dont have your meteorological datain GRIB1 or GRIB2 format.

Please check pywinter documentation in /docs and /examples.

# Installation

# Dependencies
pywinter requires:

Fortran compiler (gfortran)

Python (>= 3.6)

Numpy 

Pywinter has been tested in linux distributions.

# User installation

pip install pywinter

or

pip3 install pywinter
