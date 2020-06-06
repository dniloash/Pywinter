.. pywinter documentation master file, created by
   sphinx-quickstart on Sat Jun  6 13:42:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pywinter
====================================
Python WRF Intermediate Files

Pywinter is a Python3 library designed to create your own files in WPS intermediate file format (one step back to execute metgrid.exe). Usually you dont need to create the intermediate files by your own because that is the function of ungrib.exe, but sometimes you dont have your meteorological data in GRIB1 or GRIB2 format, It is the case of NetCDF data. Winter allows to create intermediate files by a simple way, the only prerequisite is read the data with python and next you can convert it on intermediate file format with pywinter.


Example
====================================

.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data


	## Read data
	############################################################

	# Read Geo-data
	lat = data.variables['Latitude'][:]
	lon = data.variables['Longitude'][:]
	# Read 2D data 
	tp2m = data.variables['T2'][:,:]
	# Read 3D data 
	temp = data.variables['T'][:,:,:]
	plevs = data.variables['PLEV'][:]
	# Read 3D soil data 
	soilt = data.variables['SLTY'][:,:,:]
	sl_layer = ['000010','010040','040100','100200']


	## Pywinter
	############################################################


	# Create winter fields
	winter_geo = pyw.Geo0(lat,lon)

	winter_t2m = pyw.V2d('TT',tp2m)
	winter_t = pyw.V3dp('TT',tp2m,plevs)
	winter_soilt_layer = pyw.Vsl('ST',soilt,sl_layer)

	# Listing fields
	total_fields = [winter_t2m,
	                winter_t,
	                winter_soilt_layer]

	# Out path
	path_out = '/home/Documents/intermediate_files/'

	# Write intermediate file
	pyw.cinter('FILE','1994-05-18_06',winter_geo,total_fields,path_out)



.. toctree::
   :maxdepth: 2
   :caption: Contents:

Installation
====================================
pip install pywinter

or

pip3 install pywinter

Docs
====================================

Documentation file :download:`pdf <https://github.com/dniloash/Pywinter/blob/master/docs/PYWINTER_documentation.pdf>`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
