.. pywinter documentation master file, created by
   sphinx-quickstart on Sat Jun  6 13:42:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pywinter
====================================
Python WRF-WPS Intermediate Files

Pywinter is a Python3 library designed for handling files in WRF-WPS intermediate file format. Usually you don't need to deal with the intermediate files by your own because that is the function of ungrib.exe, but sometimes you don't have your meteorological data in GRIB format. Pywinter allows to read and create intermediate files by a simple way.



.. toctree::
   :maxdepth: 2
   :caption: Contents:


Dependencies
====================================
pywinter requires:

Pywinter has been tested in linux distributions.

Fortran compiler (gfortran)

Python (>= 3.6)

Numpy

Installation
====================================
pip install pywinter

or

pip3 install pywinter


Read intermediate files (rinter)
====================================
For reading the intermediate files information you must utilize the function \textbf{rinter}. Once you have read the information you can manipulate the data easily. The results will be a dictionary that contains the variables in the intermediate file, also every key includes general information, geo information, levels, and the data array.

Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw

	infile = '/home/allyson/Documents/files/FILE:1994-05-18_06'

	interfile = pyw.rinter(infile)

	print(interfile.keys())
	>> dict_keys(['LANDSEA', 'ST', 'SST', 'SOILHGT', 'PSFC', 'HGT', 'SKINTEMP', 'TT',
	'PMSL', 'VV', 'SM', 'UU', 'RH', 'TT2M', 'RH2M', 'UU10M', 'VV10M'])

	print(interfile['TT'].general)
	>> {'VERSION': 5, 'HDATE': '2015-07-27_12:00:00', 'XFCST': 0.0, 'MAP_SOURCE': 'ECMWF',
	'FIELD': 'TT', 'UNITS': 'K', 'DESC': 'Temperature', 'XLVL': '1000', 'NX': 441, 'NY': 329,
	'EARTH_RADIUS': 6367.47021484375, 'IS_WIND_EARTH_REL': False}

	print(interfile['TT'].geoinfo)
	>> {'IPROJ': 0, 'PROJ': 'Cylindrical Equidistant (0)', 'STARTLOC': 'SWCORNER',
	'STARTLAT': 38.0, 'STARTLON': -130.0, 'DELTALAT': -0.25, 'DELTALON': 0.25}

	print(interfile['TT'].level)
	>> [100000.  97500.  95000.  92500.  90000.  87500.  85000.  82500.  80000.
	  77500.  75000.  70000.  65000.  60000.  55000.  50000.  45000.  40000.
	  35000.  30000.  25000.  20000.  17500.  15000.  12500.  10000.   7000.
	   5000.   3000.   2000.   1000.]
	   
	print(interfile['TT'].val)
	>> [[[289.86547852 289.89868164 289.85571289 ... 294.19750977 294.20141602
	   294.19360352]
	  ...
	  [278.23071289 277.93383789 277.69360352 ... 280.58032227 280.63500977
	   280.70727539]]
	   ...
	 [[234.77418518 234.80836487 234.85231018 ... 232.22828674 232.20973206
	   232.20582581]
	   ...
	  [213.70191956 213.95289612 214.19410706 ... 211.46461487 211.38844299
	   211.31422424]]]
	   
	print(interfile['TT'].val.shape)
	>> (31, 441, 329)


Creating Intermediate files
====================================
Before create intermediate files (cinter function), you must use the Geo-information(Geo0,Geo1,Geo3,Geo4,Geo5) and the type of variable you are using (V2d,V3d,V3dp,Vsl)


Geo-Information (Geo)
====================================
This funtions ares utilized to locate the information in the space. There are several kind of geo-info, it depends on projection of the original data:

- 0:  Cylindrical Equidistant (Lat/lon)
- 1:  Mercator projection
- 3:  Lambert conformal conic
- 4: Gaussian [global only] (Transverse mercator)
- 5:  Polar-stereographic projection

Geo0
-------------
**Cylindrical Equidistant (Lat/lon)**

:Geo0(stlat,stlon,dlat,dlon):

- stlat: SOUTH-WEST corner latitude of data (degrees north)
- stlon: SOUTH-WEST corner longitude of data (degrees north)
- dlat: latitude increment (degrees)
- dlon: longitude increment (degrees)


Geo1
-------------
**Mercator**

:Geo1(stlat,stlon,dx,dy,tlat1):

- stlat: SOUTH-WEST corner latitude of data (degrees north)
- stlon: SOUTH-WEST corner longitude of data (degrees north)
- dx: Grid spacing in x (Km)
- dy: Grid spacing in y (Km)
- tlat1: True latitude 1 of projection (degrees north)

Geo3
-------------
**Lambert conformal conic**

:Geo3(stlat,stlon,dx,dy,xloc,tlat1,tlat2,iswin):

- stlat: SOUTH-WEST corner latitude of data (degrees north)
- stlon: SOUTH-WEST corner longitude of data (degrees north)
- dx: Grid spacing in x (Km)
- dy: Grid spacing in y (Km)
- xloc: Center longitude of projection
- tlat1: True latitude 1 of projection (degrees north)
- tlat2: True latitude 2 of projection (degrees north)
- iswin: Earth[False] or source grid[True] rotated winds

Geo4
-------------
**Gaussian [global only] (Transverse mercator)**


:Geo4(stlat,stlon,nlats,dlon,iswin):

- stlat: SOUTH-WEST corner latitude of data (degrees north)
- stlon: SOUTH-WEST corner longitude of data (degrees north)
- nlats: Number of latitudes north of equator
- dlon: longitude increment (degrees)
- iswin: Earth[False] or source grid[True] rotated winds


Geo5
-------------
**Polar-stereographic**

:Geo5(stlat,stlon,dx,dy,xloc,tlat1,iswin):

- stlat: SOUTH-WEST corner latitude of data (degrees north)
- stlon: SOUTH-WEST corner longitude of data (degrees north)
- dx: Grid spacing in x (Km)
- dy: Grid spacing in y (Km)
- xloc: Center longitude of projection
- tlat1: True latitude 1 of projection (degrees north)
- iswin: Earth[False] or source grid[True] rotated winds


Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data

	# For Cylindrical equidistant
	# Read Geo-data (Latitudes and longitudes)
	lat = data.variables['Latitude'][:] # degrees north
	lon = data.variables['Longitude'][:] # degrees east

	# create winter Geo-information
	winter_geo = pwy.Geo0(lat,lon)



2D Field (V2d)
====================================
Surface variables

:V2d(name,field):

- name: WPS field name (see table)
- field: 2D array [lat,lon]


2D avalaible name fields
--------------------------
==================   ============   =============================       ============
Field name           Units          Description                         Notes
==================   ============   =============================       ============
PSFC                 Pa             Surface pressure
PMSL                 Pa             Mean sea-level pressure             
SKINTEMP             K              Skin temperature
SOILHGT              m              Soil height
TT                   K              2m air temperature
RH                   %              2m relative humidity                Not needed if SPECHUMD is avalaible
SPECHUMD             kg/kg          2m specific humidity                Not needed if RH is avalaible
UU                   m/s            10m wind u component 
VV                   m/s            10m wind v component 
LANDSEA              fraction       Land-sea mask                       0=water, 1=land
SST                  K              Sea surface temperature
SEAICE               fraction       Sea-ice fraction
SNOW                 kg/m^2         Water equivalent snow depth
TAVGSFC              K              Daily mean of surface air
==================   ============   =============================       ============

Some 2D fields are masked fields, so you must make sure to convert the missing values to numpy nan before create the pywinter 2d field.

Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data

	# Read 2D data (2m temperature, 10m U wind, 10m V wind)
	tp2m = data.variables['T2'][:,:]
	u10m = data.variables['U10'][:,:]
	v10m = data.variables['V10'][:,:]

	# Create winter 2D fields
	winter_t2m = pyw.V2d('TT',tp2m)
	winter_u10 = pyw.V2d('UU',u10m)
	winter_v10 = pyw.V2d('VV',u10m)


3D non-isobaric Field (V3d)
====================================
Vertical non-isobaric atmospehere variables

:V3d(name,field):

- name: WPS field name (see table)
- field: 3D array [lev,lat,lon]


3D non-isobaric avalaible name fields
--------------------------------------

==================   ============   =========================    ============
Field name           Units          Description                  Notes
==================   ============   =========================    ============
TT                   K              3D air temperature
RH                   %              3D relative humidity         Not needed if SPECHUMD is avalaible
SPECHUMD             kg/kg          3D specific humidity         Not needed if RH is avalaible
UU                   m/s            3D wind u component 
VV                   m/s            3D wind v component 
GHT                  m              3D geopotential height                  
PRESSURE             Pa             3D pressure                  Only needed for non-isobaric datasets
==================   ============   =========================    ============

Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data

	# Read 3D non-isobaric data (Pressure, U wind, V wind)
	press = data.variables['P'][:,:,:]
	uwind = data.variables['U'][:,:,:]
	vwind = data.variables['V'][:,:,:]

	# Create winter 3D non-isobaric fields
	winter_p = pyw.V3d('PRESSURE',press)
	winter_u = pyw.V3d('UU',u10m)
	winter_v = pyw.V3d('VV',u10m)


3D isobaric Field (V3dp)
====================================
Vertical isobaric atmospehere variables

:V3dp(name,field,plevs):

- name: WPS field name (see table)
- field: 3D array [plev,lat,lon]
- plevs: 1D arra of pressure levels (Pa)


3D isobaric avalaible name fields
----------------------------------

==================   ============   =========================    ============
Field name           Units          Description                  Notes
==================   ============   =========================    ============
TT                   K              3D air temperature
RH                   %              3D relative humidity         Not needed if SPECHUMD is avalaible
SPECHUMD             kg/kg          3D specific humidity         Not needed if RH is avalaible
UU                   m/s            3D wind u component 
VV                   m/s            3D wind v component 
GHT                  m              3D geopotential height                  
==================   ============   =========================    ============

Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data

	# Read 3D isobaric data (temperature, U wind, V wind)
	temp = data.variables['T'][:,:,:]
	uwind = data.variables['U'][:,:,:]
	vwind = data.variables['V'][:,:,:]

	# Read 3D pressure levels (hPa)
	plevs = data.variables['PLEV'][:] 

	# Create winter 3D isobaric fields
	winter_t = pyw.V3dp('TT',tp2m,plevs)
	winter_u = pyw.V3dp('UU',u10m,plevs)
	winter_v = pyw.V3dp('VV',u10m,plevs)



Soil Field (Vsl)
====================================
Soil level variables

:Vsl(name,field,levs):

- name: WPS field name (see table)
- field: 3D array [lev,lat,lon]
- levs: 1D list of string soil levels in cm


3D soil avalaible name fields
----------------------------------

==================   ============   =========================    ============================================
Field name           Units          Description                  Notes
==================   ============   =========================    ============================================
SM                   m^3/m^3        Soil moisture                'ttt' is layer top, 'bbb' is layer bottom
ST                   K              Soil temperature             'ttt' is layer top, 'bbb' is layer bottom
SOILM                kg^3/m^3       Soil moisture                'mmm' is tthe level depth
SOILT                K              Soil temperature             'mmm' is tthe level depth
==================   ============   =========================    ============================================

Is important to know that if you have ST or SM, you must indicate the layer as [*bbbttt*] (top layer - bottom layer in cm) and if you have SOILT or SOILM, you must indicate the level as [*mmm*] (level depth in cm). 

Aditionally if you have ST then SOILT is needed. if you have SM then SOILM is needed. 

Soil fields are used to be masked fields, so you must make sure to convert the missing values to numpy nan before create the pywinter soil field.

Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data

	# Read 3D soil data (temperature and moisture)
	soilt_lay = data.variables['SLTY'][:,:,:]
	soilm_lev = data.variables['SLTL'][:,:,:]

	# 3D soil layers and levels 
	slt_layer = ['000010','010040','040100','100200']
	slm_level = ['010','040','100','200']

	# Create winter 3D soil fields
	winter_soilt_layer = pyw.Vsl('ST',soilt_lay,sl_layer)
	winter_soitm_level = pyw.Vsl('SOILM',soilt_lev,sl_level)


Create intermediate files (cinter)
====================================
When you have all the pywinter meteorological variables, now is time to create the intermediate files, you must use **cinter** function and it is very easy.


:cinter(filen,date,geoinfo,varias,rout):

- filen: File name or prefix.
- date: Datetime of the file
- geoinfo: pywinter Geo object.
- varias: List of pywinter fields. 
- rout: Path where files will be created. ([optional] if not rout, files will be created in current folder)

Actually you can write just one file per time step. For creating several files you can use a loop, just remember that parameters filen, geoinfo, and rout will be static.



Example
-------------
.. code-block:: python

	import numpy as np
	import pywinter.winter as pyw
	import data_example as data


	## Read data
	############################################################

	# Read Geo-data
	lat = data.variables['Latitude'][:]
	lon = data.variables['Longitude'][:]
	dlat = np.abs(lat[1]-lat[0])
	dlon = np.abs(lon[1]-lon[0])
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
	winter_geo = pyw.Geo0(lat[0],lon[0],dlat,dlon)

	winter_t2m = pyw.V2d('TT',tp2m)
	winter_t = pyw.V3dp('TT',temp,plevs)
	winter_soilt_layer = pyw.Vsl('ST',soilt,sl_layer)

	# Listing fields
	total_fields = [winter_t2m,
	                winter_t,
	                winter_soilt_layer]

	# Out path
	path_out = '/home/Documents/intermediate_files/'

	# Write intermediate file
	pyw.cinter('FILE','1994-05-18_06',winter_geo,total_fields,path_out)




IMPORTANT WARNING
------------------

.. note::


	You can choose any prefix for the intermediate files though usually are called *FILE*, however the datetime must be chosen carefully, the files will be created with success but metgrid.exe sometimes can't read them:


	- If your original data has an exact hourly time resolution, for example: [ **05:00:00, 06:00:00** ], your intermediate files will be [ **FILE:YYYY-MM-DD_05**, **FILE:YYYY-MM-DD_06**} ]. **Minutes and seconds are not needed**, if you put them, it is probably metgrid.exe can't find the intermediate files.
	    
	- If your original data has less than one hour time resolution, for example: [ **05:00:00, 05:30:00** ], your intermediate files will be [ **FILE:YYYY-MM-DD\_05:00** , **FILE:YYYY-MM-DD\_05:30** ]. **Seconds are not needed**, if you put them, it is probably metgrid.exe can't find the intermediate files.

	- If your original data has a not-exact hourly time resolution, for example: [ **05:30:00, 06:30:00, 07:30:00** ], metgrid.exe won't find the intermediate files. You can trick the program if you set the **interval\_seconds** namelist parameter as **3601** and your intermediate files in this case will be [ **FILE:YYYY-MM-DD_05:30:00 , FILE:YYYY-MM-DD_06:30:01, FILE:YYYY-MM-DD_07:30:02** ]. This just works for short runs, because you must remember that in 60 hours the delay will be 1 minute with respect to your original data. This could be a bug or a simple restriction from metgrid.exe.


Final notes
====================================

Actually pywinter can't check if you write the file with all the necessary fields to run WRF, neither cannot check if your information is consistent, therefore is important you make sure the data is well before you create the files, you can get more information of how to create good files in the WRF User guide and WRF web tutorial, also you can check the intermediate files with WPS util programs:

- .../WRF/WPS/util/rd\_intermediate.exe: It reads fields into intermediate files and show them
- .../WRF/WPS/util/int2nc.exe: it converts intermediate files to netCDF format.


Documentation file
====================================

Full Documentation File :download:`pdf <https://github.com/dniloash/Pywinter/blob/master/docs/PYWINTER_documentation.pdf>`


More examples
====================================

More examples here :download:`link <https://github.com/dniloash/Pywinter/tree/master/examples>`


Contact
====================================

Pywinter have been tested with success for many cases but its yet an early version and it can be improved. if you have any problem our suggestion please contact to dniloash@gmail.com

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
