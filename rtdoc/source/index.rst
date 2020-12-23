.. pywinter documentation master file, created by
   sphinx-quickstart on Sat Jun  6 13:42:57 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pywinter
====================================
Python WRF Intermediate Files

Pywinter is a Python3 library designed to create files in WPS intermediate file format (one step back to execute metgrid.exe). Usually you don't need to create the intermediate files by your own because that is the function of ungrib.exe, but sometimes you don't have the meteorological data in GRIB1 or GRIB2 format, this is the case of NetCDF data. Pywinter allows to create intermediate files by a simple way, the only prerequisite is to read the data with python and next you can convert it on intermediate file format with pywinter.


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

:Geo0(lats,lons):

- lats: 1D array of latitudes (degrees north)
- lons: 1D array of latitudes (degrees east)


Geo1
-------------
**Mercator**

:Geo1(lats,lons,dx,dy,tlat1):

- lats: 1D vector of latitudes (degrees north)
- lons: 1D vector of latitudes (degrees east)
- dx: Grid spacing in x (m)
- dy: Grid spacing in y (m)
- tlat1: True latitude 1 of projection (degrees north)

Geo3
-------------
**Lambert conformal conic**

:Geo3(lats,lons,dx,dy,xloc,tlat1,tlat2,iswin):

- lats: 1D vector of latitudes (degrees north)
- lons: 1D vector of latitudes (degrees east)
- dx: Grid spacing in x (m)
- dy: Grid spacing in y (m)
- xloc: Center longitude of projection
- tlat1: True latitude 1 of projection (degrees north)
- tlat2: True latitude 2 of projection (degrees north)
- iswin: Earth[False] or source grid[True] rotated winds

Geo4
-------------
**Gaussian [global only] (Transverse mercator)**


:Geo4(lats,lons,nlats,iswin):

- lats: 1D vector of latitudes (degrees north)
- lons: 1D vector of latitudes (degrees east)
- nlats: Number of latitudes north of equator
- iswin: Earth[False] or source grid[True] rotated winds


Geo5
-------------
**Polar-stereographic**

:Geo5(lats,lons,dx,dy,xloc,tlat1,iswin):

- lats: 1D vector of latitudes (degrees north)
- lons: 1D vector of latitudes (degrees east)
- dx: Grid spacing in x (m)
- dy: Grid spacing in y (m)
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
- plevs: 1D arra of pressure levels (hPa)


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
Vertical isobaric atmospehere variables

:Vsl(name,field,levs):

- name: WPS field name (see table)
- field: 3D array [lev,lat,lon]
- levs: 1D list of string soil levels in cm


3D isobaric avalaible name fields
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
	pwy.cinter('FILE','1994-05-18_06',winter_geo,total_fields,path_out)



IMPORTANT WARNING
------------------

.. note::


	You can choose any prefix for the intermediate files though usually are called *FILE*, however the datetime must be chosen carefully, the files will be created with success but metgrid.exe sometimes can't read them:


	- If your original data has an exact hourly time resolution, for example: [ **05:00:00, 06:00:00** ], your intermediate files will be [ **FILE:YYYY-MM-DD_05**, **FILE:YYYY-MM-DD_06**} ]. **Minutes and seconds are not needed**, if you put them, it is probably metgrid.exe can't find the intermediate files.
	    
	- If your original data has less than one hour time resolution, for example: [ **05:00:00, 05:30:00** ], your intermediate files will be [ **FILE:YYYY-MM-DD\_05:00** , **FILE:YYYY-MM-DD\_05:30** ]. **Seconds are not needed**, if you put them, it is probably metgrid.exe can't find the intermediate files.

	- If your original data has a not-exact hourly time resolution, for example: [ **05:30:00, 06:30:00, 07:30:00** ], metgrid.exe won't find the intermediate files. You can trick the program if you set the **interval\_seconds** namelist parameter as **3601** and your intermediate files in this case will be [ **FILE:YYYY-MM-DD_05:30:00 , FILE:YYYY-MM-DD_06:30:01, FILE:YYYY-MM-DD_07:30:02** ]. This just works for short runs, because you must remember that in 60 hours the delay will be 1 minute with respect to your original data. This could be a bug or a simple restriction from metgrid.exe.



Documentation file
====================================

Full Documentation File :download:`pdf <https://github.com/dniloash/Pywinter/blob/master/docs/PYWINTER_documentation.pdf>`


More examples
====================================

More examples here :download:`link <https://github.com/dniloash/Pywinter/tree/master/examples>`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
