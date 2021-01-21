
import numpy as np
import os
import KreatE_inter_m_f as creattee_inter


#COMPILE
#f2py -c KreatE_inter_m_f.f90 -m KreatE_inter_m_f


class Interm:
    
    def __init__(self,filen,hdat,rout,geos,varto,gene={},geoinfo={},megamat=None,fields=None,units=None,descs=None,levels=None):
        
        self.filen = filen
        self.rout = rout
        self.geos = geos
        self.varto = varto
        self.gene = gene
        self.hdat = hdat
        self.geoinfo = geoinfo
        self.megamat = megamat

        self.fields = fields
        self.units = units
        self.descs = descs
        self.levels = levels


    def set_genatr(self):

        self.gene['VERSION'] = 5.0
        self.gene['XFCST'] = 0.0
        self.gene['MAP_SOURCE'] = 'PYWINTER' 


    def set_geoatr(self):

        self.geoinfo['IPROJ'] = self.geos.iproj
        self.geoinfo['NX'] = self.geos.nx
        self.geoinfo['NY'] = self.geos.ny
        self.geoinfo['STARTLOC'] = self.geos.stloc
        self.geoinfo['STARTLAT'] = self.geos.stlat
        self.geoinfo['STARTLON'] = self.geos.stlon
        self.geoinfo['EARTH_RADIUS'] = self.geos.earad
        self.geoinfo['IS_WIND_EARTH_REL'] = self.geos.iswin

        if isinstance(self.geos,Geo00):
            
            self.geoinfo['DELTALAT'] = self.geos.dlat
            self.geoinfo['DELTALON'] = self.geos.dlon

        elif isinstance(self.geos,Geo01):
            
            self.geoinfo['DX'] = self.geos.dx
            self.geoinfo['DY'] = self.geos.dy
            self.geoinfo['TRUELAT1'] = self.geos.tlat1

        elif isinstance(self.geos,Geo03):
            
            self.geoinfo['DX'] = self.geos.dx
            self.geoinfo['DY'] = self.geos.dy
            self.geoinfo['XLONC'] = self.geos.xlon
            self.geoinfo['TRUELAT1'] = self.geos.tlat1
            self.geoinfo['TRUELAT2'] = self.geos.tlat2
            
        elif isinstance(self.geos,Geo04):
            
            self.geoinfo['NLATS'] = self.geos.nlats
            self.geoinfo['DELTALON'] = self.geos.dlon

        elif isinstance(self.geos,Geo05):
            
            self.geoinfo['DX'] = self.geos.dx
            self.geoinfo['DY'] = self.geos.dy
            self.geoinfo['XLONC'] = self.geos.xlon
            self.geoinfo['TRUELAT1'] = self.geos.tlat1


    def set_varatr(self):

        self.megamat = np.zeros([len(self.varto),self.geoinfo['NX'],self.geoinfo['NY']])*np.nan
        self.levels = np.zeros(len(self.varto))*np.nan

        self.fields = []      
        self.units = []
        self.descs = []


        for i in range(len(self.varto)):
            
            variab = self.varto[i]
            self.fields.append(variab.fnam)
            self.units.append(variab.unit)
            self.descs.append(variab.desc)
            
            self.levels[i] = float(variab.levl)
            variab.field[~np.isfinite(variab.field)] = -1.e30
            self.megamat[i,:,:] = np.transpose(variab.field)


    def create_intermf(self):

        fname = self.filen
        hdate = self.hdat
        nx = self.geoinfo['NX']
        ny = self.geoinfo['NY']
        startlat = self.geoinfo['STARTLAT']
        startlon = self.geoinfo['STARTLON']

        if isinstance(self.geos,Geo00):
            
            deltalat = self.geoinfo['DELTALAT']
            deltalon = self.geoinfo['DELTALON']

        elif isinstance(self.geos,Geo01):
            
            dx = self.geoinfo['DX']
            dy = self.geoinfo['DY']
            tlat1 = self.geoinfo['TRUELAT1']


        elif isinstance(self.geos,Geo03):

            dx = self.geoinfo['DX']
            dy = self.geoinfo['DY']
            xlonc = self.geoinfo['XLONC']
            tlat1 = self.geoinfo['TRUELAT1']
            tlat2 = self.geoinfo['TRUELAT2']
            iswin = self.geoinfo['IS_WIND_EARTH_REL']
            


        elif isinstance(self.geos,Geo04):
            
            nlats = self.geoinfo['NLATS']
            deltalon = self.geoinfo['DELTALON']
            iswin = self.geoinfo['IS_WIND_EARTH_REL']


        elif isinstance(self.geos,Geo05):

            dx = self.geoinfo['DX']
            dy = self.geoinfo['DY']
            xlonc = self.geoinfo['XLONC']
            tlat1 = self.geoinfo['TRUELAT1']
            iswin = self.geoinfo['IS_WIND_EARTH_REL']
            
        
        ns = self.rout+self.filen+':'+self.hdat
        
        fci = self.hdat
        
        fnt = self.fields 
        fun = self.units 
        fds = self.descs
        flv = np.array(self.levels)
        
        nvars = len(self.megamat)
        
        va = self.megamat

        file1 = open('.inTer1','w')
        file2 = open('.inTer2','w')
        file3 = open('.inTer3','w')

        for h in range(len(fnt)):
            
            if h == len(fnt)-1:
                file1.write(fnt[h])
                file2.write(fun[h])
                file3.write(fds[h])
            else:
                file1.write(fnt[h]+ os.linesep)
                file2.write(fun[h]+ os.linesep)
                file3.write(fds[h]+ os.linesep)

        file1.close()
        file2.close()
        file3.close()

        fnt = '.inTer1'
        fun = '.inTer2'
        fds = '.inTer3'

        if isinstance(self.geos,Geo00):
            creattee_inter.crear_int0(startlat,startlon,deltalat,deltalon,fnt,fci,fun,fds,flv,ns,va)

        elif isinstance(self.geos,Geo01):
            creattee_inter.crear_int1(startlat,startlon,dx,dy,tlat1,fnt,fci,fun,fds,flv,ns,va)

        elif isinstance(self.geos,Geo03):
            creattee_inter.crear_int3(startlat,startlon,dx,dy,xlonc,tlat1,tlat2,iswin,fnt,fci,fun,fds,flv,ns,va)

        elif isinstance(self.geos,Geo04):
            creattee_inter.crear_int4(startlat,startlon,nlats,deltalon,iswin,fnt,fci,fun,fds,flv,ns,va)

        elif isinstance(self.geos,Geo05):
            creattee_inter.crear_int5(startlat,startlon,dx,dy,xlonc,tlat1,iswin,fnt,fci,fun,fds,flv,ns,va)

        print(fname+':'+hdate)


##########################################################

class Var:

    def __init__(self,fnam,unit,desc,levl,field):

        self.fnam = fnam
        self.unit = unit
        self.desc = desc
        self.levl = levl
        self.field = field
        

class Var2d(Var):

    def __init__(self,var2d,fnam='',unit='',desc='',levl='',field=[]):

        Var.__init__(self,fnam,unit,desc,levl,field)
        self.var2d = var2d

    def set_atr(self):

        ffield = self.var2d.field
        ffield[~np.isfinite(ffield)]= np.nan

        try:
            ffield[ffield.mask]= np.nan
        except:
            pass

        atrib = self.var2d.__dict__

        if (len(atrib['des']) != 0) and (len(atrib['uni']) != 0) and (len(atrib['lev']) != 0):

            self.fnam = self.var2d.name
            self.unit = self.var2d.uni
            self.desc = self.var2d.des
            self.levl = self.var2d.lev
            self.field = ffield

        else:

            nom,uni,lvl = self.var2d.idvar()
            
            self.fnam = self.var2d.name
            self.unit = uni
            self.desc = nom
            self.levl = lvl
            self.field = ffield
        



class Var3d(Var):

    def __init__(self,var3d,fnam='',unit='',desc='',levl=[],field=[]):

        Var.__init__(self,fnam,unit,desc,levl,field)
        self.var3d = var3d

    def set_atr(self):

        ffield = self.var3d.field
        ffield[~np.isfinite(ffield)]= np.nan

        try:
            ffield[ffield.mask]= np.nan
        except:
            pass

        nom,uni = self.var3d.idvar()

        nlevs = len(self.var3d.field)

        levs = np.arange(1.0,nlevs+1,1.0) 

        self.fnam = self.var3d.name
        self.unit = uni
        self.desc = nom
        self.levl = levs[::-1]
        self.field = ffield




class Var3dp(Var):

    def __init__(self,var3dp,fnam='',unit='',desc='',levl=[],field=[]):

        Var.__init__(self,fnam,unit,desc,levl,field)
        self.var3dp = var3dp


    def set_atr(self):

        ffield = self.var3dp.field
        ffield[~np.isfinite(ffield)]= np.nan

        try:
            ffield[ffield.mask]= np.nan
        except:
            pass     

        nom,uni = self.var3dp.idvar()

        self.fnam = self.var3dp.name
        self.unit = uni
        self.desc = nom
        self.levl = self.var3dp.plevs
        self.field = ffield

        
        

class Varsl(Var):

    def __init__(self,varsl,fnam='',unit='',desc='',levl='',slev=[],field=[]):

        Var.__init__(self,fnam,unit,desc,levl,field)
        self.varsl = varsl
        self.slev = slev

    def set_atr(self):

        ffield = self.varsl.field
        ffield[~np.isfinite(ffield)]= np.nan

        try:
            ffield[ffield.mask]= np.nan
        except:
            pass      

        nom,uni,lvl = self.varsl.idvar()

        self.fnam = self.varsl.name
        self.unit = uni
        self.desc = nom
        self.levl = lvl
        self.slev = self.varsl.levs
        self.field = ffield


##########################################################

class Vuser:
    def __init__(self):
        pass
        
class V2d(Vuser):

    def __init__(self,name,field,des='',uni='',lev=''):
        self.name = name
        self.field = field
        self.des = des
        self.uni = uni
        self.lev = lev

        self.verifdim()

    def verifdim(self):
         
        if len(self.field.shape) != 2:
            raise Error('field must a 2D array')

    def idvar(self):

        if self.name == 'PSFC':
            nom = 'Surface pressure'
            uni = 'Pa'
            lev = '200100'

        elif self.name == 'PMSL':
            nom = 'Mean sea-level pressure'
            uni = 'Pa'
            lev = '201300'
            
        elif self.name == 'SKINTEMP':
            nom = 'Skin temperature'
            uni = 'K'
            lev = '200100'       

        elif self.name == 'SOILHGT':
            nom = 'Soil height'
            uni = 'm'
            lev = '200100'

        elif self.name == 'TT':
            nom = '2-meter air temperature'
            uni = 'K'
            lev = '200100'        

        elif self.name == 'RH':
            nom = '2-meter relative humidity'
            uni = '%'
            lev = '200100'

        elif self.name == 'SPECHUMD':
            nom = '2-meter specific humidity'
            uni = 'kg kg-1'
            lev = '200100'

        elif self.name == 'UU':
            nom = '10-meter wind u-component'
            uni = 'kg kg-1'
            lev = '200100'

        elif self.name == 'VV':
            nom = '10-meter wind v-component'
            uni = 'kg kg-1'
            lev = '200100'

        elif self.name == 'LANDSEA':
            nom = 'Land sea mask'
            uni = 'fraction'
            lev = '200100'

        elif self.name == 'SST':
            nom = 'Sea surface temperature'
            uni = 'K'
            lev = '201300'

        elif self.name == 'SEAICE':
            nom = 'Sea-ice-fraction'
            uni = 'fraction'
            lev = '200100'

        elif self.name == 'SNOW':
            nom = 'Water equivalent of Accum snow depth'
            uni = 'kg m-2'
            lev = '200100'
            
        elif self.name == 'TAVGSFC':
            nom = 'Daily mean of surface air temperature'
            uni = 'K'
            lev = '200100'

        return nom,uni,lev

            
class V3d(Vuser):

    def __init__(self,name,field):
        self.name = name
        self.field = field

        self.verifdim()

    def verifdim(self):
         
        if len(self.field.shape) != 3:
            raise Error('field must a 3D array')

    def idvar(self):

        if self.name == 'TT':
            nom = '3-d air temperature'
            uni = 'K'
        
        elif self.name == 'RH':
            nom = '3-d relative humidity'
            uni = '%'

        elif self.name == 'SPECHUMD':
            nom = '3-d specific humidity'
            uni = 'kg kg-1'

        elif self.name == 'UU':
            nom = '3-d wind u-component'
            uni = 'm s-1'

        elif self.name == 'VV':
            nom = '3-d wind v-component'
            uni = 'm s-1'

        elif self.name == 'GHT':
            nom = '3-d geopotential height '
            uni = 'm'

        elif self.name == 'PRESSURE':
            nom = '3-d pressure '
            uni = 'Pa'

        return nom,uni


class V3dp(Vuser):

    def __init__(self,name,field,plevs):
        self.name = name
        self.field = field
        self.plevs = plevs

        self.verifdim()

    def verifdim(self):
         
        if len(self.field.shape) != 3:
            raise Error('field must a 3D array')

        if self.plevs.shape[0] != self.field.shape[0]:
            raise Error('number of pressure levels must match with 3D field')
            

    def idvar(self):
        
        if self.name == 'TT':
            nom = '3-d air temperature'
            uni = 'K'
            
        elif self.name == 'RH':
            nom = '3-d relative humidity'
            uni = '%'

        elif self.name == 'SPECHUMD':
            nom = '3-d specific humidity'
            uni = 'kg kg-1'

        elif self.name == 'UU':
            nom = '3-d wind u-component'
            uni = 'm s-1'

        elif self.name == 'VV':
            nom = '3-d wind v-component'
            uni = 'm s-1'

        elif self.name == 'GHT':
            nom = '3-d geopotential height '
            uni = 'm'

        return nom,uni


class Vsl(Vuser):

    def __init__(self,name,field,levs):
        self.name = name
        self.field = field
        self.levs = levs

        self.verifdim()
        self.veriflevs()


    def verifdim(self):
         
        if len(self.field.shape) != 3:
            raise Error('field must a 3D array')

        if len(self.levs) != self.field.shape[0]:
            raise Error('number of levels must match with 3D field')


    def veriflevs(self):

        for lev in self.levs:
            if isinstance(lev,str):
                for i in lev:
                    try:
                        a = int(i)
                    except:
                        raise Error('soil levels must be a numerical string')
                        
            else:
                raise Error('soil levels must be strings')
                

        if self.name == 'ST' or self.name == 'SM':
            for i in self.levs:
                if len(i) != 6:
                    raise Error('soil levels for ST or SM must have bbbttt format')
                
        elif self.name == 'SOILT' or self.name == 'SOILM':
            for i in self.levs:
                if len(i) != 3:
                    raise Error('soil levels for SOILT or SOILM must have mmm format')


    def idvar(self):

        if self.name == 'SM':
            nom = 'Soil moisture '
            uni = 'm3 m-3'
            lev = '200100'
            
        elif self.name == 'ST':
            nom = 'Soil temperature '
            uni = 'K'
            lev = '200100'
            
        elif self.name == 'SOILM':
            nom = 'Soil moisture '
            uni = 'kg m-3'
            lev = '200100'
            
        elif self.name == 'SOILT':
            nom = 'Soil temperature '
            uni = 'K'
            lev = '200100'

        return nom,uni,lev


##########################################################
    
class Geoinfo:

    def __init__(self,iproj,nx,ny,stloc,stlat,stlon,earad,iswin):

        self.iproj = iproj
        self.stloc = stloc
        self.stlat = stlat
        self.stlon = stlon
        self.earad = earad
        self.iswin = iswin


class Geouser:
    def __init__(self):
        pass

        
class Geo00(Geoinfo):

    def __init__(self,geoin,iproj=0,nx=0,ny=0,stloc='',stlat=0.0,stlon=0.0,earad=0.0,
                 iswin=False,dlat=0.0,dlon=0.0):

        Geoinfo.__init__(self,iproj,nx,ny,stloc,stlat,stlon,earad,iswin)
        self.geoin = geoin
        self.dlat = dlat
        self.dlon = dlon


    def set_atr(self):

        cnx = len(self.geoin.lons)
        cny = len(self.geoin.lats)

        mid1 = int(len(self.geoin.lats)/2)
        mid2 = int(len(self.geoin.lons)/2)
        
        cstlat = self.geoin.lats[0]
        cstlon = self.geoin.lons[0]
        cdlat = np.abs(self.geoin.lats[mid1+1] - self.geoin.lats[mid1])
        cdlon = np.abs(self.geoin.lons[mid2+1] - self.geoin.lons[mid2])

        self.iproj = 0
        self.nx = cnx
        self.ny = cny
        self.stloc = 'SWCORNER'
        self.stlat = cstlat
        self.stlon = cstlon
        self.earad = 6367.470215
        self.iswin = False
        self.dlat = cdlat
        self.dlon = cdlon


class Geo0(Geouser):

    def __init__(self,lats,lons):
        
        self.lats = lats
        self.lons = lons

        self.verifdim()


    def verifdim(self):

        if len(self.lats.shape) != 1 and len(self.lats.shape) != 1:
            raise Error('latitude and longitude must de 1D arrays')


#-----------------------------------------------------------------------------------------------


class Geo01(Geoinfo):

    def __init__(self,geoin,iproj=0,nx=0,ny=0,stloc='',stlat=0.0,stlon=0.0,earad=0.0,
                 iswin=False,dx=0.0,dy=0.0,tlat1=0.0):

        Geoinfo.__init__(self,iproj,nx,ny,stloc,stlat,stlon,earad,iswin)
        self.geoin = geoin
        self.dx = dx
        self.dy = dy
        self.tlat1 = tlat1



    def set_atr(self):

        cnx = len(self.geoin.lats)
        cny = len(self.geoin.lons)

        
        cstlat = self.geoin.lats[0]
        cstlon = self.geoin.lons[0]
        
        cdx = self.geoin.dx/1000.
        cdy = self.geoin.dy/1000.

        tlat1 = self.geoin.tlat1

        self.iproj = 0
        self.nx = cnx
        self.ny = cny
        self.stloc = 'SWCORNER'
        self.stlat = cstlat
        self.stlon = cstlon
        self.earad = 6367.470215
        self.iswin = False
        self.dx = cdx
        self.dy = cdy
        self.tlat1 = tlat1



class Geo1(Geouser):

    def __init__(self,lats,lons,dx,dy,tlat1):

        self.lats = lats
        self.lons = lons
        self.dx = dx
        self.dy = dy
        self.tlat1 = tlat1
        
        self.verifdim()


    def verifdim(self):

        if len(self.lats.shape) != 1 and len(self.lats.shape) != 1:
            raise Error('latitude and longitude must be 1D arrays')


#-----------------------------------------------------------------------------------------------

class Geo03(Geoinfo):

    def __init__(self,geoin,iproj=0,nx=0,ny=0,stloc='',stlat=0.0,stlon=0.0,earad=0.0,
                 iswin=False,dx=0.0,dy=0.0,xlon=0.0,tlat1=0.0,tlat2=0.0):

        Geoinfo.__init__(self,iproj,nx,ny,stloc,stlat,stlon,earad,iswin)
        self.geoin = geoin
        self.dx = dx
        self.dy = dy
        self.xlon = xlon
        self.tlat1 = tlat1
        self.tlat2 = tlat2


    def set_atr(self):

        cnx = len(self.geoin.lats)
        cny = len(self.geoin.lons)

        
        cstlat = self.geoin.lats[0]
        cstlon = self.geoin.lons[0]
        
        cdx = self.geoin.dx/1000.
        cdy = self.geoin.dy/1000.

        xl = self.geoin.xloc

        tlat1 = self.geoin.tlat1
  
        tlat2 = self.geoin.tlat2

        iswin = self.geoin.iswin

        self.iproj = 0
        self.nx = cnx
        self.ny = cny
        self.stloc = 'SWCORNER'
        self.stlat = cstlat
        self.stlon = cstlon
        self.earad = 6367.470215
        self.iswin = iswin
        self.dx = cdx
        self.dy = cdy
        self.xlon = xl
        self.tlat1 = tlat1
        self.tlat2 = tlat2




class Geo3(Geouser):

    def __init__(self,lats,lons,dx,dy,xloc,tlat1,tlat2,iswin):

        self.lats = lats
        self.lons = lons
        self.dx = dx
        self.dy = dy
        self.xloc = xloc
        self.tlat1 = tlat1
        self.tlat2 = tlat2
        self.iswin = iswin

        self.verifdim()


    def verifdim(self):

        if len(self.lats.shape) != 1 and len(self.lats.shape) != 1:
            raise Error('latitude and longitude must de 1D arrays')
#-----------------------------------------------------------------------------------------------


class Geo04(Geoinfo):

    def __init__(self,geoin,iproj=0,nx=0,ny=0,stloc='',stlat=0.0,stlon=0.0,earad=0.0,
                 iswin=False,nlats=0.0,dlon=0.0):

        Geoinfo.__init__(self,iproj,nx,ny,stloc,stlat,stlon,earad,iswin)
        self.geoin = geoin
        self.nlats = nlats
        self.dlon = dlon


    def set_atr(self):

        cnx = len(self.geoin.lons)
        cny = len(self.geoin.lats)

        mid1 = int(len(self.geoin.lats)/2)
        mid2 = int(len(self.geoin.lons)/2)
        
        cstlat = self.geoin.lats[0]
        cstlon = self.geoin.lons[0]
        nlats = self.geoin.nlats
        cdlon = np.abs(self.geoin.lons[mid2+1] - self.geoin.lons[mid2])

        iswin = self.geoin.iswin

        self.iproj = 0
        self.nx = cnx
        self.ny = cny
        self.stloc = 'SWCORNER'
        self.stlat = cstlat
        self.stlon = cstlon
        self.earad = 6367.470215
        self.iswin = iswin
        self.nlats = nlats
        self.dlon = cdlon


class Geo4(Geouser):

    def __init__(self,lats,lons,nlats,iswin):
        
        self.lats = lats
        self.lons = lons
        self.nlats = nlats
        self.iswin = iswin

        self.verifdim()


    def verifdim(self):

        if len(self.lats.shape) != 1 and len(self.lats.shape) != 1:
            raise Error('latitude and longitude must de 1D arrays')  
#-----------------------------------------------------------------------------------------------

class Geo05(Geoinfo):

    def __init__(self,geoin,iproj=0,nx=0,ny=0,stloc='',stlat=0.0,stlon=0.0,earad=0.0,
                 iswin=False,dx=0.0,dy=0.0,xlon=0.0,tlat1=0.0):

        Geoinfo.__init__(self,iproj,nx,ny,stloc,stlat,stlon,earad,iswin)
        self.geoin = geoin
        self.dx = dx
        self.dy = dy
        self.xlon = xlon
        self.tlat1 = tlat1


    def set_atr(self):

        cnx = len(self.geoin.lats)
        cny = len(self.geoin.lons)

        
        cstlat = self.geoin.lats[0]
        cstlon = self.geoin.lons[0]
        
        cdx = self.geoin.dx/1000.
        cdy = self.geoin.dy/1000.

        xl = self.geoin.xloc

        tlat1 = self.geoin.tlat1
  
        tlat2 = self.geoin.tlat2

        iswin = self.geoin.iswin

        self.iproj = 0
        self.nx = cnx
        self.ny = cny
        self.stloc = 'SWCORNER'
        self.stlat = cstlat
        self.stlon = cstlon
        self.earad = 6367.470215
        self.iswin = iswin
        self.dx = cdx
        self.dy = cdy
        self.xlon = xl
        self.tlat1 = tlat1




class Geo5(Geouser):

    def __init__(self,lats,lons,dx,dy,xloc,tlat1,iswin):

        self.lats = lats
        self.lons = lons
        self.dx = dx
        self.dy = dy
        self.xloc = xloc
        self.tlat1 = tlat1
        self.iswin = iswin

        self.verifdim()


    def verifdim(self):

        if len(self.lats.shape) != 1 and len(self.lats.shape) != 1:
            raise Error('latitude and longitude must de 1D arrays')  
        
##########################################################

class Error(Exception):
    def __init__(self,expre):
        self.expre = expre

##########################################################


def chek_date(nome,dato):
    name = nome
    date = dato
    maxlen = 19
    
    if isinstance(name,str) or isinstance(name,str):
        pass
    else:
        raise Error('name and datetime must be strings')
    
    if len(date) > maxlen:
        raise Error('datetime format must be YYYY-MM-DD_hh:mm:ss') 
        


def cinter(filen,date,geoinfo,varias,rout=''):

    chek_date(filen,date)
    
    if rout == '':
        pass
    elif rout[-1] != '/':
        rout = rout + '/'

    for i in varias:
        if isinstance(i,Vuser):
            pass
        else: ('Variable elements must be winter.V type')
            

    if isinstance(geoinfo,Geo0):
        geoo = Geo00(geoinfo)

    elif isinstance(geoinfo,Geo1):
        geoo = Geo01(geoinfo)

    elif isinstance(geoinfo,Geo3):
        geoo = Geo03(geoinfo)

    elif isinstance(geoinfo,Geo4):
        geoo = Geo04(geoinfo)

    elif isinstance(geoinfo,Geo5):
        geoo = Geo05(geoinfo)

    else:
        if isinstance(geoinfo,Geouser):
            pass
        else:
            raise Error('geo-information element is not winter.Geo type')

    geoo.set_atr()

    variasf = []
    
    for i in range(len(varias)):
        varia = varias[i]

        if isinstance(varia,V2d):

            varia2d = Var2d(varia)
            varia2d.set_atr()

            variasf.append(varia2d)




        elif isinstance(varia,V3d):

            varia3d = Var3d(varia)
            varia3d.set_atr()


            nlevs = len(varia3d.levl)

            for j in range(nlevs):

                lev = varia3d.levl[j]
                levi = str(int(lev))
                field2 = varia3d.field[j,:,:]
                
                vv2d = V2d(varia3d.fnam,field2)
                
                vvv2d = Var2d(vv2d,varia3d.fnam,varia3d.unit,varia3d.desc,levi,field2)

                variasf.append(vvv2d)

                

        elif isinstance(varia,V3dp):

            varia3dp = Var3dp(varia)
            varia3dp.set_atr()


            nlevs = len(varia3dp.levl)

            for j in range(nlevs):

                lev = varia3dp.levl[j]
                levi = str(int(lev*100))
                field2 = varia3dp.field[j,:,:]

                
                vv2d = V2d(varia3dp.fnam,field2)
                
                vvv2d = Var2d(vv2d,varia3dp.fnam,varia3dp.unit,varia3dp.desc,levi,field2)

                variasf.append(vvv2d)



        elif isinstance(varia,Vsl):

            variasl = Varsl(varia)
            variasl.set_atr()

            soilplev = variasl.levl


            nlevs = len(variasl.slev)


            for j in range(nlevs):

                lev = variasl.slev[j]

                slfname = variasl.fnam+lev 
                field2 = variasl.field[j,:,:]

                
                vv2d = V2d(slfname,field2)
                
                vvv2d = Var2d(vv2d,slfname,variasl.unit,variasl.desc,soilplev,field2)

                variasf.append(vvv2d)



    interfile = Interm(filen,date,rout,geoo,variasf)

    interfile.set_geoatr()
    interfile.set_varatr()
    interfile.create_intermf()
