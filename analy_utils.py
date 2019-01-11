# -*- coding: utf-8 -*-
from netCDF4 import Dataset
import numpy as np
import os
import re
from alt_colormaps import viridis
import matplotlib.pylab as plt
#collection of utility functions for Analysis of COSMO model output
#PRE: 	valid full path of nc-file FULLPATH and valid VAR (list), which exists in FULLPATH
#POST: 	opens FULLPATH.nc, extracts VAR , closes FULLPATH.nc, returns dict with dict.keys() equal (time,VAR).
#		i.e. time is extracted in any case
nhalo=3

vcoordvec = [22000.0, 21000.0, 20028.6, 19085.4, 18170.0, 17282.1, 16421.4, 15587.5,
        14780.0, 13998.6, 13242.9, 12512.5, 11807.1, 11126.4, 10470.0,  9837.5,
        9228.6,  8642.9,  8080.0,  7539.6,  7021.4,  6525.0,  6050.0,  5596.1,
        5162.9,  4750.0,  4357.1,  3983.9,  3630.0,  3295.0,  2978.6,  2680.4,
        2400.0,  2137.1,  1891.4,  1662.5,  1450.0,  1253.6,  1072.9,  907.5,
        757.1,   621.4,   500.0,   392.5,   298.6,   217.9,   150.0,   94.6,
        51.4,    20.0,     0.0]
""" Routine to remove axes: left and top.
    http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html  
"""
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
def getha(oroname):
    ah=oroname.replace('a','h').split('h')[1:]
    h,a=int(ah[0]),int(ah[1])
    return h,a 


def add_customcolorbar(f,cf,mylabel,orient='vertical',pos='left'):
    """ Add space in figure and add custom colorbar there

    Arguments:
    -- f        figure object
    -- cf       contour field
    -- mylabel  label for colorbar
    
    Returns:
    -- cb       colorbar object
    
    """
    if orient=='vertical':
        if pos=='right':
            cbaxes = f.add_axes([0.91,0.1,0.03,0.8])    
        else:
            cbaxes = f.add_axes([0.05,0.1,0.03,0.8]) 
    else:
        cbaxes = f.add_axes([0.2,0.04,0.5,0.03])
        
    cb = plt.colorbar(cf,orientation=orient,cax=cbaxes,
                         label=mylabel)
    return cbaxes

def decorate_ax(ax,mytextl,mytextm,mytextr,fnsz=6,ypos=0.94):
    """ Adds three text boxes at the top of the axes frame (left, middle, right)
    
    Arguments:
    -- ax:         matplotlib.plt.axes object
    -- mytextl:    text string to be placed at top left
    -- mytextm:    text string to be placed at top center
    -- mytextr:    text string to be placed at top right
    
    """    

    
    ax.text(0.1,ypos,mytextl,ha='center', va='center', 
            transform=ax.transAxes, fontsize = fnsz,fontweight='bold',
            backgroundcolor = 'white') 
    ax.text(0.5,ypos,mytextm,ha='center', va='center', 
            transform=ax.transAxes, fontsize = fnsz,fontweight='bold',
            backgroundcolor = 'white')  
    ax.text(0.8,ypos,mytextr,ha='center', va='center', 
            transform=ax.transAxes, fontsize = fnsz,fontweight='bold',
            backgroundcolor = 'white') 


def interpolatecosmofield(field,zmodel,zinterp):
    """ Interpolates the field defined at levels of vertical grid of the model (in zmodel) 
    to values of zinterp."""
    nx,ny,nz = field.shape
    



class Ensemble:
    def __init__(self,PATH,*args):
        allruns=os.listdir(PATH)
        self.analysispath=PATH
        self.runnames=(run for run in allruns if (re.search(args[0],run)))
        self.members={run:() for run in self.runnames}
        self.number=len(self.runnames)
        
    def calc_statistics(self, var='TOT_PREC'):
        for run in self.runnames:
            totprecfields=ncvar2pyvar(self.path+run,'TOT_PREC')
            self.members[run]=np.cumsum(totprecfields[:,nhalo:-nhalo,nhalo:-nhalo],axis=0)
            
        self.ens_stat_funs={'min':np.min,'max':np.max,'mean':np.mean} #'median':np.median,
        self.ens_stats={q:() for q in (self.ens_stat_funs).keys()}
        
            
            
    def plot_envelope(self):
        pass
            
        
        



class time:
	def __init__(self, cform):
		if not cform.startswith('lfff'):
			assert 0
		cform=cform[4:-3]
		self.dd=cform[0:2]
		self.hh=cform[2:4]
		self.mm=cform[4:6]
		self.ss=cform[6:8]
	def printtime(self,format='hhmm' ):
		if format=='hhmm':
			print self.hh+':'+self.mm



def ncvars2pyvars(FULLPATH,varnames,appendtime=False):
	if not type(varnames)==list:
		varnames=list(varnames)
	ncfile=Dataset(FULLPATH, 'r')
	#make sure that 'time' is extracted in any case
	
	if appendtime and ('time' not in varnames):
		varnames.append('time')
	#varsdict={varname: ncfile.variables[varname][:] for varname in varnames}
	outvar=ncfile.variables[varname][:]
	ncfile.close()
	return outvar
	
#PRE: like ncvars2pyvars, but only single vars
	
def ncvar2pyvar(FULLPATH,varname):
	ncfile=Dataset(FULLPATH, 'r')#format='NETCDF43_64BIT')
	outvar=ncfile.variables[varname][:]
	ncfile.close()
	return outvar


#PRE: valid ncfile containing var
#POST:mean vertical profile of vertprof; OUT (z,vertprof)

def vertprofile(infullpathnc,var,vertcoor='z',ilocal=False,clipx=100):
	if vertcoor=='z':
		z=ncvar2pyvar(infullpathnc,'HHL')[0,::-1,:,:]
	if vertcoor=='p':
		z=ncvar2pyvar(infullpathnc,'P')[0,::-1,:,:]
	z=np.mean(z,axis=(1,2))			
	if ilocal:
		vertprof=ncvar2pyvar(infullpathnc,var)[0,::-1,(clipx+nhalo):-(clipx+nhalo),(clipx+nhalo):-(clipx+nhalo)]
	else:
		vertprof=ncvar2pyvar(infullpathnc,var)[0,::-1,nhalo:-nhalo,nhalo:-nhalo]
	vertprof=np.mean(vertprof,axis=(1,2))
	return z,vertprof
	
	
		
	
#PRE: string with cosmo output format lfffddhhmmss.nc
#POST: returns dictionary {'dd':dd,'hh':hh,'mm':mm,'ss':ss}
	
	
def cosoutform2time(cosform): 
    cosform=cosform[4:-3] #clip lfffddhhmmss.nc
    format='ddhhmmss'
    return {format[i:i+2]:int(cosform[i:i+2]) for i in range(0,8,2)}
    
#PRE:	takes string containing a n[N,S]m[E,W] expression
#POST:	returns coordinates [n,m]

def casename2coor(case):
	ofs=re.compile('offset')
	expr=ofs.split(case)[1] #returns n[N,S]m[E,W]...
	ofs=re.compile('_')
	expr=ofs.split(expr)[0] #returns n[N,S]m[E,W]
	nr=re.compile('[N,S,E,W]')
	coor=nr.split(expr)[0:2] #returns n,m
	sgn1=sgn2=1
	if re.match('S',expr):
		sgn1=-1 
	if re.match('W',expr):
		sgn2=-1 
	coor[0]=sgn1*int(coor[0])
	coor[1]=sgn2*int(coor[1])
	return coor




def get_vertcoord(hsurf,nz=51):
    "pass surface height, return vertical coordinates"
    ivctype = 2
    vcflat = 6000.
    zz_top = 9999.

    _ks = np.arange(0,nz)
    zak = _ks * 0.
    zbk = _ks * 0.


    kflat = 0
    # Inverse coordinate transfromation to obtain zak and zbk
    for k in _ks:
        if vcoordvec[k] >= vcflat:
            zak[k] = vcoordvec[k]
            zbk[k] = 0.
            kflat = k
        else:
            zak[k] = vcoordvec[k]
            zbk[k] = (vcflat -vcoordvec[k])/vcflat

    # Calcualte hsurf
    vcoordath = zak+zbk*hsurf
    return vcoordath


def HHL_creator(Hm=250.,nx=300,ny=300,nz=51,a=20.,ay=False,surftopo='cos2'):
    """ Port of Fortran code in vgrid_refatm_utils.f90
    HHL gives AGL levels of every grid point.

    Namelist parameters:
    ----
    ivctype = 2,
    zspacing_type = 'vcoordvec', ! sub-type of coordinate spec.
    exp_galchen = 2.6,        ! exponent in the Gal-Chen formula
    vcflat = 6000.0,          ! height, above which coordinate levels become flat [m]
    zz_top = 9999.0,          ! height of model top, if it has to be specified explicitly [m]
    vcoordvec =
    22000.0, 21000.0, 20028.6, 19085.4, 18170.0, 17282.1, 16421.4, 15587.5,
    14780.0, 13998.6, 13242.9, 12512.5, 11807.1, 11126.4, 10470.0,  9837.5,
     9228.6,  8642.9,  8080.0,  7539.6,  7021.4,  6525.0,  6050.0,  5596.1,
     5162.9,  4750.0,  4357.1,  3983.9,  3630.0,  3295.0,  2978.6,  2680.4,
     2400.0,  2137.1,  1891.4,  1662.5,  1450.0,  1253.6,  1072.9,  907.5,
      757.1,   621.4,   500.0,   392.5,   298.6,   217.9,   150.0,   94.6,
       51.4,    20.0,     0.0]
    ----
    Original code:    
    CASE ( 2, 3 )
    ! Height-based hybrid vertical coordinate on input
    ! Vertical grid specified in terms of hhl
    ! here hhl depends only on the zak, zbk and vcflat

    IF     (vc_type%ivctype == 2) THEN
      ! "standard" coordinate with zak, zbk

      ! Calculate the inverse coordinate transformation, i.e. the zak's and zbk's
      vc_type%kflat = 0
      DO k = 1, ke+1
        IF( vc_type%vert_coord(k) >= vc_type%vcflat ) THEN
          zak(k) = vc_type%vert_coord(k)
          zbk(k) = 0.0_ireals
          vc_type%kflat = k
        ELSE
          zak(k) = vc_type%vert_coord(k)
          zbk(k) = (vc_type%vcflat - vc_type%vert_coord(k))/ vc_type%vcflat
        ENDIF
      ENDDO

      IF (lnew_hhl) THEN
        ! Compute the height of the model half-levels
        hhl(:,:,ke+1) = hsurf(:,:) 
        DO  k = 1, ke
          hhl(:,:,k) = zak(k) + zbk(k)*hhl(:,:,ke+1)
        ENDDO
      ENDIF
    ----      
    Returns:
        HHL field as in COSMO (z',rlon,rlat)
    """
    
    # Namelist constants    
    ivctype = 2
    vcflat = 6000.
    zz_top = 9999.

  
    x = np.arange(1,nx+1)
    y = np.arange(1,ny+1)
    _ks = np.arange(0,nz)
    zak = _ks * 0.
    zbk = _ks * 0.

    HHL = np.zeros((nz,nx,ny))
    
    X,Y = np.meshgrid(x,y)
    c = nx/2.+0.5
    R2=(X-c)**2+(Y-c)**2
    hsurf=np.zeros(R2.shape)
    
    
    if surftopo=='gauss':
        #hsurf = Hm*(2**(-(R2)/(a**2)))  #symmetric Gaussian
        if not ay: ay = a
        hsurf = Hm*(2.**(-((X-c)/a)**2-((Y-c)/ay)**2))
    if surftopo=='bell':
        hsurf = Hm/(1+(R2)/(a**2)) 
    if surftopo=='bell1.5':
        hsurf = Hm/(1+(R2)/(a**2))**1.5     
    if surftopo=='cos2':
        hsurf = Hm*np.cos(np.pi/4.*(R2)/(a**2))**2
        hsurf[R2>2*a**2] = 0.
        
    #hsurf[np.where(hsurf < 0.001)] = 0.
    HHL[-1,:,:] = hsurf    

    kflat = 0
    # Inverse coordinate transfromation to obtain zak and zbk
    for k in _ks:
        if vcoordvec[k] >= vcflat:
            zak[k] = vcoordvec[k]
            zbk[k] = 0.
            kflat = k
        else:
            zak[k] = vcoordvec[k]
            zbk[k] = (vcflat -vcoordvec[k])/vcflat

    # Calcualte HHL
    for k in range(nz-1):
        HHL[k,:,:] = zak[k]+zbk[k]*HHL[-1,:,:]

    print kflat
    return HHL  
    


   
    


# Tableau 20 Colors
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
             
# Tableau Color Blind 10
tableau20blind = [(0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
             (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
             (255, 188, 121), (207, 207, 207)]
  
# Rescale to values between 0 and 1 
for i in range(len(tableau20)):  
    r, g, b = tableau20[i]  
    tableau20[i] = (r / 255., g / 255., b / 255.)
    
for i in range(len(tableau20blind)):  
    r, g, b = tableau20blind[i]  
    tableau20blind[i] = (r / 255., g / 255., b / 255.)
# Use with plt.plot(…, color=tableau[0],…)

def f2daytime(filename):
	return int(filename[4:10])%10000 #COSMO output format assumed
	
def f2ddhhmm(filename):
	return int(filename[4:10])
 
 
if __name__ == '__main__':
    nx,ny=300,300
    x = np.arange(1,nx+1)
    y = np.arange(1,ny+1)


        
    X,Y = np.meshgrid(x,y)
    c = nx/2.+0.5
    R=np.sqrt((X-c)**2+(Y-c)**2)
    
    
    hcos2=HHL_creator(Hm=1000,a=20,surftopo='cos2')[-1,:]
    hgauss=HHL_creator(Hm=1000,a=20,surftopo='gauss')[-1,:]
    hbell=HHL_creator(Hm=1000,a=20,surftopo='bell1.5')[-1,:]

    bell = HHL_creator(surftopo='bell')
    gauss = HHL_creator(surftopo='gauss')
    cos2 = HHL_creator(surftopo='cos2')
    f,ax=plt.subplots(1,3)
    cf=ax[0].contour(gauss[-1,:,:],cmap='terrain')
    cf=ax[1].contour(bell[-1,:,:],cmap='terrain')
    cf=ax[2].contour(cos2[-1,:,:],cmap='terrain')
    plt.colorbar(cf)
    #plt.contour(HHL[:,150,::-1].transpose(),cmap=viridis)
    
