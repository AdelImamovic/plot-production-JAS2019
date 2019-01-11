import matplotlib.pylab as plt
import numpy as np
from numpy import ma
from netCDF4 import Dataset
from alt_colormaps import viridis_inv,viridis,plasma,magma,plasma_inv,inferno,inferno_inv,magma_inv
from analy_utils import decorate_ax, add_customcolorbar, get_vertcoord, adjust_spines, HHL_creator
import os,sys
sys.path.append('/home/adeli/scripts/python/project_A/')
from plot_figures_projectA import add_circle


nhalo=3
SCRATCH = '/scratch/snx3000/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
BASE = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/'
BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/'
cols=['blue','green','red','black','brown']

BASE3D = BASE+'ALLVAR_3D/' 
BASETP = BASE+'TOT_PREC/'


SIMRES = 1.
CONV = (2.*np.pi *6370.*SIMRES / 360)

#filenames 
CIRCMEAN = 'circmean_day_d0d5.nc'
ENSMEAN ='ensmean_day_3d_d0d5.nc'
LWPN = 'LWP_ensmean_day_3d_d0d5.nc'

sm='60_homo'

""" Function to convert time_index (in range 0,48) 
    to local time format, e.g. 11:30."""
daytime = lambda tinst: str(tinst/2) + (':00' if tinst%2==0 else ':30')


Lv = 2264705. #J/kg
c_p = 1005.7 #J/kg
g_g = 9.81 # m/s2


#add expnames:
hrange = [250,500,750,1000,1500,3000] #range(250,1001,250)
arange = [5,7,10,14,20,25,30]
nh = len(hrange)
na = len(arange)

expnames=np.zeros((nh,na),dtype=object)

for i in range(nh):
    for j in range(na):
        expnames[i,j] = 'h'+str(hrange[i])+'a'+str(arange[j])



thermodynvar_ls = {'hPBLr':'-',
                   'LFCr':':',
                   'LCLr':'..'}


def radialwindconv(f=False,ax=False,oro='h1000a20',t=[12,37],BASE=BASE3D,
                   plotvarn='u'): 
    if not f:
        f,ax = plt.subplots(1,1)
    ncflname = BASE3D+oro+'/60_homo/'+CIRCMEAN
    print ncflname
    ncfl = Dataset(ncflname)
    
    if oro.endswith('cos2') or oro.endswith('bell'):
        oron = oro[:-5] # clip off suffix
    else:
        oron = oro
    if oro.endswith('000a14'):
        #ah=expn.replace('a','h').split('h')[1:]
        h,a=getha(oro) #int(ah[0]),int(ah[1])
    else:    
        i,j = np.where(expnames==oron)
        
        h,a = hrange[i],arange[j]
    rlim = 128
     
    if plotvarn=='u':
        plotvar = np.mean(ncfl.variables['Urz'][t[0]:t[1],45:50,:rlim],axis=1)
    else:
        plotvar = np.mean(ncfl.variables['Wrz'][t[0]:t[1],30:31,:rlim],axis=1)
    
    tp = ncfl.variables['TOT_PRECr'][t[0]:t[1],:rlim]
    Z = ncfl.variables['Z'][:]
    ncfl.close()
        
    tplevs=[0.1,0.5,1,2,4,8,16]
    #levs = np.linspace(-5,5,9)
    ulevs = [-6.4,-3.2,-1.6,-0.8,-0.4,0.4,0.8,1.6,3.2,6.4]
    wlevs = np.array(ulevs)/4
    
    if plotvarn=='u':
        levs=ulevs
    else:
        levs=wlevs
    
    mycmap = 'bwr' if plotvarn=='u' else 'BrBG'    
    
    cf1=ax.contourf(plotvar,cmap=mycmap,levels=levs,extend='both',alpha=1)

    cf2=ax.contour(tp,cmap=plasma,levels=tplevs,extend='both',linewidths=0.4)
    ax.clabel(cf2, fmt = '%.1f',inline=False,fontsize=7)
    ax.vlines(a,t[0]-12,t[1]-12,'green',linestyle='dashed')          
    ax.set_xlim([0,40])
    #decorate_ax(ax,"",oro+' (V/V$_0$='+'%.1f'%(volume[i,j]/volume[2,3])+')',"",fnsz=5,ypos=1.05)
    #decorate_ax(ax,"",oro,"",fnsz=10,ypos=1.05)
    return f,ax,cf1,cf2


def plot_uwrz_field(f=False,ax=False, oro='h500a7', vm=[-2,2],tinst=18,var='QVrz',ik=0,jk=0):
    """ Plots the radial wind field for the configuration sm and oro at output
    time instant tinst. Local time is tinst/2.
    
    """

    tinst = int(tinst)
    
    # in case f is not passed, f and ax are created
    if not f:
        f,ax=plt.subplots()

    # Name the axes 
    decorate_ax(ax,'',daytime(tinst),'')
    

    # specify full path of the data

    iensmean = 1
        
    if iensmean:
        srcpath = BASE3D + oro + '/' + sm +'/' + CIRCMEAN
        print srcpath
    else: # 1 single ensemble member first day
        srcpath = BASE3D + oro + '/' + sm +'/'
        seedname= filter(lambda x: x.startswith('circmean_seed'), os.listdir(srcpath))[0]
        srcpath +=seedname

    # Open the data in netcdf format and read it only 
    ncfl = Dataset(srcpath,'r')
    
    # Load the fields of interest
    # COSMO fields dimensions are (time, z,rlat,rlon)
    #
    # Field dimensions in the following are (time,z,r)
    #
    # By specifing '::-1' the second timension is inverted
    urz = ncfl.variables['Urz'][tinst,::-1,:] # U component of wind
    wrz = ncfl.variables['Wrz'][tinst,::-1,:] # W component of wind

    speedrz = ncfl.variables['Speedrz'][tinst,::-1,:] # wind speed
    rhorz = ncfl.variables['RHOrz'][:,::-1,:] # density in rz coordinates
    TOTP = ncfl.variables['TOT_PRECr'][:,:]/2. # half-hourly output
    
    LCLr= ncfl.variables['LCL_MLr'][tinst,:]/1000. # LCL
    LFCr = ncfl.variables['LFC_MLr'][tinst,:]/1000. # LFC
    hPBLr = ncfl.variables['HPBLr'][tinst,:]/1000.   
    
    CAPEr = ncfl.variables['CAPE_MLr'][tinst,:] # CAPE
    CINr = ncfl.variables['CIN_MLr'][tinst,:] # CIN
    
    TOTPsum = np.sum(TOTP[:tinst,:],axis=0)
    TOTPsum = TOTP[tinst,:]
    
    TPthres = 0.5
    rainyr=np.where(TOTPsum>TPthres)


    if var: #select a specific humidity var
        qvrz = ncfl.variables[var][:,::-1,:] # specficif humidity
        qcrz = ncfl.variables['QCrz'][:,::-1,:]
    else: # otherwise sum it all up:
        qxrzvars = filter(lambda x: x.startswith('Q'),ncfl.variables.keys())
        
        qvrz = np.zeros(ncfl.variables['QVrz'].shape)
        for var in qxrzvars:
            print var
            qvrz+=ncfl.variables[var][:,::-1,:]
    
    tref = 12 # 06:00
    ipltdiff= 0
    if ipltdiff:
        plotvar = rhorz[tinst,:]*qvrz[tinst,:]-rhorz[tref,:]*qvrz[tref,:]
        qcrz = rhorz[tinst,:]*qcrz[tinst,:]-rhorz[tref,:]*qcrz[tref,:]
    else:
        plotvar = rhorz[tinst,:]*qvrz[tinst,:]
        

    rhorz = rhorz[tinst,:]
    
    
    # Load the grid in r, z coordinate
    X = ncfl.variables['X'][:]
    
     
    # nasty and temporary workaround as Z is empty in single seed data - just
    # take it from the circmean field
    srcpath2 = BASE3D + oro + '/' + sm +'/' + CIRCMEAN
    ncfl2=Dataset(srcpath2)
    Z = ncfl2.variables['Z'][::-1,:]/1000. # highest index <-> ground
    ncfl2.close()    
    
    # Add mountain contour
    r=X[0,:]
    ax.plot(r,Z[0,:],color='grey',linewidth=2)
    
#    ax.plot(r,hPBLr,'k-.',linewidth=1)
#    ax.plot(r,LCLr,'k-',linewidth=1)    
    
    ax.fill_between(X[0,:],Z[0,:],0,color='lightgrey')
    
    #indicate rainy region
    
    ifillrainy = False
    if ifillrainy:
        ax.fill_between(X[0,rainyr][0],Z[0,rainyr][0],color='lightblue',alpha=0.8) 
        ax.plot(X[0,rainyr][0],Z[0,rainyr][0],color='lightblue',linewidth=2)
    else:
        pass#ax.plot(X[0,rainyr][0],Z[0,rainyr][0],color='blue',linewidth=2,linestyle='dashed')
    
    # Velocity vectors with speed (length) smaller than speedthres are not 
    # shown
    oroorig = oro
    if oro.endswith('bell') or oro.endswith('cos2'):
        
        oro = oro[:-5]
        print oro
        
    i,j=np.where(expnames==oro)

    if i==[]:
        print 'HH'
        
    oro = oro.split('_')[0]
    h,a = getha(oro)
    #h,a = hrange[i],arange[j]
    print h,a
    xrang = [0,40]
    yrang = [0,10.]
    delx = xrang[1]-xrang[0]
    dely = yrang[1]-yrang[0]
    
    speedthres = 0.1# m/s
    M = speedrz[:] < speedthres
    U = ma.masked_array(urz[:,:],mask=M)
    scalefac_w = 1. #delx/dely
    W = ma.masked_array(scalefac_w*wrz[:,:],mask=M) 
    
    
    
    if 0: # REMOVE
        fact_agl = -dely/delx
        sin = np.arcsin(fact_agl)
        cos = np.arccos(fact_agl)
        U=cos*U-sin*W
        W=sin*U+cos*W
    # Plot radial velocity
    urzlevs = np.linspace(vm[0],vm[1],10)


    levs = np.array([2,4,6,8,10,12,14,16])
    #levs = np.array([0.25,0.5,1,2,4,8,16,32])
    #levs = np.array([0.5,1,2,4,8,16,32])
    if ipltdiff:
        levs=levs/8. #np.linspace(0.001,0.004,11)d
    #else:
    #    levs=np.array([0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004])*50
    #if var.startswith('QC'): levs=np.linspace(0.1,1.5,9)/1000.
    
    #cf = ax.contourf(X,Z,plotvar*1000., levels=levs*1000.,cmap=viridis_inv,extend='both')
    
    ipltmoistening = False    
    if ipltmoistening:
        cf = ax.contourf(X,Z,plotvar*1000.,cmap=viridis_inv,extend='both',
                    levels=levs)
    else:
        cf = 0
                    
    #cfc = ax.contourf(X,Z,qvrz*1000.,cmap='Blues',extend='both', alpha=1,
    #             levels=levs/4.)              

    iplt_thermodynvar = False
    if iplt_thermodynvar:
        for var in [LFCr]:
            print var.shape
            print X.shape
            ax.plot(var)
            
        
    qstx = 2
    qstz = 2
    Q=ax.quiver(X[::qstz,::qstx],(Z+0.02)[::qstz,::qstx],U[::qstz,::qstx],W[::qstz,::qstx],
                units='inches',scale=10,angles='xy', color='red',pivot='mid',
                headwidth=4,width=0.009)
    #Q=ax.quiverkey(Q,a/2,h/2000.,1,r'$1\frac{m}{s}$',labelpos='N',coordinates='data')

    if jk==2 or oroorig.endswith('cos2'):#(ik,jk)==(1,0):
        ax.quiverkey(Q,-2,0.05,1,r'$1 \frac{m}{s}$',labelsep=0.05,
                     labelpos='N',coordinates='data')


    ax.set_xticks(np.arange(0,80,20)) 

    ax.set_xlim(xrang)
    ax.set_ylim(yrang)

    if ik==2:
        ax.set_xlabel('r in km')
    if jk==0:
        ax.set_ylabel('z in km')    
    
    ax.grid(True)

    ncfl.close()
    
    return f,ax ,cf # plot_uwrz_field


sq2 = np.sqrt(2)

hrange = [250,500,750,1000,1500] # range(250,1001,250)
arange = [5,7,10,14,20,25,30] 
arangeval = [5,5*sq2,10,10.*sq2,20,25,30]
nh = len(hrange)
na = len(arange)

expnames=np.zeros((nh,na),dtype=object)
volume=np.zeros((nh,na))

slope=np.zeros((nh,na))

for i in range(nh):
    for j in range(na):
        volume[i,j] = hrange[i]*(arangeval[j]**2)
        slope[i,j] = hrange[i]/arangeval[j]

for i in range(nh):
    for j in range(na):
        expnames[i,j] = 'h'+str(hrange[i])+'a'+str(arange[j])

nx = 512
ny = 512
x = np.arange(1,nx+1)
y = np.arange(1,ny+1)
X,Y = np.meshgrid(x,y)
c = 256
R = np.sqrt((X-c)**2+(Y-c)**2)


nxs = 512/2
nys = 512/2
xs= np.arange(1,nxs+1)
ys = np.arange(1,nys+1)
Xs,Ys = np.meshgrid(xs,ys)
cs = 256/2
Rs = np.sqrt((Xs-cs)**2+(Ys-cs)**2)

colorsbyh = np.copy(expnames)
colnamesbyh = ['black','blue','green','red','violet']

for j in range(na):
    colorsbyh[:,j] = np.array(colnamesbyh)

colorsbya = np.copy(expnames)
colnamesbya = ['grey','brown','orange','sienna','lightblue','darkblue','black']

for i in range(nh):
    colorsbya[i,:] = np.array(colnamesbya)


for e in expnames:
    for h in e:
        h = '1'
#mask = R < 20.  
        
def plt_LWP_E_contplt(f=False,ax=False,expn='h3000a14', inorm=0):#,i=3,j=2):
    if not f:
        f,ax=plt.subplots(1,1)

    #expn = expnames[i,j]


        
    
    fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/circmean_LWP_ensmean_day_3d_d0d5.nc'
    nc = Dataset(fln)
    TWPr = nc.variables['TWPr'][12:37,0:41]
    EFLUXr = nc.variables['EFLUXr'][12:37,0:41]
    
    
    if inorm:
        #TWPrsat = nc.variables['TWPrsat']
        # calculate saturation water volume Wtot_sat
        a = 17.67
        b = 243.5       
        c = 6.112
        esat = lambda T: c*np.exp(a*(T-273.15)/(b+(T-273.15)))    
        fln2 = BASE+'/ALLVAR_3D/'+expn+'/60_homo/circmean_day_d0d5.nc'
        nc2 = Dataset(fln2)
        
        Trz = nc2.variables['Trz'][12:37,:,0:41]
        Prz = nc2.variables['Prz'][12:37,:,0:41]
        RHOrz = nc2.variables['RHOrz'][12:37,:, 0:41]
        QVrz = nc2.variables['QVrz'][12:37,:, 0:41]
        Z = nc2.variables['Z'][:,0:41]
        h,a = getha(expn)
        Zgrid = HHL_creator(Hm=h,nx=512,ny=512,nz=51,a=a,ay=False,surftopo='gauss')
        
                
        
        Ztemp = Zgrid[:,256,256:256+41]    
        deltaz = Ztemp[:-1,:]-Ztemp[1:,:]   
        
        
               
        X = nc2.variables['X'][:,0:41]
        qvsat = 0.622*esat(Trz)/Prz*1000 # converg to g/kg
        
        print QVrz[12,-1,:]/(qvsat[12,-1,:]), deltaz[-1,:]
        
        
        WTOTqvsat = np.zeros([25,41])
      
        for var in (WTOTqvsat,qvsat,RHOrz,deltaz): 
            print '0' #var.shape
            

        for ti in range(12,37,1):
            # f,ax2=plt.subplots()
            # ax2.contourf(qvsat[i-12,:]*RHOrz[i-12,:]*deltaz,cmap=viridis)
#            f,ax = plt.subplots(3,1)
#            ax[0].contourf(qvsat[i-12,:])
#            ax[1].contourf(RHOrz[i-12,:])
#            ax[2].contourf(deltaz)
#            
#            continue
            WTOTqvsat[i-12,:] = np.sum(qvsat[i-12,:]*RHOrz[i-12,:]*deltaz,axis=0)
        
        #WTOTtemp = np.sum(WTOTqvsat,axis=1)
        f,ax=plt.subplots(1,2)
        levs = np.arange(0,0.2,0.01)
        cf = ax[0].contourf(qvsat[:,-1,:],cmap=viridis,levels=levs,extend='both')
        
        cf = ax[1].contourf(QVrz[:,-1,:],
                            cmap=viridis,extend='both')
        plt.colorbar(cf)
         
        
        nc2.close()
    
        assert 0    
    
    levs = np.arange(10,40,2)
    a = int(expn.split('a')[-1])
    
    cf=ax.contourf(TWPr,levels=levs,cmap=inferno_inv,extend='both',)
    cf2=ax.contour(EFLUXr*1000*1800/Lv,cmap='Blues',extend='both')
    ax.vlines(a,0,24,color='green',linestyle='dashed')

    ttcks = np.arange(0,25,6)
    ttckl = map(daytime,ttcks+12)
    
    rtcks = np.arange(0,41,10)
    ax.set_yticks(ttcks)
    ax.set_yticklabels(ttckl)
    
    ax.set_xticks(rtcks)
    
    ax.clabel(cf2, fmt = '%.1f',inline=False,fontsize=7)
    ax.grid(True)
        
    nc.close() 
    return cf,cf2       
        
def plot_theta(f=False,ax=False,oro='h500a14',t=20,ipltsym=1):
    """print docstring."""
    if not f:
        f,ax=plt.subplots()
        
    sm='60_homo'
        
    srcpath = BASE3D + oro + '/' + sm +'/' + CIRCMEAN #'circmean_seed_d0.nc'
    print srcpath

    data = Dataset(srcpath)
    
    RL = 287. #J/kg/K
    cp = 1005. #J/kg/K

    g = RL/cp
    
    T = data.variables['Trz'][t,::-1,:]
    p = data.variables['Prz'][t,::-1,:]
    rho = data.variables['RHOrz'][t,::-1,:]
    qc = data.variables['QCrz'][t,::-1,:]

    qwc = rho*qc*1000.
    print np.max(qwc)
    qwc[qwc<0.0047] = 0.#True

    p0 = p[0,-1]


    X = data.variables['X'][:]
    
     
    # nasty and temporary workaround as Z is empty in single seed data - just
    # take it from the circmean field
    srcpath2 = BASE3D + oro + '/' + sm +'/' + CIRCMEAN
    data2=Dataset(srcpath2)
    Z = data2.variables['Z'][::-1,:]/1000. # highest index <-> ground
    

    # Add mountain contour

    if not ipltsym:
        r=X[0,:]
        ax.plot(r,Z[0,:],color='grey',linewidth=2)
        
    #    ax.plot(r,hPBLr,'k-.',linewidth=1)
    #    ax.plot(r,LCLr,'k-',linewidth=1)    
        
        ax.fill_between(X[0,:],Z[0,:],0,color='lightgrey')
   

    i,j = np.where(expnames==oro)
    #h,a=hrange[i],arange[j]
    a=14.
    
    Theta = T*(p0/p)**g
    
    xloglow=np.log10(0.01)
    xloghi=np.log10(0.2)
    qc_levs = np.logspace(xloglow,xloghi,8)*1000
    qc_levs  =[0.01,0.05,0.1]
    #qwc*=1000

    levs=np.arange(293.,300.1,0.5)    
    plevs = np.linspace(95000,75000,10)
    

    if ipltsym:
        #copy and mirror

        Thetasym = np.hstack((Theta[:,::-1],Theta))
        qwcsym = np.hstack((qwc[:,::-1],qwc))
        Xsym = np.hstack((X[:,::-1]*-1.,X))
        Zsym = np.hstack((Z[:,::-1],Z))

        cf=ax.contour(Xsym,Zsym,Thetasym,cmap=plasma,levels=levs,extend='both')

        ax.contour(Xsym,Zsym,qwcsym,cmap='Blues',levels=qc_levs)
        #cb=plt.colorbar(cf)

        tfac = 5
        xtcks = np.arange(-tfac*a,tfac*a,a)
        xtcklabs = map(lambda x: str(int(x/a)),xtcks)
        ax.set_xticks(xtcks)
        ax.set_xticklabels(xtcklabs)
        r=Xsym[0,:]
        ax.plot(r,Zsym[0,:],color='grey',linewidth=2)
        ax.set_xlim([-60,60.])
        ax.set_xlabel('x in a')
        ax.fill_between(Xsym[0,:],Zsym[0,:],0,color='lightgrey')

    else:
        cf=ax.contourf(X,Z,Theta,cmap=plasma_inv,linewidths=0,levels=levs,extend='both',manual=True)
        cfc=ax.contour(X,Z,qwc,cmap='Greys',levels=qc_levs,extend='both',linewidths=2.5)
        #ax.contour(X,Z,p,cmap='Blues',levels=plevs)
        ax.set_xlim([0,40.])
    #ax.clabel(cfc, fmt = '%.1f',inline=False,fontsize=8,color='k')
    ax.set_ylim([0,2.5])
    
    
    
    ax.grid(True)
    data.close()
    data2.close()
        
        
    return f,ax,cf
    
    


def u_csect(f=False,ax=False,oro='h1500a14',BASE=BASEnow,t=20):
    if not f:
        f,ax=plt.subplots()

    sm='60_homo'
    
    srcpath = BASE+'ALLVAR_3D/' + oro + '/' + sm +'/ensmean_day_3d_d0d5.nc'
    print srcpath

    data = Dataset(srcpath)
    
    RL = 287. #J/kg/K
    cp = 1005. #J/kg/K

    g = RL/cp
    
    T = data.variables['T'][t,::-1,nhalo:-nhalo,nhalo:-nhalo]
    p = data.variables['P'][t,::-1,nhalo:-nhalo,nhalo:-nhalo]
    rho = p/T
    u = data.variables['U'][t,::-1,nhalo:-nhalo,nhalo:-nhalo]
    qc = data.variables['QC'][t,::-1,256,nhalo:-nhalo]
    qi = data.variables['QI'][t,::-1,256,nhalo:-nhalo]
    qg = data.variables['QG'][t,::-1,256,nhalo:-nhalo]
    qs = data.variables['QS'][t,::-1,256,nhalo:-nhalo]
    
    qt = (qc + qi+qg+qs)*rho[:,256,:]
        
    
    w = data.variables['W'][t,::-1,nhalo:-nhalo,nhalo:-nhalo]
    x = data.variables['rlon'][nhalo:-nhalo]*CONV-(256+2.)
    z = np.arange(50)    
    
    X,Ztemp = np.meshgrid(x,z)    
    
    hacode = {'h500a14':(500,14),'h1500a14':(1500,14),'h3000a14':(3000,14)}
    h,a = hacode[oro]
    
    ulevs = np.arange(-4.4,5.,0.8)
    p0=95000. #Pa        
    Theta = T*(p0/p)**g
    from analy_utils import HHL_creator
    HHL = HHL_creator(h,nx=512,ny=512,a=a,surftopo='gauss')[::-1,256,:]/1000.
    
    iplt_vec = 0
    if iplt_vec:
        xt = X 
        ht = HHL[1:,:]
        ut = u[:,256,:]
        wt = w[1:,256,:]
        
        sp = ut**2+wt**2
        masksp = sp < 0.01
        ut = ma.masked_array(ut,mask=masksp)
        wt = ma.masked_array(wt,mask=masksp)
      
        qmi, qma = np.min(qt),np.max(qt)
        print qmi,qma
        qlevs = np.logspace(qma/1000.,qma,10)
        qlevs = np.log10(qlevs)
        ax.contourf(xt,ht,qt,levels=qlevs,extend='both',cmap=viridis_inv)
        ax.quiver(xt,ht,ut,wt,angles='xy')
        
    

    Thetacross = Theta[:20,256,:]
    
    l = 40
    u = u[:l,256,:]
    uabs = np.copy(u)
    w = w[:l,256,:]
    HHL = HHL[:l,:]
    X = X[:l,:]
    print HHL.shape, Thetacross.shape, X.shape
    
    ipltuanom=1
    if ipltuanom:
        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
        BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/'
        BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/'
        ucode = {BASEnow:0.,BASEw05:1.,BASEw1:2.,BASEw2:4.,BASEw4:8.}
        usurf = ucode[BASE]
        
        print usurf
        anom = map(lambda y: y*(4*usurf)/11.5+usurf,HHL+0.02)
# Test on correction of anomaly  
#        print HHL
#        cf=plt.contourf(X,HHL,anom,cmap=viridis)
#        plt.colorbar(cf)        
#        assert 0
        u-=anom
        #ulevs = np.arange(-4.25,4.261,0.5)
        ulevs = np.arange(-3.5,3.51,1)
    

        
    
    cf = ax.contourf(X,HHL,u,levels=ulevs,cmap='bwr',extend='both')
    
    ax.plot(X[0,:],HHL[0,:],color='grey')
    ax.fill_between(X[0,:],HHL[0,:],0,color='lightgrey')
    
    
    iplt_wf = 0    
    if iplt_wf:
       s = 4
       sp = u**2+w**2
       masksp = sp < 0.01
       um = ma.masked_array(uabs,mask=masksp)
       ax.quiver(X[:,::s],HHL[:,::s]+0.02,um[:,::s],w[:,::s],angles='xy',
                 color='black',pivot='mid',alpha=0.8,headwidth=3)     
    
    ax.set_xlim([-20,8*a])
    ax.set_ylim([0.,12.])
    ax.grid(True)
    data.close()
    return f,ax,cf

def plot_radflux(oro='h500a20',ipltdry=0):
    h,a=getha(oro)  
    ri=30.
    h_atr = h*(2**(-(ri/a)**2))

    from analy_utils import get_vertcoord
    
    vcoord=get_vertcoord(hsurf=h_atr)[::-1]
    deltaz = vcoord[1:]-vcoord[:-1]
    ztot = np.sum(deltaz)
    zcoord = 0.5*(vcoord[1:]+vcoord[:-1])
    #zcoord = zcoord[::-1]
    #deltaz = deltaz[::-1]

    srcpath = BASE +'ALLVAR_3D/' + oro + '/60_homo/circmean_day_d0d5.nc'
    data = Dataset(srcpath)
    
    urz = data.variables['Urz'][:,:,ri]
    rhorz = data.variables['RHOrz'][:,:,ri]
    qvrz = data.variables['QVrz'][:,:,ri]
    data.close()
    
    
    srcpath = BASE +'ALLVAR_3D/' + oro + '/60_homo/ensmean_day_3d_d0d5.nc'

    ensmdata = Dataset(srcpath)  
    
    analyregion = R<30.
    Es = ensmdata.variables['EFLUX'][:,-1,nhalo:-nhalo,nhalo:-nhalo]
    Hs = ensmdata.variables['HFLUX'][:,-1,nhalo:-nhalo,nhalo:-nhalo]
    ensmdata.close()
    
    #Es = Es[:,mask]
    #Hs = Hs[:,mask]
    
    Esmean = np.mean(Es[:,analyregion],axis=1)
    Hsmean = np.mean(Hs[:,analyregion],axis=1)
    
    HEtot = np.mean(Hsmean[12:36]),np.mean(Esmean[12:36])
    
    if ipltdry:
        qvrz=1.
    flux = urz*rhorz*qvrz
    intflux = np.copy(flux)
    
    #weight with cell height:
    for i in range(48):
        intflux[i,:]*=(deltaz/ztot)
    
    #cf=plt.contourf(intflux,cmap='bwr',vmin=-0.003,vmax=0.003)
    
    flxevol=np.sum(intflux,axis=1)
    
    
    
    convtimes=np.where(flxevol[:28]<0.)
    
    convflux =flxevol[convtimes]    
    
    totalflux = -np.sum(convflux)*1800#.*2*np.pi*(ri*1000)*ztot/1e9

    
    return flxevol,totalflux,(Hsmean,Esmean), HEtot


def getha(oroname):
    ah=oroname.replace('a','h').split('h')[1:]
    h,a=int(ah[0]),int(ah[1])
    return h,a 

def calc_MSE(f=False,ax=False,oro='h500a14', icalc=False):
    c = 256
    wi = 30
    iplt = 0
    if not f:
        f,ax=plt.subplots(1,1)
    if icalc:
        BASE= '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        datapath = BASE+'ALLVAR_3D/'+oro+'/60_homo/ensmean_day_3d_d0d5.nc'
    
        data=Dataset(datapath)
    

    

#    if 0:
#        MSEint = data.variables['MSEint'][:]   
#        MSEmtn_ts = MSEint[:] 
#        MSEfar = 0
#        return MSEint
#    else:
#    data.createVariable('MSE',datatype=float,dimensions=('time','level','rlat'))
#    data.createVariable('MSEint',datatype=float,dimensions=('time'))

        from analy_utils import HHL_creator
        h,a = getha(oro)
        Zfull = HHL_creator(h,nx=512,ny=512,a=a,surftopo='gauss')[:,256]
        Z = 0.5*(Zfull[0:-1]+Zfull[1:])
        deltaz = Zfull[1:]-Zfull[0:-1]
        x=np.arange(-256,256)
        z=np.arange(50)
        X,b = np.meshgrid(x,z)
    
    
        Zdraw = Z[::-1,:]
        
            
        trange = np.arange(0,48,2) 
        MSEmtn_ts = np.zeros(trange.shape)
        MSEfar_ts = np.zeros(trange.shape)
    
        
        if iplt:
            nstri=nhalo
        else:
            nstri = 256-30+nhalo
        deltaz = deltaz#[:,nstri-nhalo:-nstri+nhalo]
            
        
        # cross section through mountain
    
        r = np.arange(-0.256,0.256,0.001)
        
        Arings = np.pi*r**2
        
    
        Rd = 287.
        
    
        evallev = 1 #corresponds to 6287.5 m, i.e. above the height where zlines are equal
    
        
        
        Zz=Z[evallev:,:]#:,226:-226]
        deltaz/=1000.
        deltaz=deltaz[evallev:,:]
        for ti in trange:
            Ta = data.variables['T'][ti,evallev:,257,nhalo:-nhalo]
            QVa = data.variables['QV'][ti,evallev:,257,nhalo:-nhalo]
            Pa = data.variables['P'][ti,evallev:,257,nhalo:-nhalo]
            RHOa = Pa/(Ta*Rd)
     
            
            iDSE=0
            if iDSE:
                MSE = (c_p*Ta+g_g*Zz)/c_p
            else:            
                MSE = (c_p*Ta+Lv*QVa+g_g*Zz)/c_p
            
        
        
            MSEint = np.mean(MSE*RHOa*deltaz*Arings,axis=0)/np.mean(deltaz*RHOa*Arings,axis=0)      



                
        
        
#        plt.contourf(MSE[::-1,c-wi:c+wi],cmap=magma,levels=np.linspace(290,320,10))
#        
#        assert 0 
        
            MSEmtn = np.mean(MSEint[c-wi:c+wi])
            MSEfar = np.mean(MSEint[c+4*wi:])
    
            MSEmtn_ts[ti/2] = MSEmtn
            MSEfar_ts[ti/2] = MSEfar
            data.close()
            return MSEmtn_ts,MSEfar_ts
    else:
        BASE= '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        datapath = BASE+'ALLVAR_3D/'+oro+'/60_homo/interp_ensmean_day_3d_d0d5.nc'           
        data = Dataset(datapath)
        it = 24
        MSE = data.variables['DSE'][:,::-1,nhalo:-nhalo,nhalo:-nhalo]
        h,a=getha(oro)
        if h>1000:
            MSE/=c_p
        f,ax=plt.subplots(1,1)
        
        Z = data.variables['grid_int'][::-1,nhalo:-nhalo,nhalo:-nhalo]
        X,a = np.meshgrid(np.arange(512),np.arange(29))
        levs = np.arange(293,300,0.5)
        
        icalcbulk = 1
        icalcgrad = 0
        if icalcbulk:        
            mask = np.logical_and(R>25,R<30.)
            maskout =np.logical_and(R>30.,R<33.)
            
            
            
            kz = 21
            # bad fix 
            mtnpoints = MSE[0,:]<100.
            MSE[:,mtnpoints] = np.nan
    
            ax.set_title(str(Z[kz,0,0]))
    
            tempin = MSE[:,:kz,mask]
            tempout = MSE[:,:kz,maskout]
    
            MSEmean = np.nanmean(tempin,axis=(1,2))
            MSEmean_out = np.nanmean(tempout,axis=(1,2))
            ax.plot(MSEmean-MSEmean_out,label=oro)
            
        if icalcgrad:
            mask = np.logical_and(R<30.,R>28)
            maskout = np.logical_and(R>30.,R<32)
            MSEtempin = np.mean(MSE[:,:,mask],axis=2)
            MSEtempout = np.mean(MSE[:,:,maskout],axis=2)
            
            for ti in range(12,37,6):
                #plt.plot(MSE[ti,:,256,285]-MSE[ti,:,256,286],Z[:,256,256],label=daytime(ti))
                plt.plot(MSEtempin[ti,:]-MSEtempout[ti,:],Z[:,256,256],label=daytime(ti))
            plt.legend()
        
        data.close()
    
    if iplt:
        f,ax=plt.subplots(1,1)  
        levs = np.linspace(310,330,11)
        cf=ax.contourf(X,Zdraw,MSE[::-1,:]/c_p,cmap=magma,levsels=levs,extend='both')
        ax.set_ylim([0,10000])
        ax.set_title(daytime(ti))
        add_customcolorbar(f,cf,'$\Theta_e$ in K',pos='right')
        ax.set_xlim([-40,40])
        f.savefig('/home/adeli/figures/project_B2/figoutline/MSEMov/'+str(ti)+'.png',dpi=100)

    pass
    
    


def calc_P(BASE,oro):
    wi =30.
    datapath = BASE+'ALLVAR_3D/'+oro+'/60_homo/ensmean_day_3d_d0d5.nc'
    data = Dataset(datapath)
    mask = R<wi
    TP = data.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo]
    idbg=0
    if idbg:
        plt.contourf(np.sum(TP,axis=0))
        assert 0
    TP = np.mean(TP[:,mask],axis=1)
    P = np.sum(TP)
    data.close()
    return P/2.
    
if __name__=="__main__":
    FIGPATH = '/home/adeli/figures/project_B2/figoutline/submission_2/'
    mydpi = 300
    isavefig = 0
    
    
    BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
    BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
    BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
    BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/'
    BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/'
  
  

    i_nodc_plot = 0
    if i_nodc_plot:
        srcpath = BASEnow+'ALLVAR_3D/h1000a10_nodc/60_homo/'+CIRCMEAN
        dset = Dataset(srcpath)

        urz = dset.variables['Urz']
        f,ax=plt.subplots(1,2)
        axc,axe=ax
        axc.contourf(urz[:,-1,:],cmap='bwr',levels=np.arange(-3.8,3.9,0.4),extend='both')
        axc.set_xlim([0,30])
        
        
        axc.set_yticks(np.arange(0,240,20))
        axc.set_yticklabels(map(daytime,np.arange(0,240,20)))
        axc.grid(True)
        axc.set_xlabel('distance from mtn top in km')
        axc.set_ylabel('time in h')        
        
        axe.set_xlabel('time in h')
        axe.set_ylabel('radial speed at r=5 km in m/s')
        axe.set_xticks(np.arange(0,240,20))
        axe.set_xticklabels(map(daytime,np.arange(0,240,20)))
        axe.grid(True)
        axe.plot(urz[:96,-1,5])
        figname=FIGPATH+'nodc_HOVM.pdf'
        f.savefig(figname)
        dset.close()        


    if 0:
        f,ax = plt.subplots(1,2,sharex=True,sharey=True,figsize=(10,5))
        plt_LWP_E_contplt(f,ax[0],expn='h1000a14',inorm=0)
        cf,cfE=plt_LWP_E_contplt(f,ax[1],expn='h4000a14',inorm=0)
        add_customcolorbar(f,cf,'$W_{tot}$ in mm',orient='vertical',pos='right')
        for a in ax:
            a.set_xlabel('r in km')
        ax[0].set_ylabel('time')
        figname = FIGPATH+'wtot_h1000h3000.pdf'
        f.savefig(figname,bbox_inches='tight')
        assert 0
    
#    calc_P(BASEnow,'h250a20')
#    
#    assert 0
    iplt_MSE = 0
    if iplt_MSE:
        f,ax = plt.subplots()
        oros = ['h250a14','h750a14','h1000a14','h3000a14']
        for oro in oros:
            calc_MSE(f,ax,oro)
            ax.legend()
            assert 0
        #orons = ['h250a20','h500a14','h1000a10'] 
        
        orons = ['h250a20','h500a14','h1000a10']#,'h500a20','h1000a20','h1500a20']
        
        cols = ['r','b','g','k']    
        
        i=0
        f,axmtn=plt.subplots(1,1)#,sharey=True)
        #axmtn,axfar = ax
    
        test = expnames.reshape((35,))    
        
        resultsMSEmtn = {x: np.zeros(24) for x in test}
        resultsMSEfar = {x: np.zeros(24) for x in test}
        
        i=0
        for oro in test:
            print oro,i
            MSEmtn,MSEfar=calc_MSE(oro)
            
            resultsMSEmtn[oro] = MSEmtn 
            resultsMSEfar[oro] = MSEfar
            i+=1
            
        assert 0
        
        
        iplt_pic_titlepage = 0
        if iplt_pic_titlepage:
            u_csect(oro='h1500a14')
            





    iplt_VARh_radflux = 0
    
    if iplt_VARh_radflux:
        myhrange = np.array([250,500,750,1000,1500,2000,3000,4000])
        myexpnames = ['h'+str(hi)+'a14' for hi in myhrange]
        
        nhs = len(myexpnames)
        tf = np.zeros([nhs])
        tp = np.zeros([nhs])
        tfdry = np.zeros([nhs])
        flxevols = np.zeros((nhs,48))
        HEsurfflxevols = np.zeros((2,nhs,48)) #EH, h, time
        Htf = np.copy(tf)
        Etf = np.copy(tf)
        for i in range(nhs):
                oro = myexpnames[i]
                print oro,i
                flxevols[i,:],tf[i],HE,HEmn = plot_radflux(oro=oro,ipltdry=0)
                temp,tfdry[i],temp2,temp3 = plot_radflux(oro=oro,ipltdry=1)
                HEsurfflxevols[0,i,:],HEsurfflxevols[1,i,:] = HE[0],HE[1]
                
                tp[i] = calc_P(BASEnow,oro)
                Htf[i] = HEmn[0]
                Etf[i] = HEmn[1]
        flxevols*=1000.
        area = np.pi*(30000.)**2
        tf*=1e9/area    
        
        f,ax=plt.subplots(1,2,figsize=(7,3))
        axradfl,axtftp = ax
        
        #for i in range(nhs):
        #    axradfl.plot(-flxevols[i,:])
        
        vols = myhrange/500.
        #axradfl.set_ylim([0,60])
        axradfl.plot(vols,tf,marker='o',color='green')
        
        ipltEH = False
        if ipltEH:
            axEH.plot(vols,Etf,marker='s',color='blue')
            axEH.plot(vols,Htf,marker='v',color='red')
            axEH.set_ylabel('surface heat flux in W m$^{-2}$')
            axEH.set_xlabel('$V/V_0$')
        axradfl.set_xlabel('$V/V_0$')
       
        
        axradfl.set_ylabel('$F_{v,tot}$ in g m$^{-2}$ s$^{-1}$')
        
        
        #axtftp.scatter(tfdry/tfdry[1],tp/tp[1],color='brown')
        axtftp.scatter(tf/tf[1],tp/tp[1],color='k')
        
        axtftp.set_xlabel('$F_{v,tot}/F_{v,tot}(h500a14)$')
        axtftp.set_ylabel('$P/P_0$')
        
        
        axtftp.plot(np.linspace(0,4.1,10),np.linspace(0,4.1,10),'k')
 
        
        f.subplots_adjust(wspace=+0.5)
        for a in ax: 
            a.grid()
            adjust_spines(a,['left','bottom'])
        for a in ax[:1]:
            atop = a.twiny()
            atop.set_xticks(vols)
            atop.set_xticklabels(myhrange, fontsize=8,rotation=90)#rotation=90)
            atop.set_xlabel('h in m')
        axtftp.set_xticks(range(5))    
        axtftp.set_yticks(range(5))
        axtftp_top = axtftp.twiny()
        axtftp_top.set_xticks(tf/tf[1]/4.)
        axtftp_top.set_xticklabels(myhrange,fontsize=8,rotation=90)
        axtftp_top.set_xlabel('h in m')
        axtftp.set_xlim(([0,4.]))
        axtftp.set_ylim(([0,4.]))
        figname = FIGPATH +'expl.pdf'
        f.savefig(figname,dpi=300,bbox_inches='tight')
        assert 0
    
    iplt_radflux = 1
    if iplt_radflux:
        oro = 'h500a20'
        tf = np.zeros(expnames.shape)
        pt = np.zeros(expnames.shape)
        
        nh,na = expnames.shape
        flxevols = np.zeros((nh,na,48))
        HEsurfflxevols = np.zeros((2,5,7,48)) #EH, h, a, time
        Htf = np.copy(tf)
        Etf = np.copy(tf)         
        
        for i in range(5):
            for j in range(7):
                print oro
                oro = expnames[i,j]
                flxevols[i,j,:],tf[i,j],HE,HEmn = plot_radflux(oro=oro)
                HEsurfflxevols[0,i,j,:],HEsurfflxevols[1,i,j,:] = HE[0],HE[1]
                Htf[i,j] = HEmn[0]
                Etf[i,j] = HEmn[1]
                pt[i,j] = calc_P(BASEnow,oro)
                #scale to mm        
        flxevols*=1000.
        area = np.pi*(30000.)**2
        tf*=1e9/area
        
        
        # plot results:
        f,ax2=plt.subplots(2,2,figsize=(10,10))
        ax1,ax=ax2[1,:]
        ax1EH,axEH=ax2[0,:]
        vscal = volume[3,2]
        
        fp,axp=plt.subplots(1,1,figsize=(4,4))
        
#        for i in range(nh):
#           for j in range(na):
#                axp.text(tf[i,j]/tf[1,3],pt[i,j]/pt[1,3],expnames[i,j])
                
        axp.scatter(tf/tf[1,3],pt/pt[1,3],color='k')
        axp.plot(np.arange(0,8,1),np.arange(0,8,1),color='k')
        axp.grid(True)
        adjust_spines(axp,['left','bottom'])
        axp.set_xlim([0,7])
        axp.set_ylim([0,7])


        axp.set_xlabel('$F_{v,tot}/F_{v,tot}(h500a14)$')
        axp.set_ylabel('$P/P_0$')
                
        figname = FIGPATH +'totfl_vs_p.pdf'
        fp.savefig(figname,bbox_inches='tight')
        #  
        for i in range(5):
            
            ax.scatter(volume[i,:]/vscal,tf[i,:],color=colorsbyh[i,0])
            axEH.scatter(volume[i,:]/vscal,Etf[i,:],marker='s',color=colorsbyh[i,0])
            axEH.scatter(volume[i,:]/vscal,Htf[i,:],marker='v',color=colorsbyh[i,0])
       
          
        for i in range(5):
            leglabel = 'h=' +str(hrange[i])+ 'm'
            ax.scatter(volume[i,0]/vscal,tf[i,0],marker='o',color=colorsbyh[i,0],label=leglabel)
        ax.legend(loc='lower right',frameon=False)
        ax.set_xlim([0,14])
        ax.set_ylim([0,100])
        ax.grid(True)
                
                

        
        # plot lateral flux evolution of selected mountains
        ax1.plot(-flxevols[0,4,:],'o-',label=expnames[0,4],markeredgecolor=None,markersize=3,color='k')
        ax1.plot(-flxevols[1,3,:],'o-',label=expnames[1,3],markeredgecolor=None,markersize=3,color='b')
        ax1.plot(-flxevols[3,2,:],'o-',label=expnames[3,2],markeredgecolor=None,markersize=3,color='r')
        
        ax1EH.plot(HEsurfflxevols[0,0,4,:],'v-',label=expnames[0,4],markeredgecolor=None,markersize=3,color='k')
        ax1EH.plot(HEsurfflxevols[0,1,3,:],'v-',label=expnames[1,3],markeredgecolor=None,markersize=3,color='b')
        ax1EH.plot(HEsurfflxevols[0,3,2,:],'s-',label=expnames[3,2],markeredgecolor=None,markersize=3,color='r')
        ax1EH.legend(loc='upper left',frameon=False)
        ax1EH.plot(HEsurfflxevols[1,0,4,:],'s-',label=expnames[0,4],markeredgecolor=None,markersize=3,color='k')
        ax1EH.plot(HEsurfflxevols[1,1,3,:],'s-',label=expnames[1,3],markeredgecolor=None,markersize=3,color='b')
        ax1EH.plot(HEsurfflxevols[1,3,2,:],'s-',label=expnames[3,2],markeredgecolor=None,markersize=3,color='r')
                
        
        ax1.fill_between(range(16,28),0,-flxevols[0,4,16:28],alpha=0.4,color='k')            
        ax1.fill_between(range(15,26),0,-flxevols[1,3,15:26],alpha=0.4,color='b')   
        ax1.fill_between(range(15,27),0,-flxevols[3,2,15:27],alpha=0.4,color='r')    
            
            
        
        
        ax1.set_ylabel('moisture influx $F_v$ \n in g m$^{-2}$ s$^{-1}$')
        ax1.set_xlabel('time of day')
        ax1EH.set_ylabel('surface heat flux in W m$^{-2}$')
        ax1EH.set_xlabel('time of day')
        ax1EH.grid(True)
        
        #ax1.fill_between(range(15,27),0,-flxevols[1,4,15:27],alpha=0.4)
        
        for a in [ax1,ax1EH]:
            a.grid(True)
                
            a.set_xticks(range(0,48,12))
            a.set_xticklabels(map(daytime,range(0,48,12)))
            
        ax1.legend(loc='lower left',frameon=False)
        
        #second column
        for a in [ax,axEH]:
            a.set_xlabel('volume in $V_0$')
            a.grid(True)
        axEH.set_xlim([0,14])
        
        axEH.set_ylabel('mean surface heat flux in W m$^{-2}$')    
        ax.set_ylabel('total moisture influx $F_{v,tot}$ \n in kg m$^{-2}$')
        
        if ax2.shape==(2,2):
            figname = FIGPATH + 'new_moistinflow.pdf'
        else:
            figname = FIGPATH + 'moistinflow.pdf'
        f.tight_layout()
        if isavefig: f.savefig(figname, dpi=mydpi,bbox_inches='tight')
        

    assert 0
    
   
    iplt_radwind_byh = 0
    if iplt_radwind_byh:    
        hs = ['h250a20','h500a20','h750a20','h1000a20','h1500a20']
        nhs = len(hs)
        f,ax = plt.subplots(3,nhs,sharex=True,sharey=True,figsize=(13,8))
        i=0
        for oro in hs:
            a,b,cu,c=radialwindconv(f,ax[0,i],oro=oro,plotvarn='u')
            a,b,cw,c=radialwindconv(f,ax[1,i],oro=oro,plotvarn='w')
            
            mi,mj = np.where(expnames==oro)
            mi,mj = mi[0],mj[0]
            cftwp,b=plt_LWP_E_contplt(f,ax[2,i],mi,mj)
            
            ax[0,i].set_ylim([0,24])
            i+=1
            
        cbax_u=f.add_axes([0.92,0.65,0.02,0.2])
        cbax_w=f.add_axes([0.92,0.4,0.02,0.2])
        cbax_twp=f.add_axes([0.92,0.1,0.02,0.2])
        plt.colorbar(cftwp,cax=cbax_twp,orientation='vertical',
                         label='TWP in mm')
        plt.colorbar(cu,cax=cbax_u,orientation='vertical',label='u in m/s')
        plt.colorbar(cw,cax=cbax_w,orientation='vertical',label='w in m/s')
            
        for i in range(3): ax[i,0].set_ylabel('local time')
        for j in range(nh): ax[-1,j].set_xlabel('r in km')
            
        figname = FIGPATH + 'all_uwtwp.pdf'
        f.savefig(figname,bbox_inches='tight')
    
    iplt_moist_byh = 1
    if iplt_moist_byh:
        hs = ['h1000a10_nodc','h1000a14_nodc_dry', 'h1500a14_nodc_dry']
        nhs = len(hs)
        
        tinsts = [0,0]
        nts = len(tinsts)
                
        
        f,ax = plt.subplots(nts,nhs,sharex=True,sharey=True,figsize=(8,8))
        i=0
        for i in range(nts):
               
            for j in range(nhs):
                ax[0,j].set_title(hs[j],fontweight='bold')
                f,ax[i,j],cf=plot_uwrz_field(f,ax[i,j],oro=hs[j],tinst=tinsts[i],ik=i,jk=j)
        
        if cf:
            add_customcolorbar(f,cf,'$\Delta q_v$ in g/kg',pos='right')
        figname = FIGPATH + 'all_moist.pdf'
        for j in range(nhs): ax[-1,j].set_xlabel('r in km')
        f.savefig(figname,bbox_inches='tight')
    
    assert 0
    
    iplt_wind_all = 0
    if iplt_wind_all:
        hs = ['h250a14','h500a14','h750a14','h1000a14','h1500a14']
        
        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'

        wsf=3

        ipltvecs = False
        evalrad = 30. 
        i_evalmethodnew = True
        
        BASEps=[BASEnow, BASEw05,BASEw1,BASEw2,BASEw4][:3]
        
        
        nsets=len(BASEps)
        simsetname = ['U0','U5','U10','U20','U40'][:nsets]
        
        tinsts = [10,12]
        nts = len(tinsts)
        levs=[1,5,10,15,20,25,30]
        f,ax=plt.subplots(nts,nsets,figsize=(5,6))


        for i in range(nts):
            for j in range(nsets):
                BASE = BASEps[j]
                ti = tinsts[i]
                expn='h1500a14'
                fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/ensmean_day_3d_d0d5.nc'
                nc = Dataset(fln)
                TPfld = np.sum(nc.variables['TOT_PREC'][:ti,nhalo:-nhalo,nhalo:-nhalo],axis=0)
                
        
                #frainsignal[j,i]=np.mean(TPfld[np.where(TPfld>10.)])
                if i==0 and j==0:
                    mintp=np.mean(TPfld)
    
                h,a=1500,14
    

                
                cf = ax[i,j].contourf(X,Y,TPfld,cmap=viridis_inv,levels=levs,extend='both')

                titlename = expn+'_'+simsetname[j] #+' dommean='+str(np.mean(TPfld[:]))
                print titlename
                
                
                
                c = 256.5
                add_circle(ax[i,j],'h1500m','60_homo',c,a)
                # add analysis evaluation circle
                i_evalmethodnew=True
                if i_evalmethodnew:
                    strw = evalrad/3.
                    centr=256
                    TPfldt = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)
                    tpstrip = TPfldt[centr-strw:centr+strw,:]
                    tpdist = np.mean(tpstrip,axis=0)
                    cmov = np.argmax(tpdist)
                    cx = X[centr,cmov]
                    cy = Y[centr,cmov]
                    circoff = plt.Circle((cx,cy),30.,color='red',fill=False,linestyle='dashed',linewidth=1.5)
                  
                    #ax[i,j].add_artist(circoff)    
                    #ax[i,j].plot(cx,cy,'x',markersize=9,color='red')        
                nc.close()
                fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/ensmean_day_3d_d0d5.nc'
                
                nc = Dataset(fln)
                print fln
                
                U=nc.variables['U'][ti,-1,nhalo:-nhalo,nhalo:-nhalo]#-np.mean(nc.variables['Uinterp'][20,-6,nhalo:-nhalo,nhalo:-nhalo])
                V=nc.variables['V'][ti,-1,nhalo:-nhalo,nhalo:-nhalo]
                W=nc.variables['W'][ti,-1,nhalo:-nhalo,nhalo:-nhalo]
                #X=nc.variables['lon'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                #Y=nc.variables['lat'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                print X.shape
                st = 4. #int(a/5)
                qv = ax[i,j].quiver(X[::st,::st],Y[::st,::st],U[::st,::st],V[::st,::st],
                                    units='inches',scale=20, color='grey',pivot='mid',angles='xy',
                                    headwidth=2.5,width=0.008)
                
                
                tcks = np.arange(c-10*a,c+10*a+1.,a)
                ax[i,j].set_xticks(tcks)
                ax[i,j].set_yticks(tcks)
                
                tckl = np.arange(-10,10+1,1)
                tckl = np.array(tckl,dtype=object)
                tckl = map(lambda x: x if x%2==1 else '',tckl)
                ax[i,j].set_xticklabels(tckl)
                ax[i,j].set_yticklabels(tckl)
                decorate_ax(ax[i,j],"",daytime(ti),"")
                ax[i,j].tick_params(axis='both',which='major',labelsize=6)
                ax[i,j].grid(True)
                wsfl = 1
                
                cxnew = cmov #(c+cmov)/2.
                ax[i,j].set_xlim([cxnew-wsf*a,cxnew+wsf*a])
                ax[i,j].set_ylim([c-wsf*a,c+wsf*a])
                ax[i,j].set_aspect('equal')
                #ax[i,j].tick_params(axis='both', which='major', labelsize=6)

        ax[0,0].quiverkey(qv,0.5,1,2,r'$1 \frac{m}{s}$',labelsep=0.05,labelpos='N')
        for i in range(nts): ax[i,0].quiverkey(qv,-1.3,0.5,2,r'$1 \frac{m}{s}$',
                                                labelsep=0.05,labelpos='N')
        for i in range(nsets):
            ax[-1,i].set_xlabel('x/a')
        for i in range(nts): ax[i,0].set_ylabel('y/a')
        add_customcolorbar(f,cf,'accumulated rain amount in mm',pos='right',orient='vertical')
        
        for j in range(nsets):
            ax[0,j].set_title(simsetname[j],fontsize=10,fontweight='bold')
        
        figname=FIGPATH+'Ufld'+expn+'.pdf'
        #f.savefig(figname,dpi=100,bbox_inches='tight')
        

    iplt_windcs_all = 0
    if iplt_windcs_all:

        
        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'

        wsf=3

        ipltvecs = False
        evalrad = 30. 
        i_evalmethodnew = True
        
        BASEps=[BASEnow, BASEw05,BASEw1,BASEw2,BASEw4][:3]
        
        nsets=len(BASEps)
        simsetname = ['U0','U5','U10','U20','U40'][:nsets]
        
        tinsts = [12,18,22]
        nts = len(tinsts)
        levs=[1,5,10,15,20,25,30]
        f,ax=plt.subplots(nts,nsets,sharex=True,sharey=True,figsize=(12,8))


        for i in range(nts):
            for j in range(nsets):
                BASE = BASEps[j]
                ti = tinsts[i]
                expn='h3000a14'

                
                fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/ensmean_day_3d_d0d5.nc'
                
#                nc = Dataset(fln)
#                print fln
#                nc.close()
                a,b,cf = u_csect(f,ax[i,j],oro=expn,BASE=BASEps[j],t=tinsts[i])
                decorate_ax(ax[i,j],"",daytime(tinsts[i]),"")

        av = 14.

        for j in range(nsets):
            ax[-1,j].set_xlabel('x in km')
            
            ax[0,j].set_title(simsetname[j],fontweight='bold')
        for i in range(nts):
            ax[i,0].set_ylabel('z in km')  
        
        add_customcolorbar(f,cf,mylabel='$u^{\prime}$ in m/s',pos='right')                
        
        figname=FIGPATH+'U_xzfld'+expn+'.pdf'
        if isavefig: f.savefig(figname,dpi=300,bbox_inches='tight')
        
    

    
    iplt_convvariab = 0
        
    if iplt_convvariab:
        expn = 'h1500a14'
        
        ishade = 1
        if ishade:
            f,ax2=plt.subplots(1,2,gridspec_kw = {'width_ratios':[3, 1]},sharey=True)   
            ax,axs = ax2
            ensmemdata = np.zeros([5,10])
        else:
            f,ax2=plt.subplots(1,3,sharey=True)
            ax,axh,axmn=ax2
       
        nos= 5
        expns = ['h1000a14','h2000a14','h3000a14'] #expnames[:nos,4]
        orocol = dict(zip(expns,['r','g','b','k','cyan'][:nos]))
        iexpn=0
        for expn in expns:
            srcpath=BASEnow+'TOT_PREC/'+expn +'/60_homo/'
            datalist = filter(lambda x: x.endswith('daysum_TP.nc'),os.listdir(srcpath))[:10]
            
            
            #if iexpn==0:
            #    mask = R>287.
            #else:
            mask = R<30.       
            j=0
            daymean = np.zeros(5)
            for datan in datalist:
                filepath = srcpath +datan
                data = Dataset(filepath)
                TP = data.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo]

                print TP.shape
                TPmn = np.mean(TP[:,mask],axis=1)
                print TPmn.shape
                
                if ishade:

                    ensmemdata[:,j] = TPmn
                    
                else:
                    ax.plot(TPmn,label=str(j),markersize=8,color=orocol[expn],alpha=0.5)
               
                    daymean+=TPmn
                j+=1
                data.close
            if ishade:
                datamin = np.min(ensmemdata,axis=1)
                datamax = np.max(ensmemdata,axis=1)
                ensmn = np.mean(ensmemdata,axis=1)
                ensmed = np.median(ensmemdata,axis=1)
                
                
                ax.fill_between(range(5),datamin,datamax,color=orocol[expn],alpha=0.3)                
                ax.plot(range(5),ensmn,color=orocol[expn],marker='s',linewidth=2,markeredgewidth=0)
                ax.plot(range(5),ensmed,color=orocol[expn],marker='s',
                        markeredgewidth=0,linestyle='--',linewidth=2)
                ndata=np.ones(ensmemdata.shape)*iexpn
                
                iplt_reduced = 1
                if iplt_reduced:
                    if iexpn==0: fred,axred = plt.subplots(1,1,figsize=(3,4))
                    axred.scatter(ndata,ensmemdata,color=orocol[expn],alpha=0.3)
                    axred.plot([iexpn-0.2,iexpn+0.2],np.ones(2)*np.mean(ensmn),linewidth=3,
                             color=orocol[expn])
  
                
                axs.scatter(ndata,ensmemdata,color=orocol[expn],alpha=0.3)                
                #axs.scatter(ndata,ensmemdata,color='white')
                #for ida in range(5):
                #    for im in range(10):
                #        axs.text(1,ensmemdata[ida,im],str(ida+1),color=orocol[expn])
                axs.plot([iexpn-0.2,iexpn+0.2],np.ones(2)*np.mean(ensmn),linewidth=3,
                         color=orocol[expn])
                iexpn+=1
                
                # mean line
                ax.text(1.25,np.mean(ensmn),expn,
                         color=orocol[expn],verticalalignment='center')
                
                
            else:
                axh.plot(daymean/(j+1),marker='o',color=orocol[expn])
                axmn.plot(1,np.mean(daymean/(j+1)),marker='s',color=orocol[expn])
                axmn.plot(1,np.median(daymean/(j+1)),linestyle='--',marker='s',color=orocol[expn])
                axmn.text(1.01,np.mean(daymean/(j+1)),expn,
                          verticalalignment='center')
        if ishade:
            #ax.plot(2,15.5,color='white')
            for a in ax2: a.set_ylim([0,16])
            adjust_spines(ax,['left','bottom'])
            ax.set_xticks(range(5))        
            ax.set_xticklabels(range(1,6))
            ax.set_xlabel('day')
            ax.set_ylabel('rain amount in mm/d')
            ax.grid(True)      
            
            #stat axis
            axs.set_xticks(range(nos))
            
            #axs.set_xlabel('')
            axs.grid(True)
            adjust_spines(axs,['left','bottom'])
            #axs.set_xlim([-0.5,3])
        
        
        else:
            for a in [ax,axh]:
                a.plot(2,15.5,color='white')
                adjust_spines(a,['left','bottom'])
                a.set_xticks(range(5))        
                a.set_xticklabels(range(1,6))
                a.set_xlabel('day')
                a.grid(True)
                if a==ax: a.set_ylabel('rain amount in mm/d')
            
            axmn.set_xticks([1.])
            axmn.set_xlim([0.9,1.1])
            ax.set_title('all members')
            axh.set_title('ensemble mean')
            axmn.set_title('day mean')        
        
        if iplt_reduced:
            adjust_spines(axred,['left','bottom'])
            axred.set_xticks(range(3))        
            axred.set_xticklabels(expns,rotation=45,position=(-0.5,0.01))
            axred.set_ylabel('rain amount in mm/d')
            axred.grid(True) 
            axred.set_ylim([0,16])
            axred.grid(False)
            figname = FIGPATH +'vconst_ensspread_estimate_2.pdf'
            fred.savefig(figname, bbox_inches='tight',dpi=mydpi)
                
        axs.set_xticklabels(expns,rotation=45,position=(-0.5,0.01))
        
                
        figname = FIGPATH +str(nos)+'_vconst_ensspread_estimate.pdf'
        print figname
        if isavefig: f.savefig(figname,bbox_inches='tight',dpi=mydpi)
   
    
    iplt_u_cs = 0 
    if iplt_u_cs:
        revno = 'rev1'
            
        orobyrevs = {'rev0':['h500a14','h1500a14'],'rev1':['h500a14','h1500a14','h3000a14']}
        oros = orobyrevs[revno]
        
        wndexps = [BASEnow,BASEw05,BASEw1]  
        Ucode = ['U0','U5','U10']
        no,nw = len(oros),len(wndexps)
        f,ax=plt.subplots(no,nw,figsize=(5,5),sharex=True, sharey=True)
        for i in range(no):
            for j in range(nw):
                a,b,cf=u_csect(f,ax[i,j],oro=oros[i],BASE=wndexps[j],t=24)
                simname = oros[i]+' '+Ucode[j]
                ax[i,j].set_title(simname,fontweight='bold',fontsize=8)
        av=14.
        for j in range(nw):
            ax[-1,j].set_xlabel('x/a')
            ax[-1,j].set_xticks(np.arange(-2*av,10*av,2*av))
            ax[-1,j].set_xticklabels(map(int,np.arange(-2*av,10*av,2*av)/av))
        for i in range(no):
            ax[i,0].set_ylabel('z in km')  
        
        add_customcolorbar(f,cf,mylabel='$u^{\prime}$ in m/s',pos='right')
        
        figname = FIGPATH+revno+'_u_fld.pdf'
        
        if isavefig: f.savefig(figname,bbox_inches='tight')


    iplt_oro1 = 0
    if iplt_oro1:
        f,ax12=plt.subplots(1,2,sharex=True,sharey=True,figsize=(6,3))
        ax,axw = ax12
        x = np.arange(-5,5,0.01)
        noise = 0.5*np.cos(2*np.pi*x/1)**2+0.5
        #noise[np.where(np.abs(x)>2)] = 0.   
        h1 = np.exp(-x**2/0.7) 
        h1*=noise
        ax.set_ylabel('height in m')
        ax.fill_between([-2,2],0,2,color='red',alpha=0.2)
        #ax.plot((-2,-2),(0,3),color='red',linestyle='dashed',linewidth=2)
        #ax.plot((2,2),(0,3),color='red',linestyle='dashed',linewidth=2)
        rect1=plt.Rectangle((-2,0),4,2,color='k',alpha=1,linestyle='dashed',
                            linewidth=2.5,fill=False)         
        ax.add_artist(rect1)
         
        x1 = np.array([-2,3])
        axw.fill_between([-2,2],0,2,color='red',alpha=0.1)
        rect2=plt.Rectangle((0,0),4,2,color='k',alpha=1,linestyle='dashed',
                            linewidth=2.5,fill=False)   
        axw.add_artist(rect2)
        
        for a in ax12: 
            a.plot(x,h1,color='grey')
          
            a.text(0,1.5,'$V_{c}-V_{m}$',color='red',fontweight='bold',fontsize=12,horizontalalignment='center')
            a.text(0,0.25,'$V_{m}$',color='black',fontweight='bold',fontsize=12,horizontalalignment='center')
            a.fill_between(x,0,h1,color='grey')
            a.set_ylim([0,2])
            a.set_xticks([0])
            a.set_yticks([])
            a.set_xlabel('r in km')
            a.set_yticklabels([])
            a.set_xticklabels([])
            #ax.set_title(titles[i],fontweight='bold')
            a.set_ylabel('height in m')
              
            a.set_xlim([-3,5])
            

        figname = FIGPATH + 'vscal_sketch.pdf'
        f.tight_layout()
        f.savefig(figname)


    iplt_sketch_small_big_scales = 0
    if iplt_sketch_small_big_scales:
        f,ax12=plt.subplots(1,2,sharex=True,sharey=True,figsize=(6,3))
        ax,axw = ax12
        x = np.arange(-5,5,0.01)
        noise = 0.5*np.cos(2*np.pi*x/0.8)**2+0.5
        #noise[np.where(np.abs(x)>2)] = 0.   
        h1 = np.exp(-x**2/0.7) 
        h2 = np.copy(h1)*0.8
        
        h1*=noise
        
        

        ax.fill_between(x,0,h1,color='grey')
        axw.fill_between(x,0,h2,color='grey')

        for a in ax12:
            a.set_ylim([0,2])
            a.set_xlim([-2,2])
            a.set_xticks([])
            a.set_yticks([])
            adjust_spines(a,[])

        figname = FIGPATH + 'vscal_smallscale_features_sketch.pdf'
        f.tight_layout()
        f.savefig(figname)
      

    iplt_orosketsch = 0
    if iplt_orosketsch:
        f,ax=plt.subplots(1,4,sharex=True,sharey=True,figsize=(10,2))
        x = np.arange(-5,5,0.001)
        h1 = np.exp(-x**2/0.7) 
        h2 = 0.5*np.exp(-x**2/2)
        hcos = 0.7*np.cos(-x**2/1.3)**2.
        hcos[np.where(np.abs(x)>1.4)]=0.
        
        hmulti = 0.5*(np.exp(-(x-1.2)**2/0.7)+np.exp(-(x+1.2)**2/0.7))
        
        noise = 0.01*np.cos(2*np.pi*x/0.2)
        noise[np.where(np.abs(x)>2)] = 0.        
        
        titles = ['(a)','(b)','(c)','(d)']        
        
        h3 = h2+noise
        hs = [h1,h2,hcos,hmulti]
        for i in range(4):
            ax[i].plot(x,hs[i],color='grey')
            ax[i].fill_between(x,0,hs[i],color='grey')
            ax[i].set_ylim([0,2.])
            ax[i].set_xticks([0])
            ax[i].set_yticks([])
            ax[i].set_xlabel('r in km')
            ax[i].set_yticklabels([])
            ax[i].set_xticklabels([])
            ax[i].set_title(titles[i],fontweight='bold')
        ax[0].set_ylabel('height in m')
        figname = FIGPATH + 'sketch.pdf'
        f.tight_layout()
        f.savefig(figname)

    iplt_theta = 0
    if iplt_theta:
        ipltsym=1
        tinsts=[12,20,24,28]
        nt = len(tinsts)
        oro='h500a14'
        #plt.rcParams['axes.facecolor']='grey'
        f,ax=plt.subplots(1,nt,figsize=(10,4),sharex=True,sharey=True)
        i=0
        for t in tinsts:
            print t
            a,b,cf=plot_theta(f,ax[i],oro=oro,t=t,ipltsym=ipltsym)
            decorate_ax(ax[i],"",daytime(t),"")
            
            ax[i].set_xticks(np.arange(0,41,10))
            ax[i].set_xlabel('r in km')
            i+=1
        ax[0].set_ylabel('z in km')
        add_customcolorbar(f,cf,'$\Theta$ in K',orient='vertical',pos='right')
        figname=FIGPATH+oro+'theta_field.pdf'        
        #f.tight_layout()
        if isavefig: f.savefig(figname,bbox_inches='tight',dpi=mydpi)
    plt.rcParams['axes.facecolor']='white'


    iplt_E_comp = 0
    if iplt_E_comp:
        BASEs = [BASEnow,BASEw05,BASEw1,BASEw2,BASEw4]
        times = [12, 20, 24, 28]
        nt = len(times)
        i,j = 3,0
        nsets=len(BASEs)
        f,axi = plt.subplots(1,nt)
        #for i in range(nsets):
        for j in range(nt):
            ax = axi[j]            
            ti = times[j]        
            oron = 'h500a20'
            srcpath = BASEs[i]+'/ALLVAR_3D/'+oron+'/60_homo/ensmean_day_3d_d0d5.nc'        
            print srcpath        
            data = Dataset(srcpath)
    
    
            offset = 128
            ncut = nhalo+offset
    
            U = data.variables['U'][ti,-1,ncut:-ncut,ncut:-ncut]
            V = data.variables['V'][ti,-1,ncut:-ncut,ncut:-ncut]
    
            X = data.variables['lon'][ncut:-ncut,ncut:-ncut]*CONV-256.
            Y = data.variables['lat'][ncut:-ncut,ncut:-ncut]*CONV-256.
            data.close()
    
    
            print np.mean(U**2+V**2)
    
    
            stz = 10
            Q = ax.quiver(X[::stz,::stz],Y[::stz,::stz],U[::stz,::stz],V[::stz,::stz],
                       scale=60,angles='xy', color='red',pivot='mid',
                       headwidth=2,width=0.005)
                    
            ax.quiverkey(Q,0,-70,1,r'$1 \frac{m}{s}$',labelsep=0.05,
                         labelpos='N',coordinates='data')
            ax.set_aspect('equal')   

        
        


    iplt_tpfld_U0_selection = 0
    if iplt_tpfld_U0_selection:
        
        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASE = BASEnow
        oros = ['h250a20','h500a14','h1000a10']
        ha_vals = {'h250a20':(250,20),'h500a14':(500,14),'h1000a10':(1000,10)}
        nsets = 3 
        levs=np.linspace(1,40,7)
        f,ax=plt.subplots(1,nsets)

        i = 0
        evalrad = 30.
        
        for oro in oros:
            expn = oro
            fln = BASE+'/TOT_PREC/'+expn+'/60_homo/ensmean.nc'
            nc = Dataset(fln)
            TPfld = np.mean(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/2.
            nc.close()
            h,a=ha_vals[oro]

            cf = ax[i].contourf(TPfld,levels=levs,cmap=viridis_inv,extend='both')
            c = 255.5
            add_circle(ax[i],'h500m','60_homo',c,a)
            # add analysis evaluation circle
            
            strw = evalrad/3.
            centr=c#s256
            tpstrip = TPfld[centr-strw:centr+strw,:]
            tpdist = np.mean(tpstrip,axis=0)
            cmov = np.argmax(tpdist)
            cmov = c
            cx = X[centr,cmov]
            cy = Y[centr,cmov]
            circoff = plt.Circle((cx,cy),30.,color='red',fill=False,linestyle='dashed',linewidth=1.5)
            ax[i].add_artist(circoff)    
            ax[i].plot(cx,cy,'x',markersize=9,color='red')


            wsf = 2
            tcks = np.arange(c-wsf*a,c+wsf*a+1.,a)
            tckl = np.arange(-wsf,wsf+1,1)
            ax[i].set_xticks(tcks)
            ax[i].set_yticks(tcks)
            ax[i].set_xticklabels(tckl)
            ax[i].set_yticklabels(tckl)
            
            ax[i].set_xlim(c-30,c+30)
            ax[i].set_ylim(c-30,c+30)
            ax[i].set_aspect('equal')
            ax[i].grid(1)
            
            if 0:
                wsf=6
    
                
                
                tcks = np.arange(c-wsf*a,c+wsf*a+1.,a)
                ax[i].set_xticks(tcks)
                ax[i].set_yticks(tcks)
                
                tckl = np.arange(-wsf,wsf+1,1)
                tckl = np.array(tckl,dtype=object)
                tckl = map(lambda x: x if x%2==1 else '',tckl)
                ax[i].set_xticklabels(tckl)
                ax[i].set_yticklabels(tckl)
               
                ax[i].tick_params(axis='both',which='major',labelsize=6)
                ax[i].grid(True)
                wsfl = 1
                
                cxnew = (c+cmov)/2.
                ax[i].set_xlim([cxnew-wsf*a,cxnew+wsf*a])
                ax[i].set_ylim([c-wsf*a,c+wsf*a])
                
            i+=1



        for i in range(nsets):
            ax[i].set_xlabel('x/a')
        #for i in range(na): ax[i].set_ylabel('y/a')
        add_customcolorbar(f,cf,'rain amount in mm',pos='right',orient='horizontal')
        
        figname = FIGPATH+'Fx_TPfld_orosel.pdf'
        
        print figname
    
        if isavefig: f.savefig(figname,dpi=mydpi,bbox_inches='tight')  



        

    
    iplt_circ_sel = 0
    if iplt_circ_sel:
        """ Snippet for the production of figure2.pdf """
        
              
        
        iconstvol = 0
        if iconstvol:
            oros = ['h250a20','h500a14','h1000a10']
            h_s = [0.25,0.5,1]
            a_s = [20,14,10]
            tinsts = [18,22,26]
        else:
            oros = ['h500a14','h500a14_cos2']
            oronames = ['h500a14 (Gauss)','h500a14 (cos$^2$)']
            tinsts = [20,22,24]
            tinsts = [18,22,26]
        noro = len(oros)
        
        
        ntinsts = len(tinsts)

        imatr = True
        if imatr:
            f,ax=plt.subplots(ntinsts,noro,sharex=True,sharey=True,figsize=(8,8))
        else:
            f,ax=plt.subplots(1,noro,sharex=True,sharey=True)
        
        if imatr:
            for i in range(ntinsts):
                j=0            
                for o in oros:
                    f,ax[i,j],cf=plot_uwrz_field(f,ax[i,j],oro=o,tinst=tinsts[i],ik=i,jk=j)
                    j+=1
            for j in range(noro):
                if iconstvol:
                    ax[0,j].set_title(oros[j])
                else:
                    ax[0,j].set_title(oronames[j])
             
            if iconstvol:
                ax[0,0].set_xticks([0,20,40])
                ax[0,0].set_yticks([0.,1,2,3,4,5])
        else:
            for i in range(noro):
                f,ax[i],cf=plot_uwrz_field(f,ax[i],oro=oros[i],tinst=tinsts[i],jk=i)
                ax[i].set_title(oros[i],fontweight='bold')
                ax[i].set_xlabel('r in km')
                
    
        
        if isavefig: 
            #f.tight_layout() 
            add_customcolorbar(f,cf,'$\Delta q_v$ in g/kg',pos='right')
            pref = 'matr_' if imatr else 'seltime_'
            pref2 = '' if iconstvol else 'gauss_cos2_'
            figname = FIGPATH+pref+pref2+'F2.pdf'
            f.savefig(figname,dpi=mydpi,bbox_inches='tight')
            


    iplt_EKIN_evol = False
    if iplt_EKIN_evol:
        
        i,j=3,2
        f,ax=plt.subplots(1,1)
        expn = expnames[i,j]        
        basepath = BASE3D+expn+'/60_homo/circmean_day_d0d5.nc'
        
        print basepath        
        
        ncfl = Dataset(basepath,'r')
        s = ncfl.variables['Speedrz'][:]
        rho = ncfl.variables['RHOrz'][:]  
        r = ncfl.variables['X'][:]
        z = ncfl.variables['Z'][:]
        
        ncfl.close()
        
        Ek = 0.5*rho*(s**2)
        Eksum = np.sum(Ek[:]*z*r,axis=(1,2))
        plt.plot(Eksum)
        
        
        
    
    

    iplt_xt_of_r = False
    
    if iplt_xt_of_r:
        BASE = BASETP  
        f,ax1 = plt.subplots(1,4,sharex=True,sharey=True)
        for i in range(nh):
            for j in range(4):
        
                ax = ax1[j]
                expn = expnames[i,j]
                basepath = BASE+expn+'/60_homo/'      
                print expn
                
                ensmems = filter(lambda x: x.endswith('series_TP.nc'), os.listdir(basepath))        
                print ensmems  
                
                rs = np.arange(0,128,5)
                tp_xt = np.zeros(len(rs)-1)
                       
                rvals = 0.5*(rs[1:]+rs[0-1])
                
                for testmem in ensmems[:1]:
                    testmem = 'dc_TP.nc'
                    ncfl = Dataset(basepath+testmem)
                    tp = ncfl.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo]    
                     
                                
                    for k in range(len(rs)-1):
                        hi = rs[k+1]
                        lo = rs[k]
                        mask = np.logical_and(R>=lo,R<hi)
                        
                        tp_xt[k] = np.max(tp[:,mask])
                        
           
                    ax.plot(rvals,tp_xt)
                    ax.set_title(expnames[i,j])
                    ncfl.close()
    


    iplt_ref = 0
    if iplt_ref:
        orosel = 'h1500a14'
        for t in [23]:#range(36,49,2):
            selvar='QVrz'
            f,ax=plt.subplots(1,1)
            a,b,cf=plot_uwrz_field(f,ax,oro=orosel,tinst=t,var=selvar)
            add_customcolorbar(f,cf,'q_v in g/kg',pos='right')
            tstamp = str(t) if t>9 else '0'+str(t)
            figname=FIGPATH+'/'+orosel+'/abs_'+selvar+'_t'+tstamp+'.pdf'
            f.savefig(figname)
                
      
    iplt_timevol_circs = False
    if iplt_timevol_circs:
        orosel = 'h500a20'
        trange = range(22,27,2)
        for ti in trange:        
            nts = 1 #len(trange)
            
            
            #i = 0
            #for t in trange: #range(36,49,2):
            selvar='QVrz'
            f,ax = plt.subplots(1,nts,sharex=True,sharey=True)
            f,ax,cf=plot_uwrz_field(f,ax,oro=orosel,tinst=ti,var=selvar)
            add_customcolorbar(f,cf,'$\Delta q_v$ in g/kg',pos='right')
            tstamp = str(t) if t>9 else '0'+str(t)
            
                
            figname=FIGPATH+'movie/demo_'+str(ti)+'_'+orosel+'.png'   
            print figname     
            if isavefig: f.savefig(figname,dpi=mydpi)

    
        
    # Plot the surface flux evolution
    iplt_srfcflxevol = False
    if iplt_srfcflxevol:
        
        # sep by h
        fa,axh = plt.subplots(1,nh,sharey=True)
        for i in range(nh):
            for j in range(na):
                expn=expnames[i,j]
                flname = BASE3D+expn+'/60_homo/'+CIRCMEAN
                ncfl = Dataset(flname,'r')
                E = ncfl.variables['EFLUXrz'][20,-1,:]
                RHO = ncfl.variables['RHOrz'][20,-1,:]
                #E = E*RHO
                Etot = np.sum(E,axis=0)
                axh[i].plot(Etot,label=expn)
                ncfl.close()
                axh[i].plot()
            axh[i].legend()
        
        # sep by a
        fa,axa = plt.subplots(1,na,sharey=True)
        for i in range(nh):
            for j in range(na):
                print 'temp'
                

    
    # Plot circulations
    iplt_circulations = False
    if iplt_circulations:
        f,ax,cf = plot_uwrz_field(tinst=20,var='QVrz')
        add_customcolorbar(f,cf,'qc',pos='right')
        mynh=nh
        myna=na

        for tinst in [20]:
            f,ax = plt.subplots(mynh,myna,sharex=True,sharey=True)   
            for i in range(mynh):
                for j in range(myna):
                    oro = expnames[i,j]
                    a,b,cf = plot_uwrz_field(f,ax[i,j],oro=oro,tinst=tinst,var='QVrz')
                ax[i,0].set_ylabel('z in km')
                ax[-1,j].set_xlabel('r in km')
            
            add_customcolorbar(f,cf,'qv in kg/kg',pos='right')
            figname = FIGPATH+str(tinst)+'_'+str(mynh)+'_x_'+str(myna)+'_qvrz.pdf'
            print figname
            if isavefig: f.savefig(figname,dpi=mydpi,bbox_inches='tight',papertype='a3')



    
    
    iplt_EP_EmP = False
    if iplt_EP_EmP:       
        # Evap minus precip plot:
        Lv = 2.264705*1e6 # latent heat of evaporation for water
        delt_out = 1800 # output time step
        
        # E, P, and MF (mass flux) matrix
        Ematrix = np.zeros([nh,na])
        Pmatrix = np.zeros([nh,na])
        MFmatrix = np.zeros([nh,na])
        
        # iterable from here        
        i,j=4,5        
        expn = expnames[i,j]  
        print expn
        a = arange[j]
        
        
        # Extract the fields from Data
        fln = BASE3D+expn+'/60_homo/'+CIRCMEAN
        ncfl = Dataset(fln,'r')
        P = np.sum(ncfl.variables['TOT_PRECr'][:],axis=0)
        Z = ncfl.variables['Z'][:,a]
        efluxrz = ncfl.variables['EFLUXrz'][:,-1,:]
        urz = ncfl.variables['Urz'][:,:,a]
        rhorz = ncfl.variables['RHOrz'][:]
        qvrz = ncfl.variables['QVrz'][:]
        ncfl.close()
        
        # vertical coordinates at r = a
        Zfull = Z
        Zhalf = 0.5*(Z[1:]+Z[:-1])
        Zthick = 0.5*(Z[1:]-Z[:-1])
        print Z.shape,qvrz.shape,urz.shape,rhorz.shape
        
        #moistflux = rhorz[:,:,a]*qvrz[:,:,a]*urz*delt_out
        
        E = efluxrz[:48]*1800/Lv
        E = np.sum(E,axis=0)
        
        rs = np.arange(0.5,a+0.5,1) 
        
        Emean = np.sum(E[:a]*rs)
        Pmean = np.sum(P[:a]*rs)
        Ematrix[i,j] = Emean
        Pmatrix[i,j] = Pmean     
        
        
        plt.plot(E,label='E')
        plt.hold(True)
        plt.plot(P,label='P')
        plt.legend()
        
        
    iplt_windvsrain = 1
    if iplt_windvsrain:
        """ Documentation!"""
        
        
        t0,t1=[12,36]
        plotvarn='u'     
        
        iorosel=True
        
        
        if iorosel:
            sel = 'hi'
            orogbc = ['h500a20','h500a20_bell','h500a20_cos2']
            
            oroexpsels = {'hi':['h1000a14','h3000a14'],
                            'norm':['h250a20','h500a14','h1000a10'], 
                          'gbc':['h500a20','h500a20_bell','h500a20_cos2'],
                          'ext':['h750a30','h1500a20','h3000a14']}            
            
            oros = oroexpsels[sel]
            norosel = len(oros)
            nr = 2
            ipltELWP = 1            
            if ipltELWP: nr+=1
            
            f,ax = plt.subplots(nr,norosel,sharex=True,sharey=True,figsize=(8,8))
            i,j = 0,0
            
            for plotvarn in ['u','w']:
                j=0
                for oro in oros:
                    trange=[12,37]
                    f,b,cf1,cf2 = radialwindconv(f=f,ax=ax[i,j],oro=oro,
                                                 t=trange,plotvarn=plotvarn)
                    if 0:
                        # plot colorbar old style                         
                        cbar_pos = 'vertical' if plotvarn=='u' else 'horizontal'
                        add_customcolorbar(f,cf1,plotvarn+" in m/s",orient=cbar_pos,pos='right')                             
                        ax[i,j].grid(True)
                    else:
                        # plot colorbar new style
                        if plotvarn=='u':
                            cbax=f.add_axes([0.92,0.65,0.02,0.2])
                        else:
                            cbax=f.add_axes([0.92,0.4,0.02,0.2])
                        plt.colorbar(cf1,cax=cbax,orientation='vertical',
                                     label=plotvarn+' in m/s')
                                     
                    ax[i,j].grid(True)
                    j+=1
                i+=1
            if ipltELWP:
                for ci in range(norosel):
                    idbg=False
                    if idbg:
                        mi,mj = np.where(expnames==oros[ci])
                        mi,mj = mi[0],mj[0]
                        print mi,mj
                    oro=oros[ci]    
                    cftwp,b=plt_LWP_E_contplt(f,ax[-1,ci],oro)
                    cbax=f.add_axes([0.92,0.1,0.02,0.2])
                    plt.colorbar(cftwp,cax=cbax,orientation='vertical',
                                     label='$W_{tot}$ in mm')
            
            for a in ax[-1,:]:
                    a.set_xlabel('$r$ in km')
                    a.set_xticks([0,10,14,20,40])
                    
            for a in ax[:,0]: a.set_ylabel('local time')
   
            ax[0,0].set_xticks(np.arange(0,128,20))
            ax[0,0].set_xlim([0,40])
            ax[0,0].set_ylim([0,24])
            ist = 6
            ax[0,0].set_yticks(np.arange(t0,t1+1,ist)-t0)
            ax[0,0].set_yticklabels(map(daytime,range(t0,t1+1,ist)))
            #add_customcolorbar(f,cf1,plotvarn+" in m/s",pos='right')            
            for i in range(len(oros)):ax[0,i].set_title(oros[i],fontweight='bold')
            
            figname= FIGPATH+plotvarn+'_'+sel+'_Vconst_rzplot.pdf'
            #figname = FIGPATH +'F3.pdf'
            print figname
            f.savefig(figname,bbox_inches='tight',dpy=mydpi)
        else:
            f,ax = plt.subplots(nh,na,sharex=True,sharey=True)
            for i in range(nh):
                for j in range(na):
                    trange=[12,37]
                    f,b,cf1,cf2 = radialwindconv(f=f,ax=ax[i,j],oro=expnames[i,j],
                                                 t=trange,plotvarn=plotvarn)
                    ax[i,j].grid(True)
        

            for i in range(nh):
                ax[i,0].set_ylabel('local time')
            
            for j in range(na):
                ax[-1,j].set_xlabel('r in km')
                ax[0,j].set_title(oros[j],fontsize=18)
            
            ax[0,0].set_xticks(np.arange(0,128,20))
            ax[0,0].set_xlim([0,40])
            ist = 6
            ax[-1,0].set_yticks(np.arange(t0,t1+1,ist)-t0)
            ax[-1,0].set_yticklabels(map(daytime,range(t0,t1+1,ist)))
            add_customcolorbar(f,cf1,plotvarn+" in m/s",pos='right')
                
            add_customcolorbar(f,cf2,"rain rate in mm/h",orient='horizontal',pos='right')
            figname = FIGPATH +orosel+plotvarn+'_windsandrain_.pdf'
            figname = FIGPATH +'F3.pdf'
            print figname
            #f.tight_layout()
            f.savefig(figname)
            


    iplt_LWPvsTP = False
    if iplt_LWPvsTP:
        
        rainamounts = np.zeros(expnames.shape)
        LWP_tots = np.zeros(expnames.shape)
        f,ax=plt.subplots(1,1)
        myhrange = range(nh)
        myarange = range(na)
        for i in myhrange: #range(nh):
            
            for j in myarange:#[1]:#[1,2,3,4,5]: #myarange:
                if not ((i,j) in [(1,3)]):#[(0,4),(1,3),(3,2)]):
                    continue
                    
                expn = expnames[i,j]
                print expn
                hm = hrange[i]          
                a = arange[j]
                
                
                
                analyloc = True
                if analyloc:
                    mask = R<=a
                    masks = Rs<=a
                else: 
                    ri = a
                    ro = a*np.sqrt(2)
                    mask = np.logical_and(R>ri,R<=ro)
                    masks = np.logical_and(Rs>ri,Rs<=ro)

                # rain amount
                flnprefix = 'seed'
                basepath = BASE3D+expn+'/60_homo/'
                fln = filter(lambda x: x.startswith(flnprefix),os.listdir(basepath))[0]
                basepath += fln
                ncfl = Dataset(basepath)    
                TP = ncfl.variables['TOT_PREC'][:]/2.
                
                print TP.shape
                
                M_EVAP = ncfl.variables['EFLUX'][:,-1,:]*1800./Lv
                #M_EVAP = np.cumsum(M_EVAP[:,masks], axis=0)
                M_EVAP = np.mean(M_EVAP[:,masks],axis=1)
                
                print M_EVAP.shape
                
                rainamount_sum = np.sum(TP[:,mask],axis=0)
                rainamount = np.mean(rainamount_sum)
        
                rainamounts[i,j] = rainamount
                ncfl.close()
               
                # LWP total
                basepath = BASE3D+expn+'/60_homo/'+LWPN
                ncfl = Dataset(basepath)
                
                
                LWP_tot=np.zeros((48,256,256))
                LWP_tot=LWP_tot[:,masks]
                
                print LWP_tot.shape
        

                # Summing up to LWP_total
                for var in ncfl.variables.values():
                    var = var.name
                    if not (var.startswith('LWP')):
                        continue
                    print var
                    LWP_temp = ncfl.variables[var][:]
                    LWP_temp = LWP_temp[:,masks]
                    #plt.plot(np.mean(LWP_temp,axis=1),label=var)
                    #plt.xlim([12,36])
                    LWP_tot += LWP_temp
                
                iplt_LWP_evol = 0
                if iplt_LWP_evol:
                    plt.plot(M_EVAP,label='E')
                    plt.plot(np.cumsum(np.mean(TP[:,mask],axis=1)),label='TP')
                    plt.legend()
                    plt.savefig(FIGPATH+'LWP_evol'+expnames[i,j]+'.pdf')
                LWP_m = np.mean(LWP_tot,axis=1)  
                
                TPcumsum = np.cumsum(np.mean(TP[:,mask],axis=1))
                #LWP_m -= (M_EVAP) #-TPcumsum)
                ax.plot(M_EVAP,'--',color=colorsbyh[i,j],label=expnames[i,j])
                ax.plot(LWP_m,color=colorsbyh[i,j],label=expnames[i,j])
                ax.set_xticks(range(0,48,6))
                ax.set_xticklabels(map(daytime,range(0,48,6)))
                ax.set_xlabel('time of day')
                ax.set_ylabel('total water path in mm')  
                ax.set_ylim([0,40])
                ax.grid(True)
                #ax.set_title(expn)
                
                LWPmax,LWPmin = np.max(LWP_m),LWP_m[-1]
                LWP_del = LWPmax-LWPmin
                LWP_tots[i,j] = LWP_del
                ncfl.close()
        ax.legend()
        figname = FIGPATH+'_LWP_evol.pdf'
        f.savefig(figname)
        
        
        f,ax=plt.subplots(1,1)
        
        if isavefig: 
            figname = FIGPATH +'/LWP_mE_evol.pdf'                
            ax.plot(np.transpose(LWP_tots/rainamounts))
            ax.legend(hrange)
            f.savefig(figname,bbox_inches='tight')
        
        
           

        
    iplt_extremes = False
    
    if iplt_extremes:
        BASE = BASETP    

        oros = expnames[:,0] # ['h250a20','h500a14','h1000a10']
        #oros = ['h250a30','h500a20','h1000a14']
        oros = ['h250a20','h500a14','h1000a10']
        f,ax = plt.subplots(1,1)
        i,j=2,2
        nbins = 30
        freq = np.zeros(nbins)
        for oro in oros:
            expn = oro
            basepath = BASE+expn+'/60_homo/'      
            print expn
            
            ensmems = filter(lambda x: x.endswith('series_TP.nc'), os.listdir(basepath))        
            i,j=np.where(expnames==oro)
            mask = R <= 30.#arange[j]        
            
            
            ipltens = True
            if ipltens:
                enslist = ensmems[:-1]
            else:
                enslist = ensmems[:1]

            
            
            for testmem in enslist:
                print testmem
                ncfl = Dataset(basepath+testmem)
                tp = ncfl.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo]            
                tp = tp[:,mask]

       
                #freq,bins,temp=plt.hist(tp.flatten(),range=[0,60],bins=60,alpha=0.5,normed=False,
                #        color=colorsbyh[i,j],label=expn,histtype='step',log=1)
                        
                memfreq,bins=np.histogram(tp.flatten(),range=[0,70],bins=nbins,normed=True)
                freq += memfreq
                
        
                ncfl.close()
            
            binscent = 0.5*(bins[:-1]+bins[1:])
            ax.semilogy(binscent,freq/np.sum(freq),alpha=0.5,label=oro)
        
        ax.set_xlabel('mm/h')
        ax.set_ylabel('frequency')
        ax.grid()
        ax.legend()
        if isavefig: f.savefig(FIGPATH+'allmem_extremes.pdf',dpi=mydpi,bbox_inches='tight')
                
                


        
    iplt_E_LWP_scaling =0
    if iplt_E_LWP_scaling:
        evalrad = 30.
        Esignal = np.zeros(expnames.shape)
        
        f,ax=plt.subplots()
        for i,j in [(0,4),(1,3),(3,2)]:
            expn = expnames[i,j]
            
            print expn
            basepath = BASE3D + expn+'/60_homo/'
            fname = basepath+filter(lambda x: x.startswith('LWP_ens'),os.listdir(basepath))[0]
            ncfl = Dataset(fname)
            
            ncflens = Dataset(basepath+ENSMEAN)
            TP = ncflens.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo]*2 # mm / h precip
            E = ncflens.variables['EFLUX'][:,-1,nhalo:-nhalo,nhalo:-nhalo]*1800./Lv # mm evaporation
            H = ncflens.variables['HFLUX'][:,-1,nhalo:-nhalo,nhalo:-nhalo] # mm evaporation
          
            LWP_t = np.zeros(ncfl.variables['LWP_QV'].shape,dtype=float)        
            
            #evalrad = arange[j]
            masks = (Rs <= evalrad)#arange[j])    
            mask = R <= evalrad
    
            
            
            for valn in filter(lambda x:x.startswith('LWP'),ncfl.variables.keys()):
                LWP_var = ncfl.variables[valn][:]  
                print LWP_var.shape
                LWP_t+=LWP_var
                dc=np.mean(LWP_var[:,masks],axis=1)           
                #ax.plot(dc,label=valn)
    
            dc_t = np.mean(LWP_t[:,masks],axis=1)
            dc_tp = np.mean(TP[:,mask],axis=1)
            dc_E = np.mean(E[:,mask],axis=1)
    
            ax.plot(dc_t-np.mean(dc_t),'--',color=colorsbyh[i,j],label=expn)
            ax.plot(dc_tp,'-',color=colorsbyh[i,j])#,label='P(r<30) '+expn)
            #ax.plot(np.cumsum(dc_E),'r',color=colorsbyh[i],label='TP')
            
            
           
            ax.grid(True)
            ncfl.close()
            ncflens.close()  
            
        xtcks = np.arange(0,49,6)
        ax.set_xticks(xtcks)
        ax.set_xlim([0,47])
        ax.set_xticklabels(map(daytime,xtcks))
        ax.legend()
        figname = FIGPATH + 'Fx_LWP_Tp_E.pdf'
        f.savefig(figname)

        
    iplt_dc_ = False
    if iplt_dc_:
        # diurnal cycle of prcecipitation
        # as a function of h
        
        
        mynh,myna = 5,3
        fh,axh = plt.subplots(1,mynh,sharey=True)
        tptot_vals = np.zeros([nh,na])
        #axh,axhhist = ax

        for i in range(mynh):
            for j in range(myna):
                print expnames[i,j]
                filename=BASETP+expnames[i,j]+'/60_homo/dc_TP.nc'
                ncfl = Dataset(filename,'r')
                mask = R <= 10 #arange[j]                
                tp = ncfl.variables['TOT_PREC'][:]
                ncfl.close()
                tp = tp[:,mask]
                #average over 5 days
                nppd = 48 # number of pts per day
                a,npts = tp.shape
                tp_day = np.zeros([nppd,npts])
                for ti in range(nppd):
                    tp_day[ti] = np.mean(tp[ti::nppd,:],axis=0)                  
                tp_dc = np.mean(tp_day[:],axis=1) 
                meantp = np.sum(tp_dc[12:36])/2.
                tptot_vals[i,j] = meantp
                axh[i].plot(tp_dc,color=colnamesbya[j],label=str(format(meantp,'.2f')))
            axh[i].set_xlim([20,36])
            axh[i].set_xticks(range(12,37,4))
            
                
                
                
                #filename = BASETP+expnames[i,j]+'/60_homo/meandc_TP.nc'
                #ncfl = Dataset(filename,'r')
                #mask = R <= arange[j]
                #print arange[j]
                #tp = ncfl.variables['TOT_PREC'][:]
                #tp_dc = np.mean(tp[:,mask],axis=1)
                
                #axh[i].plot(tp_dc,label=expnames[i,j],linewidth=2,color=colnamesbya[i])
                #axh[i].set_title('h='+str(hrange[i])+' m',fontsize=12)
                #ncfl.close()
            
        #axh[0].legend()
            
        figname = FIGPATH + 'dc_all_byh.pdf'
        if isavefig: fh.savefig(figname,dpi=mydpi,bbox_inches='tight')     
        assert 0   
        
        # as a function of a
        fa,axa = plt.subplots(1,na,sharey=True)
        tptot_vals = np.zeros([nh,na])
        #axa,axahist = ax
        for j in range(na):        
            for i in range(nh):
                print expnames[i,j]
                filename=BASETP+expnames[i,j]+'/60_homo/dc_TP.nc'
                ncfl = Dataset(filename,'r')
                mask = R <= arange[j]                
                tp = ncfl.variables['TOT_PREC'][48:]
                tp = tp[:,mask]
                #average over 5 days
                nppd = 48 # number of pts per day
                a,npts = tp.shape
                tp_day = np.zeros([nppd,npts])
                for ti in range(nppd):
                    tp_day[ti] = np.mean(tp[ti::nppd,:],axis=0)                  
                tp_dc = np.mean(tp_day[:],axis=1)        
                
                meantp = np.sum(tp_dc[12:36])/2.
                tptot_vals[i,j] = meantp
                axa[j].plot(tp_dc,color=colnamesbyh[i],label=str(format(meantp,'.2f')))
                ncfl.close()
            axa[j].grid(True)
            axa[j].set_title('a = '+str(arange[j])+' km')
            axa[j].set_xlim([12,42])
            axa[j].set_xticks(range(12,42,6))
            axa[j].legend()
        
        f,ax=plt.subplots(1,1)
        for h in tptot_vals:
            ax.plot(arange,h)
        
        #axa[0].legend()
            
        figname = FIGPATH + 'dc_all_bya.pdf'
        if isavefig: fa.savefig(figname)  
        

    iplt_vol_scalingabs_smallscale = False
    
    if iplt_vol_scalingabs_smallscale:
        expnames_temp= ['h1000a10',
                        'h1000a5_x4_offset_10',
                        'h1000a7_x2_offset_10',
                        'h250a10_x4_offset_10',
                        'h500a10_x2_offset_10']
        i = 0
        statvari = 'abs'
        evalrad = 30.
        nexps = len(expnames_temp)      
        nens = 10
        ir,jr = 2,0
        f,ax=plt.subplots()
        tpms = np.zeros([nexps,nens])
        for expn in expnames_temp:
            print expn
            
            basepath = BASETP+expn+'/60_homo/'
            fln = 'ensmean.nc'
            ensmemlist = filter(lambda x: x.endswith('daysum_TP.nc'), os.listdir(basepath))
            
            if statvari=='abs':
                mask = (R <= evalrad)
            else:
                mask = (R <= arange[j])
            kens=0
            for fln in ensmemlist:
                
                flp = basepath+fln
                nc = Dataset(flp)

                tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)   
                tpms[i,kens] = np.mean(tp[mask])
                nc.close()
                kens+=1
                if kens>9:
                    break
            i+=1
            
        tpensmean = np.median(tpms,axis=1)       
        V0 = 1
        tp0 = tpensmean[0]
        
        tpensmean/=tp0
        tpms/=tp0
        volume/=V0            
        
        
        scalevar = V0*np.ones(nexps)         
        
        for i in range(nexps):
                tpmean=np.median(tpms[i,:])
                
                ax.plot(scalevar[i]*np.ones(10),tpms[i,:],'s',
                        markersize=4,linewidth=3,
                        markeredgewidth=0,alpha=0.5)
                        
                ax.plot(scalevar[i],tpmean,'s',
                         markersize=8,markeredgewidth=0)

                ax.text(scalevar[i],tpmean,expnames_temp[i])
        f.savefig(FIGPATH+'manymtns_scaling.pdf',dpi=mydpi,bbox_inches='tight')
            

    iplt_draw_surfprof = False
    if iplt_draw_surfprof:
        r=np.linspace(0,100,1000)
        a=20
        h1=1/(1+(r/a)**2)**(1.5)
        h2=2**(-(r/a)**2)
        h3=np.cos(np.pi/4*(r/a)**2)**2
        h3[np.where(r>np.sqrt(2)*a)]=0.
        f,ax=plt.subplots()
        for h in [h1,h2,h3]:
            ax.plot(r,h)
        ax.legend(['witch','gauss','cos$^2$'])
        f.savefig(FIGPATH+'surftopo.pdf')

    iplt_draw_topo = False
    
    if iplt_draw_topo:
        f,ax=plt.subplots()
        x=np.arange(-256,256)
        y=np.copy(x)
        X,Y=np.meshgrid(x,y)       
        c=0
        
        def gauss(c1,c2,h,a):
            return h*np.exp(-((X-c1)**2+(Y-c2)**2)/a**2)
        
        dx=10
        h1=gauss(-dx,-dx,500,7.071)
        h2=gauss(-dx,dx,500,7.071)
        h3=gauss(dx,-dx,0,7.071)
        h4=gauss(dx,dx,0,7.071)
        h=h1+h2+h3+h4        
        ax.contourf(X,Y,h,cmap='terrain')
        f.savefig(FIGPATH+'oroexmpl.pdf')

    
    iplt_vol_scalingabs = 0
    
    if iplt_vol_scalingabs:
        evalrads = [20]#,30,40,50]
        evalradcols = {20:'r',30:'g',40:'b',50:'k'}
        #f,ax=plt.subplots(2,4)
        nens=10 #10 otherwise
        tpms = np.zeros([nh,na,nens])
        expnames_bell = ['h500a'+str(a)+'_bell' for a in [5,7,10,14,20,25,30]]
        statvari = 'abs'
        ir,jr = 3,2
        f2,ax2=plt.subplots(2,1)
        markerbyh=[1,2,3,4,8]
        for evalrad in evalrads:
            for i in range(nh):
                for j in range(na):
                    expn = expnames[i,j]#expnames_bell[j]   
                    print expn             
                    basepath = BASETP+expn +'/60_homo/'
                    fln = 'ensmean.nc'
                    ensmemlist = filter(lambda x: x.endswith('daysum_TP.nc'), os.listdir(basepath))
                    
                    if statvari=='abs':
                        mask = (R <= evalrad)
                    else:
                        strw = evalrad
                        nc = Dataset(basepath+fln)
                        tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)
                        nc.close()
                        tpstrip = tp[c-strw:c+strw,:]
                        tpdist = np.mean(tpstrip,axis=0)
                        cmov = np.argmax(tpdist)
                        print cmov
                        R=np.sqrt((X-cmov)**2+(Y-c)**2)
                        mask = R<=arange[j]

                        
                    kens=0
                    for fln in ensmemlist:
                        
                        flp = basepath+fln
                        nc = Dataset(flp)
                        
                        
                        tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/5.   
                        tpms[i,j,kens] = np.mean(tp[mask])
                        nc.close()
                        kens+=1
                        if kens>9:
                            break
                        
            tpensmean = np.median(tpms,axis=2)       
            V0 = volume[ir,jr]       
            tp0 = tpensmean[ir,jr]
            
            
            
            amatr = np.zeros((nh,na))
            for i in range(nh): amatr[i,:] = arangeval            
            
            tpensmean/=tp0
            tpms/=tp0
            volume/=V0            
            
            volume_bell = np.copy(volume)
            for i in range(nh):
                volume_bell[i,:]=volume[1,:]
            
            scalevar = volume
            
#            aranges = np.zeros(volume.shape)
#            for i in range(nh):
#                aranges[i,:] = arange
#                
#            scalevar = aranges
            
            #ax2.scatter(volume/V0,tpms/tp0,color=evalradcols[evalrad],label='r = '+str(evalrad)+' km')
            for i in range(nh):
                for j in range(na):
                    tpmin=np.min(tpms[i,j,:])
                    tpmax=np.max(tpms[i,j,:])
                    tpmedian=np.median(tpms[i,j,:])
                    tpmean=np.mean(tpms[i,j,:])
                    
                    #Minmax
#                    ax2.plot(volume[i,j]*np.ones(2),[tpmin,tpmax],'s-',
#                             color=colorsbyh[i,j],markersize=3,linewidth=3,
#                            markeredgewidth=0,alpha=0.5)
                    for axi in range(2):
                        ax2[axi].plot(scalevar[i,j]*np.ones(10),tpms[i,j,:],'s',
                                 color=colorsbyh[i,j],markersize=4,linewidth=3,
                                markeredgewidth=0,alpha=0.5)
                        ax2[axi].plot(scalevar[i,j],tpmedian,'s-',color=colorsbyh[i,j],
                                 markersize=8,markeredgewidth=0)
                    #ax2.plot(scalevar[i,j],tpmean,'o',color=colorsbyh[i,j],
                    #         markersize=8,markeredgewidth=0)


           
            for axi in range(2):
                for i in range(nh):
                    ax2[axi].plot(scalevar[i,:],tpensmean[i,:],'s',color=colorsbyh[i,j],
                             markersize=6, markeredgewidth=0,linewidth=3,label=str(hrange[i])+' m')

        
                
                ax2[axi].set_ylabel('P/P$_0$')
                ax2[axi].grid(True)
                #ax2.set_xlim([0,250])
                #ax2.set_ylim([0,7])
        ax2[0].legend(loc='lower right',frameon=False)
        ax2[1].set_xlabel('V/V$_0$')       
        ax2[1].set_xlim([0,2])
        ax2[1].set_ylim([0,2])
        
        
        circ = plt.Circle((1,1),0.5,color='black')
        for axi in range(2):
            ax2[axi].fill_between([0,2],0,2,alpha=0.2,color='grey')
            #ax2[axi].add_artist(circ)
        
        for a in ax2: adjust_spines(a,['left','bottom'])
        
        
        figname = FIGPATH+statvari+'_TP_vs_Vol.pdf'        
        if isavefig: f2.savefig(figname,dpy=mydpi,bbox_inches='tight')


    iplt_ha_MFETP_ =0
    if iplt_ha_MFETP_:
        # h,a-contourfield
        H,A=np.meshgrid(arange,hrange)        
   
        totalvals = np.zeros([nh,na])
        totalvals_E = np.copy(totalvals)
        totalvals_H = np.copy(totalvals)
        totalvals_MF = np.copy(totalvals)
        
        fmf,axmf = plt.subplots(nh,na,sharex=True,sharey=True)
        
        for i in range(nh):
            for j in range(na):
                                                
                name=expnames[i,j]
                mask = R <= 30. #arangeval[j]
                
                # PRECIPITATION
                nc = Dataset(BASETP+name+'/60_homo/ensmean.nc')
                tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)
                nc.close()
                totalvals[i,j] = np.mean(tp[mask])
             
                # EVAPORATION
                nc = Dataset(BASE3D+name+'/60_homo/ensmean_day_3d_d0d5.nc')
                E = np.mean(nc.variables['EFLUX'][:,-1,:],axis=0)
                H = np.mean(nc.variables['HFLUX'][:,-1,:],axis=0)
                nc.close()
                totalvals_E[i,j] = np.mean(E[mask])
                totalvals_H[i,j] = np.mean(H[mask])
                print expnames[i,j]              
                
                # MOISTURE FLUX
                nc = Dataset(BASE3D+name+'/60_homo/'+CIRCMEAN)
                ra = arange[j]
                U = nc.variables['Urz'][:,::-1,ra]
                RHO = nc.variables['RHOrz'][:,::-1,ra]
                QV = nc.variables['QVrz'][:,::-1,ra]
                               
                nc.close()
                
                #thickness of layers
                h_at_a=hrange[i]/2.
                vc = get_vertcoord(h_at_a)[::-1]
                vc_t = vc[1:]-vc[0:-1]
                vc_half = 0.5*(vc[1:]+vc[0:-1])
                
                # moisture flux calculation
                mflux = RHO*QV*U
                
                # integration of vertical profile
                mflux_int = np.sum(mflux[12:24,:],axis=0)*1800 # time mean
                
                mflux_integrand = vc_t*mflux_int
                klev = 50
                totmflux = np.sum(mflux_integrand[:klev])
                totalvals_MF[i,j]=totmflux #arange[j]**2
                
                axmf[i,j].plot(mflux_integrand,vc_half,label=str(totmflux))
                axmf[i,j].legend()
                # the integral    
                
                
          
        axmf[0,1].set_ylim([0,3500])      
        figname = FIGPATH+'mf_matrix.pdf'
        fmf.savefig(figname)
                       
          
        f,ax=plt.subplots()
        cf=ax.imshow(totalvals,cmap=viridis_inv,interpolation='nearest')

        ax.set_xticks(range(na))
        ax.set_yticks(range(nh))

        ax.set_yticklabels(['h'+str(h) for h in hrange])
        ax.set_xticklabels(['a5','a7','a10','a14','a20','a25','a30'])
        
        add_customcolorbar(f,cf,'mm',pos='right')
        if isavefig: f.savefig(FIGPATH+'hvsa_tp.pdf',dpi=mydpi,bbox_inches='tight')
        
        
        
        # T vs A vs P
        f,ax=plt.subplots(1,1)
        for i in range(nh):
            ax.plot(arange,totalvals[i,:],'s-',color=colnamesbyh[i],label=hrange[i])

        ax.legend()
        f.savefig(FIGPATH+'tp_vs_h_vs_a.pdf')
        
        
        
        
        
        # varcombination
        varnames= {'E':totalvals_E,'MF':totalvals_MF,'H':totalvals_H,
                    'TP':totalvals,'Vol':volume,'Sl':slope}
                    
        varlabs = {'E':'E/E$_0$','MF':'MF/MF$_0$','H':'H/H$_0$',
                    'TP':'P/P$_0$','Vol':'V/V$_0$','Sl':'S/S$_0$'}
        
        pn='TP'
        sn='MF'
        for pn in varnames:
            if pn=='Vol' or pn=='Sl':
                continue
            for sn in varnames:
                if pn==sn or sn in ['TP','Sl']:
                    continue
                print pn,sn
                f,ax=plt.subplots(1,1)
                plotvar = varnames[pn]
                scalevar = varnames[sn]

                for i in range(nh):
                    ax.plot(scalevar[i,:]/scalevar[2,3],plotvar[i,:]/plotvar[2,3],'s-', # volume[i,:]/volume[2,3]
                            markersize=8,label=str(hrange[i])+' m',markeredgewidth=0,
                            color=colnamesbyh[i],alpha=0.6,linewidth=2)
                    for j in range(na): 
                        ax.text(scalevar[i,j]/scalevar[2,3],plotvar[i,j]/plotvar[2,3],
                                arange[j],fontsize=6,color='k',horizontalalignment='center',
                                verticalalignment='center')
                            
                
                varxlab=varlabs[sn]
                varylab=varlabs[pn]
                ax.set_xlabel(varxlab)
                ax.set_ylabel(varylab) 
                #ax.set_ylim([80,120])
                ax.grid(True)
                #ax.set_xticks(range(11))
                #ax.set_yticks(np.arange(11)/2.)

                ax.legend(loc='lower right',frameon=False)
                
                figname=FIGPATH+pn+'_vs_'+sn+'.pdf'
                f.savefig(figname)
                del f
                
  
    iplt_calc_meanTP=False
    if iplt_calc_meanTP:
        flname = '/net/o3/hymet/adeli/project_A/256x256_7Kkmnowind_1km/postprocessing/composites/flat/60_homo/ensmean_day_d1d5.nc'
        nc=Dataset(flname)
        tpmean = np.mean(nc.variables['TOT_PREC'][:],axis=0)
        print np.mean(tpmean)
        plt.contourf(tpmean)
        nc.close()
      
      
    iplt_windsim_TP = 0
    if iplt_windsim_TP:
        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
        BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/'
        BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/'
        
        BASEU5 = '/net/o3/hymet/adeli/project_B2/512x512_7KkmU5_1km/postprocessing/composites/'
        BASEU10 = '/net/o3/hymet/adeli/project_B2/512x512_7KkmU10_1km/postprocessing/composites/'
        ipltvecs = False
        evalrad = 30. 
        i_evalmethodnew = True
        
        BASEps=[BASEnow, BASEw05,BASEw1,BASEU5,BASEU10]
        
        nsets=5#len(BASEps)
        simsetname = ['U0','U1','U2','Uni5','Uni10'][:nsets]
        
        
        arange = [20,20]
        oron = ['h500a'+str(a) for a in arange]#filter(lambda x: x.startswith('h500'), os.listdir(BASE))
        oron_bell = ['h500a'+str(a)+'_bell' for a in arange]
        oron_cos2 = ['h500a'+str(a)+'_cos2' for a in arange]
        #orons = np.array([oron,oron_bell,oron_cos2,oron,oron])
        orons = np.array([oron,oron,oron,oron,oron])
        na = len(oron)
        levs=np.linspace(5.,80,7)
        f,ax=plt.subplots(nsets,na)
        rainsignal = np.zeros([nsets,na])
        
        myarange = [0,2,4]
        for i in range(nsets):
            for j in range(na):
                BASE = BASEps[i]
                expn = orons[i,j]
                fln = BASE+'/TOT_PREC/'+expn+'/60_homo/ensmean.nc'
                nc = Dataset(fln)
                TPfld = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/2.
                nc.close()
                if i==0:
                    fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/interp_ensmean_day_3d_d0d5.nc'
                    
                    nc = Dataset(fln)
                    print fln
                    U=nc.variables['Uinterp'][20,26,nhalo:-nhalo,nhalo:-nhalo]
                    V=nc.variables['Vinterp'][20,26,nhalo:-nhalo,nhalo:-nhalo]
                    X=nc.variables['lon'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                    Y=nc.variables['lat'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                    print X.shape
        
                rainsignal[i,j]=np.mean(TPfld[np.where(TPfld>10.)])
                if i==0 and j==0:
                    mintp=np.mean(TPfld)
                    #levs[0] = 5.
    
                h,a=hrange[i],arange[j]
    
    
                if i_evalmethodnew:
                    print 0

                else:
                    c1=c+arange[j]*0.5
                    c2=c
                    R = np.sqrt((X-c1)**2+(Y-c2)**2)
                    mask= R<30. #arange[j]
                    
                print 
                #rainamnt = np.sum(TPfld[:,mask],axis=0)
                #rainsignal[i,j] = np.mean(TPfld[mask])            
                
                cf = ax[i,j].contourf(X,Y,TPfld,cmap=viridis_inv,levels=levs,extend='both')
                #ax[j].quiver(X[::5,::5],Y[::5,::5],U[::5,::5],V[::5,::5])
                st = int(a/5)
               
                titlename = simsetname[i]#+' dommean='+str(np.mean(TPfld[:]))
                print titlename
                ax[i,j].set_title(titlename,fontsize=7,fontweight='bold')
                
                
                c = 259
                a = arange[j]
                add_circle(ax[i,j],'h500m','60_homo',c,a)
                # add analysis evaluation circle
                if i_evalmethodnew:
                    strw = evalrad/3.
                    centr=256
                    tpstrip = TPfld[centr-strw:centr+strw,:]
                    tpdist = np.mean(tpstrip,axis=0)
                    cmov = np.argmax(tpdist)
                    cx = X[centr,cmov]
                    cy = Y[centr,cmov]
                    circoff = plt.Circle((cx,cy),30.,color='red',fill=False,linestyle='dashed',linewidth=1.5)
                    ax[i,j].add_artist(circoff)    
                    ax[i,j].plot(cx,cy,'x',markersize=9,color='red')
                    
                else:
                    if i in [0,1]:
                        circoff = plt.Circle((c,c),30.,color='blue',fill=False,linestyle='dashed',linewidth=1.5)
                        ax[i,j].add_artist(circoff)
                    if i == 2:
                        circoff = plt.Circle((c+0.5*a,c),30.,color='blue',fill=False,linestyle='dashed',linewidth=1.5)
                        ax[i,j].add_artist(circoff)
                    if i == 4:
                        circoff = plt.Circle((c+60.,c),30.,color='blue',fill=False,linestyle='dashed',linewidth=1.5)
                        ax[i,j].add_artist(circoff)
                wsf=3

                
                
                tcks = np.arange(259-wsf*a,259+wsf*a+1.,a)
                ax[i,j].set_xticks(tcks)
                ax[i,j].set_yticks(tcks)
                
                tckl = np.arange(-wsf,wsf+1,1)
                tckl = np.array(tckl,dtype=object)
                tckl = map(lambda x: x if x%2==1 else '',tckl)
                ax[i,j].set_xticklabels(tckl)
                ax[i,j].set_yticklabels(tckl)
               
                ax[i,j].tick_params(axis='both',which='major',labelsize=6)
                ax[i,j].grid(True)
                wsfl = 1
                ax[i,j].set_xlim([c-wsf*a,c+wsf*a])
                ax[i,j].set_ylim([c-wsf*a,c+wsf*a])
                ax[i,j].set_aspect('equal')
                #ax[i,j].tick_params(axis='both', which='major', labelsize=6)


        for j in range(na):
            ax[-1,j].set_xlabel('x/a')
        for i in range(nsets): ax[i,0].set_ylabel('y/a')
        add_customcolorbar(f,cf,'rain amount in mm',pos='right')
        figname = FIGPATH +'Fx_tpfield??_windsims.pdf'
        print figname
        f.subplots_adjust(hspace=0.25, wspace=-0.5)
        if isavefig: f.savefig(figname,dpi=mydpi,bbox_inches='tight')

    iplt_windsim_TP_selection = 0
    if iplt_windsim_TP_selection:
        """ Print only selection of hs as function of wind."""
        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
        BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/'
        BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/'
        
        BASEU5 = '/net/o3/hymet/adeli/project_B2/512x512_7KkmU5_1km/postprocessing/composites/'
        BASEU10 = '/net/o3/hymet/adeli/project_B2/512x512_7KkmU10_1km/postprocessing/composites/'        
        
        
        ipltvecs = False
        evalrad = 30. 
        i_evalmethodnew = True
        
        BASEps=[BASEnow, BASEw05,BASEw1,BASEU5, BASEU10]
        BASEps=[BASEnow, BASEw05,BASEw1]#,BASEU5, BASEU10]
        
        nsets=len(BASEps)
        simsetname = ['U0','U5','U10','Uni5','Uni10'][:nsets]
        
        
        ipltaconst = 1
        if ipltaconst:
           ipltcomb = 1
           if ipltcomb: 
               hrange = [500,1500,3000]
               #hrange = [1000,2000,3000]
           else:
               hrange = [1500]
           oron = ['h'+str(h)+'a14' for h in hrange]
           arange = [14,14,14,14,14] 
            
        else:
            arange = [14]
            oron = ['h500a'+str(a) for a in arange]

        orons = np.array([oron,oron,oron,oron,oron])

        na = len(oron)
        levs=[1,5,10,15,20,25,30,50]#np.linspace(1,30,7)
        f,ax=plt.subplots(na,nsets,figsize=(10,10))

        if na == 1:
            ax = np.array([ax,ax])
        
        
        rainsignal = np.zeros([na,nsets])

        for i in range(na):
            for j in range(nsets):
                BASE = BASEps[j]
                
                if ipltaconst:
                    expn = oron[i]
                else:
                    expn = 'h500a'+str(arange[i])
                #fln = BASE+'/TOT_PREC/'+expn+'/60_homo/ensmean.nc'
                fln = BASE +'/ALLVAR_3D/'+expn+'/60_homo/ensmean_day_3d_d0d5.nc'
                nc = Dataset(fln)
                TPfld = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)
                nc.close()
        
                #frainsignal[j,i]=np.mean(TPfld[np.where(TPfld>10.)])
                if i==0 and j==0:
                    mintp=np.mean(TPfld)
                    #levs[0] = 5.
    
                h,a=500,arange[i]
    
    
                if i_evalmethodnew:
                    print 0

                else:
                    c1=c+arange[i]*0.5
                    c2=c
                    R = np.sqrt((X-c1)**2+(Y-c2)**2)
                    mask= R<30. #arange[j]
                    
                print 
                
                cf = ax[i,j].contourf(X,Y,TPfld,cmap=viridis_inv,levels=levs,extend='both')

                titlename = expn+'_'+simsetname[j] #+' dommean='+str(np.mean(TPfld[:]))
                print titlename
                ax[i,j].set_title(titlename,fontsize=10,fontweight='bold')
                
                
                c = 256.5
                add_circle(ax[i,j],'h500m','60_homo',c,a)
                # add analysis evaluation circle
                if i_evalmethodnew:
                    strw = evalrad/3.
                    centr=256
                    tpstrip = TPfld[centr-strw:centr+strw,:]
                    tpdist = np.mean(tpstrip,axis=0)
                    cmov = np.argmax(tpdist)
                    cx = X[centr,cmov]
                    cy = Y[centr,cmov]
                    if j==0: # in the U0 case the maximum is found over the mtn top
                        cx,cy = 256.5,256.5
                    circoff = plt.Circle((cx,cy),30.,color='red',fill=False,linestyle='dashed',linewidth=1.5)
                    ax[i,j].add_artist(circoff)    
                    ax[i,j].plot(cx,cy,'x',markersize=9,color='red')
                    
                else:
                    if i in [0,1]:
                        circoff = plt.Circle((c,c),30.,color='blue',fill=False,linestyle='dashed',linewidth=1.5)
                        ax[i,j].add_artist(circoff)
                    if i == 2:
                        circoff = plt.Circle((c+0.5*a,c),30.,color='blue',fill=False,linestyle='dashed',linewidth=1.5)
                        ax[i,j].add_artist(circoff)
                    if i == 4:
                        circoff = plt.Circle((c+60.,c),30.,color='blue',fill=False,linestyle='dashed',linewidth=1.5)
                        ax[i,j].add_artist(circoff)
                wsf=3.5

                ipltvec = 1
                if ipltvec:
                    #fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/interp_ensmean_day_3d_d0d5.nc'
                    fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/interp_ensmean_day_3d_d0d5.nc'
                    
                    nc = Dataset(fln)
                    print fln
                    nt = 20
                    U=nc.variables['Uinterp'][nt,-6,nhalo:-nhalo,nhalo:-nhalo]#-np.mean(nc.variables['Uinterp'][20,-6,nhalo:-nhalo,nhalo:-nhalo])
                    V=nc.variables['Vinterp'][nt,-6,nhalo:-nhalo,nhalo:-nhalo]
                    W=nc.variables['Winterp'][nt,-6,nhalo:-nhalo,nhalo:-nhalo]
                    #X=nc.variables['lon'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                    #Y=nc.variables['lat'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                    print X.shape
                    st = 8. #int(a/5)
                    qv = ax[i,j].quiver(X[::st,::st],Y[::st,::st],U[::st,::st],V[::st,::st],
                                        units='inches',scale=10, color='grey',pivot='mid',angles='xy',
                                        headwidth=4,width=0.02)
                    #ax[i,j].contour(X,Y,W,cmap='coolwarm')
                
                
                tcks = np.arange(c-10*a,c+10*a+1.,a)
                ax[i,j].set_xticks(tcks)
                ax[i,j].set_yticks(tcks)
                
                tckl = np.arange(-10,10+1,1)
                tckl = np.array(tckl,dtype=object)
                tckl = map(lambda x: x if x%2==1 else '',tckl)
                ax[i,j].set_xticklabels(tckl)
                ax[i,j].set_yticklabels(tckl)
               
                ax[i,j].tick_params(axis='both',which='major',labelsize=6)
                ax[i,j].grid(True)
                wsfl = 1
                
                cxnew = cx #(c+cmov)/2
                ax[i,j].set_xlim([cxnew-wsf*a,cxnew+wsf*a])
                ax[i,j].set_ylim([c-wsf*a,c+wsf*a])
                ax[i,j].set_aspect('equal')
                #ax[i,j].tick_params(axis='both', which='major', labelsize=6)

        #ax[0,0].quiverkey(qv,0.5,1,2,r'$1 \frac{m}{s}$',labelsep=0.05,labelpos='N')
        
        for i in range(nsets):
            ax[-1,i].set_xlabel('x/a')
        for i in range(na): ax[i,0].set_ylabel('y/a')
        add_customcolorbar(f,cf,'rain amount in mm/d',pos='right',orient='vertical')
        
        figname = FIGPATH +'Fx_tpfield_windsims.pdf'
        if ipltaconst: figname = FIGPATH +'Fx_tpfield_windsims_aconst.pdf'

        if hrange[-1]==3000.:  figname = FIGPATH +'rev1_Fx_tpfield_windsims.pdf'
        print figname
        f.subplots_adjust(hspace=+0.5, wspace=+0.1)
        #f.tight_layout()       
        if ipltvec: ax[1,0].quiverkey(qv,-1.3,0.5,2,r'$1 \frac{m}{s}$',labelsep=0.05,labelpos='N')
        f.savefig(figname,dpi=mydpi,bbox_inches='tight')            



    
    
    
#    iplt_volscaling_windsims_OLD = 0
#    if iplt_volscaling_windsims_OLD:
#
#        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
#        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
#        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
#        BASEw2 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind2_new_1km/postprocessing/composites/'
#        BASEw4 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind4_new_1km/postprocessing/composites/'
#        cols=['blue','green','red','black','brown']
#        nsets = 3
#        BASEps=[BASEnow, BASEw05, BASEw1, BASEw2,BASEw4][:nsets]
#        nsets = len(BASEps)
#        
#        simsetname = ['500m_gauss','500m_cos2','multipeaks~V(h1000a10)',
#                      '500m_wind05 (moved window)','500m_wind2 (moved window)'][:nsets]
#        simsetname = ['U0','U5','U10','U20','U40'][:nsets]
#        
#        c=256
#        i_evalmethodnew = True
#        
#        expnames_bell = ['h500a'+str(a)+'_bell' for a in [5,7,10,14,20,25,30]]
#        expnames_cos2 = ['h500a'+str(a)+'_cos2' for a in [5,7,10,14,20,25,30]]
#        expnames_multi =['h1000a10','h1000a5_x4_offset_10','h1000a7_x2_offset_10','h250a10_x4_offset_10',
#                         'h500a10_x2_offset_10','h500a10_x2_offset_10','h500a10_x2_offset_10']
#
#
#        plth3000 = False
#        
#        pltallhs = True
#        
#        if plth3000:
#            hrange = [250,500,750,1000,1500,3000]
#            
#        if pltallhs:
#            hrange = [250,500,750,1000,1500,2000,4000]
#        else:
#            hrange = [250,500,750,1000,1500]
#        expnames_h500 = ['h500a'+str(idx) for idx in arange]
#        expnames_a14 = ['h'+str(h)+'a14' for h in hrange]
#
#        evalrad = 30.
#        #arange = [5,7,10,14,20,25,30]
#
#
#        ipltaconst = 1   
#       
#        voltest=np.zeros(len(hrange))
#        
#        voltest[0] = volume[0,3] # volume of h250a14 mountain
#        
#        
#        voltest = np.array([voltest[0]*h/250. for h in hrange])
#        
#        
#
#        f,ax2=plt.subplots(1,2,sharex=True,sharey=True,figsize=(10,4))        
#        for ipltaconst in range(2):
#            if ipltaconst:
#                    orons = expnames_a14
#            else:
#                    orons = expnames_h500
#                    
#        
#            oronames = [orons,orons,orons,orons,orons]
#    
#            volh500 = np.zeros(len(arange))
#    
#            if 0 in arange: 
#                volh500[0] = 0
#                volh500[1:] = volume[1,:]
#            else:
#                volh500 = volume[1,:]
#                
#            
#                
#            vola14 = volume[:,3]
#            if plth3000 or pltallhs: vola14=voltest
#         
#            custom_volume = vola14 if ipltaconst else volh500
#            
#            
#            vscaling = volh500[1]    
#    
#    
#            
#            na = len(orons)        
#            nens = 5
#            
#            ipltens = True            
#            
#            
#            
#            ax = ax2[ipltaconst]
#            if ipltens:
#                rainsignal = np.zeros([nsets,na,nens]) # exp x a X ens
#            else:
#                rainsignal = np.zeros([nsets,na])
#                
#            rainsignal_neg = np.copy(rainsignal) # away from the 
#                
#            
#            statvari = 'abs'
#    
#            for i in range(nsets):
#                for j in range(na):
#    
#                    BASE = BASEps[i] #BASEps[0] if i in [0,2,3] else BASEps[1]
#                    if arange[j]==0:
#                        expn = 'h500a5'
#                    else:
#                        expn = oronames[i][j]
#                    basepath= BASE+'/TOT_PREC/'+expn+'/60_homo/'
#                    #print basepath
#                    print expn
#                    
#                    if i_evalmethodnew:
#                        pass
#                    else:
#                        raise RuntimeError('not implemented')
#                        if i in [0,1,2,3]:
#                            mask = R < evalrad
#                        if i in [4]:
#                            c1=c+arange[j]*0.5
#                            c2=c
#                            Radv = np.sqrt((X-c1)**2+(Y-c2)**2)
#                            mask = Radv<evalrad #Radv < evalrad #arange[j]  
#                        if i in [5]:
#                            c1=c+50.#arange[j]*0.5
#                            c2=c
#                            Radv = np.sqrt((X-c1)**2+(Y-c2)**2)
#                            mask = Radv<evalrad #Radv < evalrad #arange[j] 
#              
#                    if ipltens:
#                        ensmemlist = filter(lambda x: x.endswith('daysum_TP.nc'), os.listdir(basepath))
#                        kens=0
#                        for fln in ensmemlist:
#                            
#                            flp = basepath+fln
#                            nc = Dataset(flp)
#                            
#                            
#                            tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)   
#                            
#                            if i_evalmethodnew: 
#                                strw = evalrad
#                                tpstrip = tp[c-strw:c+strw,:]
#                                tpdist = np.mean(tpstrip,axis=0)
#                                cmov = np.argmax(tpdist)
#                                if j==0:
#                                    print "Assumption on prec max made"
#                                    cmov=c
#                                R=np.sqrt((X-cmov)**2+(Y-c)**2)
#                                mask = R<=evalrad
#                                maskneg = R>evalrad
#                                if 0:#j==3:     
#                                    print simsetname[i],expn,cmov,c
#                                    ft,axt=plt.subplots()
#                                    axt.contour(mask,color='r')
#                                    axt.plot(cmov,c,'x')
#                                    axt.imshow(tp,#levels=[3,10,20,30,50,100],extend='both',
#                                                 cmap=viridis_inv,interpolation='nearest')
#                                
#                                if arange[j]==0:
#                                    mask = X<200.
#                            
#                            
#                            rainsignal[i,j,kens] = np.mean(tp[mask])
#                            rainsignal_neg[i,j,kens] = np.mean(tp[maskneg])
#                            
#                            nc.close()
#                            kens+=1
#                            if kens>=nens:
#                                break
#    
#                    else:    
#                        fln = BASE+'/TOT_PREC/'+expn+'ensmean.nc'
#                        nc = Dataset(fln)
#                        TPfld = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/2.
#                        nc.close()
#                    
#                        h,a=1,arange[j]
#                        
#        
#                            
#                        rainsignal[i,j]=np.mean(TPfld[mask])
#               
#                    
#            
#                          
#                         
#                          
#            pi = np.pi
#            
#            #volume_custom = np.array([volume[1,:],volume[1,:]*np.log(2),
#                                      #volume[3,2]*np.ones(len(expnames_multi)),volume[1,:],volume[1,:]])
#            
#    
#            rainsignalensm = np.mean(rainsignal,axis=2)
#            rainsignalnegensm = np.mean(rainsignal_neg,axis=2)
#            
#            if ipltaconst:
#                # take the result for h500a14
#                rscaling = rainsignalensm[0,1]
#                vscaling = vola14[1]
#                
#            else:
#                # for the 
#                rscaling = rainsignalensm[0,3]
#                vscaling = volh500[3]
#            
#            
#            #rscaling = rainsignalnegensm
#            
#            #
#            # Plotting environment
#            #
#            
#            
#            ax.hold(True)
#            for i in range(nsets):
#                ax.plot(custom_volume/vscaling,rainsignalensm[i,:]/rscaling,'s-',
#                        markeredgewidth=0,markersize=10, color=cols[i],label=simsetname[i],linewidth=1.5)
#                #ax.plot(volume_custom[i,:]/vscaling,rainsignalnegensm[i,:]/rscaling,'--',
#                #        markeredgewidth=0,markersize=8, color=cols[i],label=simsetname[i])
#                
#            if nens>1:
#                for i in range(nsets):
#                    for j in range(na):
#                        ax.plot(custom_volume[j]*np.ones(nens)/vscaling,rainsignal[i,j,:]/rscaling,'s',
#                                   markeredgewidth=0,markersize=5,color=cols[i])
#            
#            #ax.set_xticks([0,1,2,3,4,5])
#            #ax.set_yticks([0,1,2,3,4])
#            ax.plot(8,4,'white')     
#            ax.plot(0.,0.,'white') 
#            ax.legend(loc='lower right',frameon=False)
#            
#            xlenwi=7.
#            ax.plot(xlenwi,3.4,color='white')
#            
#            ax.set_ylabel('P/P$_0$')
#            ax.set_xlabel('V/V$_0$')
#            adjust_spines(ax,['left','bottom'])
#            
#            addit = '_aconst_' if ipltaconst else ''
#        
#        if plth3000:
#            figname = FIGPATH +'rev1_vscaling_VARa_VARh_new.pdf'
#        else:
#            
#            figname = FIGPATH +'vscaling_VARa_VARh.pdf'
#        print figname
#            
#        #add twinx
#        axvara = ax2[0].twiny()
#        axvarh = ax2[1].twiny()
#        axvara.set_xticks(volh500/volh500[3]/xlenwi) # to scale to the range
#        axvara.set_xticklabels(arange,fontsize=8)
##        
#        
#        for a in [axvara,axvarh]: a.grid(True)
#        axvarh.set_xticks(vola14/vola14[1]/xlenwi)
#        axvarh.set_xticklabels(np.array(hrange),fontsize=8)
##        
#        axvara.set_xlabel('a in km')
#        axvarh.set_xlabel('h in m')
#        
#        ax2[0].set_title('VARa (h=500 m)\n \n \n',fontweight='bold')
#        ax2[1].set_title('VARh (a=14 km)\n \n \n',fontweight='bold')
#        
#        f.savefig(figname,dpi=mydpi,bbox_inches='tight')


    iplt_volscaling_windsims = 0 
    if iplt_volscaling_windsims:

        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEw05 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_new_1km/postprocessing/composites/'
        BASEw1 = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind_1km/postprocessing/composites/'
        cols=['blue','green','red']
         
        BASEps=[BASEnow, BASEw05, BASEw1]
        nsets = len(BASEps)
        

        simsetname = ['U0','U5','U10']
        
        c=256
        evalrad = 30.
        i_evalmethodnew = True
        
       
        nens = 10            
        ipltens = True   
        
        VARarange = [5,7,10,14,20,25,30]
        VARhrange = [250,500,750,1000,1500,2000,3000,4000]

        expnames_VARa = ['h500a'+str(idx) for idx in VARarange]
        expnames_VARh = ['h'+str(h)+'a14' for h in VARhrange]

        expnames = {'VARa':expnames_VARa,'VARh':expnames_VARh}
       
        vscal = volume[1,3] # volume of h500a14 mountain
        
        
        
        VARhvol = np.array([vscal*h/500. for h in VARhrange])/vscal
        
        
        VARarangevals = [5,5*np.sqrt(2),10,10*np.sqrt(2),20,25,30]
        VARavol = np.array([vscal*a**2/200. for a in VARarangevals])/vscal #200=14.14**2
        
        
        volumes = {'VARa':VARavol,'VARh':VARhvol}        
        


        f,ax = plt.subplots(1,2,sharey=True,figsize=(10,5))       
        axVARa,axVARh = ax
        
        axs = {'VARa':axVARa,'VARh':axVARh}
        
        
        for simsetn in ['VARa','VARh']:
            
            orons = expnames[simsetn]
            noro = len(orons)                          
            custom_volume = volumes[simsetn]

            if ipltens:
                rainsignal = np.zeros([nsets,noro,nens]) # exp x a X ens
            else:
                rainsignal = np.zeros([nsets,noro])
                
            rainsignal_neg = np.copy(rainsignal) # away from the 
                
            
            statvari = 'abs'
    
            for i in range(nsets):
                for j in range(noro):
    
                    BASE = BASEps[i] #BASEps[0] if i in [0,2,3] else BASEps[1]
                    expn = orons[j]

                    basepath= BASE+'/TOT_PREC/'+expn+'/60_homo/'
                    print expn
              
                    if ipltens:
                        ensmemlist = filter(lambda x: x.endswith('daysum_TP.nc'), os.listdir(basepath))
                        kens=0
                        for fln in ensmemlist:
                            print fln
                            flp = basepath+fln
                            nc = Dataset(flp)
                            tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)   
                            
                            if i_evalmethodnew: 
                                strw = evalrad
                                tpstrip = tp[c-strw:c+strw,:]
                                tpdist = np.mean(tpstrip,axis=0)
                                cmov = np.argmax(tpdist)
                                if i==0:
                                    #print "Assumption on prec max made"
                                    #set equal to c if no background wind
                                    cmov=c
                                R=np.sqrt((X-cmov)**2+(Y-c)**2)
                                mask = R<=evalrad
                                maskneg = R>evalrad
#                                if 0:#j==3:     
#                                    print simsetname[i],expn,cmov,c
#                                    ft,axt=plt.subplots()
#                                    axt.contour(mask,color='r')
#                                    axt.plot(cmov,c,'x')
#                                    axt.imshow(tp,#levels=[3,10,20,30,50,100],extend='both',
#                                                 cmap=viridis_inv,interpolation='nearest')                    
                            rainsignal[i,j,kens] = np.mean(tp[mask])
                            rainsignal_neg[i,j,kens] = np.mean(tp[maskneg])
                            
                            nc.close()
                            kens+=1
                            if kens>=nens:
                                break                     
            pi = np.pi            
            # take the ensemble mean
            rainsignalensm = np.mean(rainsignal,axis=2)
            rainsignalnegensm = np.mean(rainsignal_neg,axis=2)
            
            # define the rscaling (h500a14)
            #
            # for VARh: index 2
            # for VARa: index 3
            print simsetn
            if simsetn=='VARa':
                rscaling = rainsignalensm[0,3]
            else:
                #simsetn=='VARh':
                rscaling = rainsignalensm[0,1]
                
            rainsignalensm/=rscaling
            rainsignal/=rscaling
           
            ax = axs[simsetn]
            ax.hold(True)
            for i in range(nsets):
                ax.plot(custom_volume,rainsignalensm[i,:],'s-',
                        markeredgewidth=0,markersize=10, color=cols[i],label=simsetname[i],linewidth=1.5)
                #ax.plot(volume_custom[i,:]/vscaling,rainsignalnegensm[i,:]/rscaling,'--',
                #        markeredgewidth=0,markersize=8, color=cols[i],label=simsetname[i])
                
            if nens>1:
                for i in range(nsets):
                    for j in range(noro):
                        ax.plot(custom_volume[j]*np.ones(nens),rainsignal[i,j,:],'s',
                                   markeredgewidth=0,markersize=5,color=cols[i])
            

                 
            
            ax.legend(loc='lower right',frameon=False)

            ax.set_ylabel('P/P$_0$')
            ax.set_xlabel('V/V$_0$')
            
            adjust_spines(ax,['left','bottom'])
           
            
             # to scale to the range
            ax.plot(custom_volume,np.ones(noro),alpha=0)
            axtop = ax.twiny()
            #axtop.set_xticks(custom_volume)

            if simsetn=='VARa':
                axtop.plot([0,8],[0,0],alpha=0)
                axtop.set_xticklabels(VARarange,fontsize=6)
                axtop.set_xticks(custom_volume)
                print custom_volume,simsetn
                #axtop.set_xticklabels(VARarange,fontsize=8)
                axtop.set_xlabel('a in km')
                axtop.set_title('VARa (h=500 m)\n \n \n',fontweight='bold')
    
            else:
                axtop.set_xticks(custom_volume)
                axtop.set_xticklabels(VARhrange,fontsize=6)
                axtop.set_xlabel('h in m')
                axtop.set_title('VARh (a=14 km)\n \n \n',fontweight='bold')
    
            axtop.grid(True)
            #axh.set_xticks(vola14/vola14[1]/xlenwi)
            #axh.set_xticklabels(np.array(hrange),fontsize=8)
            ax.plot(0.,0.,'white') 
            ax.plot(8,4,'white')
        
        f.subplots_adjust(wspace=0.3)
        figname = FIGPATH +'rev1_vscaling_VARa_VARh_new.pdf'
        f.savefig(figname,dpi=mydpi,bbox_inches='tight')

    iplt_volscaling_geom = 0
    if iplt_volscaling_geom:

        BASEnow = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_1km/postprocessing/composites/'
        BASEnow_NOCP = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmnowind_NOCP_1km/postprocessing/composites/'

        cols=['blue','red','g','lightblue','darkblue','black','brown','red']
        marker = ['s','^','>']
        #BASEps=[BASEnow, BASEnow,BASEnow,BASEw,BASEw2]
        BASEps=[BASEnow, BASEnow,BASEnow,BASEnow_NOCP][:3]
        nsets = len(BASEps)
        
        simsetname = ['Gauss','cos$^2$','multiple\npeaks','h500_NOCP']

        
        c=256.5
        i_evalmethodnew = True
        
        expnames_bell = ['h500a'+str(a)+'_bell1.5' for a in [5,7,10,14,20,25,30]]
        expnames_cos2 = ['h500a'+str(a)+'_cos2' for a in [5,7,10,14,20,25,30]]
        expnames_multi =['h1000a10','h1000a5_x4_offset_10','h1000a7_x2_offset_10','h250a10_x4_offset_10',
                         'h500a10_x2_offset_10','h500a10_x2_offset_10','h500a10_x2_offset_10']

        evalrad = 30.
        orons = ['h500a'+str(idx) for idx in arange]#filter(lambda x: x.startswith('h500'), os.listdir(BASE))
        oronames = [orons,expnames_cos2,expnames_multi,orons]
        na = len(orons)        
        nens = 10
        
        ipltens = True

        
        
        if ipltens:
            rainsignal = np.zeros([nsets,na,nens]) # exp x a X ens
        else:
            rainsignal = np.zeros([nsets,na])
            
        rainsignal_neg = np.copy(rainsignal) # away from the 
            
        
        statvari = 'abs'

        for i in range(nsets):
            for j in range(na):
                print 1
                BASE = BASEps[i] #BASEps[0] if i in [0,2,3] else BASEps[1]
                expn = oronames[i][j]
                basepath= BASE+'/TOT_PREC/'+expn+'/60_homo/'
                #print basepath
                
                
                if i_evalmethodnew:
                    pass
                else:
                    if i in [0,1,2,3]:
                        mask = R < evalrad
                    if i in [4]:
                        c1=c+arange[j]*0.5
                        c2=c
                        Radv = np.sqrt((X-c1)**2+(Y-c2)**2)
                        mask = Radv<evalrad #Radv < evalrad #arange[j]  
                    if i in [5]:
                        c1=c+50.#arange[j]*0.5
                        c2=c
                        Radv = np.sqrt((X-c1)**2+(Y-c2)**2)
                        mask = Radv<evalrad #Radv < evalrad #arange[j] 
          
                if ipltens:
                    ensmemlist = filter(lambda x: x.endswith('daysum_TP.nc'), os.listdir(basepath))
                    #print ensmemlist
                    
#                    if statvari=='abs':
#                        mask = mask  #(R <= evalrad)
#                    else:
#                        assert 0
#                        mask = (R <= arange[j])
                    kens=0
                    for fln in ensmemlist:
                        
                        flp = basepath+fln
                        nc = Dataset(flp)
                        
                        
                        tp = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/5   
                        
                        if i_evalmethodnew: 
                            strw = evalrad
                            tpstrip = tp[c-strw:c+strw,:]
                            tpdist = np.mean(tpstrip,axis=0)
                            cmov = np.argmax(tpdist)
                            R=np.sqrt((X-cmov)**2+(Y-c)**2)
                            mask = R<=evalrad
                            #maskneg = np.logical_and(R>evalrad,R<evalrad*2.)
                            maskneg = R>evalrad
                            if 0:#j==3:     
                                print simsetname[i],expn,cmov,c
                                ft,axt=plt.subplots()
                                axt.contour(mask,color='r')
                                axt.plot(cmov,c,'x')
                                axt.set_title(simsetname[i]+expn)
                                axt.imshow(tp,#levels=[3,10,20,30,50,100],extend='both',
                                             cmap=viridis_inv,interpolation='nearest')
                            
                           
                        # calcualted the evaluation radius   
                        
                        rainsignal[i,j,kens] = np.mean(tp[mask])
                        rainsignal_neg[i,j,kens] = np.mean(tp[maskneg])
                        
                        nc.close()
                        kens+=1
                        if kens>=nens:
                            break

                else:    
                    fln = BASE+'/TOT_PREC/'+expn+'ensmean.nc'
                    nc = Dataset(fln)
                    TPfld = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/2.
                    nc.close()
                
                    h,a=1,arange[j]
                    
    
                        
                    rainsignal[i,j]=np.mean(TPfld[mask])
           
   
        pi = np.pi
        
        volume_custom = np.array([volume[1,:],volume[1,:]*np.log(2),
                                  volume[3,2]*np.ones(len(expnames_multi))])
       
        rainsignalensm = np.mean(rainsignal,axis=2)
        rainsignalnegensm = np.mean(rainsignal_neg,axis=2)
        
        vscaling = 1000*10**2#volume_custom[0,1]        
        rscaling = rainsignalensm[0,3]
        
        #rscaling = rainsignalnegensm
        
        #
        # Plotting environment
        #
        
        f,ax=plt.subplots(figsize=(4,4))
        ax.hold(True)
        labeln = ['Gauss','cos$^2$','multiple peaks']
        for i in range(nsets):
            ax.plot(volume_custom[i,:]/vscaling,rainsignalensm[i,:]/rscaling,marker=marker[i],
                    markeredgewidth=0,markersize=8, color=cols[i],label=simsetname[i],linewidth=1.5)
            xp,yp = volume_custom[i,i+2]/vscaling,rainsignalensm[i,i+2]/rscaling
            #ax.text(xp,yp*0.9,labeln[i],color=cols[i],fontsize=12,fontweight='bold')
            
            #ax.plot(volume_custom[i,:]/vscaling,rainsignalnegensm[i,:]/rscaling,'--',
            #        markeredgewidth=0,markersize=8, color=cols[i],label=simsetname[i])
            
        if nens>1:
            for i in range(nsets):
                for j in range(na):
                    ax.plot(volume_custom[i,j]*np.ones(nens)/vscaling,rainsignal[i,j,:]/rscaling,'s',marker=marker[i],
                               markeredgewidth=0,markersize=4,color=cols[i])
                             
        

        ax.legend(loc='upper left',frameon=False)
        ax.set_ylabel('P/P$_0$')
        ax.set_xlabel('V/V$_0$')
        adjust_spines(ax,['left','bottom'])
        ax.grid(True)
        figname = FIGPATH +'geom_sims_scaling.pdf'
        if isavefig: f.savefig(figname,dpi=mydpi,bbox_inches='tight')

    iplt_elliptic_mtn = 0
    if iplt_elliptic_mtn:
        import matplotlib.patches as ptch
        BASEw = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_1km/postprocessing/composites/TOT_PREC/'
        evalrad=30
        oroh = 'h1000ax'
        exps = filter(lambda x: x.startswith(oroh),os.listdir(BASEw))
        
        nexps = len(exps)
        f,ax=plt.subplots(2,3,sharex=True,sharey=True)
        mask=R<evalrad   
        c = 256
        dr = 80
        
        expgeom = {'h500ax10ay20':(10,20),
                 'h500ax10ay40':(10,40),
                 'h500ax10ay60':(10,60),
                 'h500ax20ay10':(20,10),
                 'h500ax40ay10':(40,10),
                 'h500ax60ay10':(60,10),
                 'h1000ax10ay20':(10,20),
                 'h1000ax10ay40':(10,40),
                 'h1000ax10ay60':(10,60),
                 'h1000ax20ay10':(20,10),
                 'h1000ax40ay10':(40,10),
                 'h1000ax60ay10':(60,10)}
        
        levs=[1,5,10,15,20,25,30] #np.linspace(1.,100,8)
        
        tcks = np.arange(c-c/2,c+c/2,c/4)
        for i in range(nexps):
            fln = BASEw+exps[i]+'/60_homo/ensmean.nc'
            nc = Dataset(fln)
            
            tpfld = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/5.
            cf=ax[i/3,i%3].contourf(tpfld,cmap=viridis_inv,levels=levs,extend='both')
            tp = np.mean(tpfld[mask])
            nc.close()
            h,w = expgeom[exps[i]]
            ell = ptch.Ellipse((256,256),h*2.,w*2.,alpha=1,color='k',fill=False)
            ax[i/3,i%3].add_artist(ell)
            ax[i/3,i%3].set_aspect('equal')
            ax[i/3,i%3].set_title(exps[i]+' tp(r<'+str(evalrad)+')='+str(tp),fontsize=8)
            #if i <3:            
            #    ax[i/3,i%3].set_title(exps[i],fontsize=12,fontweight='bold')
            

            ax[i/3,i%3].set_aspect('equal')
            ax[i/3,i%3].grid(True)
            ax[i/3,i%3].set_xticks(tcks)
            ax[i/3,i%3].set_yticks(tcks)
            ax[i/3,i%3].set_xlim([c-dr,c+dr])
            ax[i/3,i%3].set_ylim([c-dr,c+dr])
            ax[i/3,i%3].set_xlabel('x in km')
            ax[i/3,i%3].set_ylabel('y in km')
            
            
            strw = evalrad/3.
            centr=256
            tpstrip = tpfld[centr-strw:centr+strw,:]
            tpdist = np.mean(tpstrip,axis=0)
            cmov = np.argmax(tpdist)
            cx = X[centr,cmov]
            cy = Y[centr,cmov]
            circoff = plt.Circle((cx,cy),evalrad,color='red',fill=False,linestyle='dashed',linewidth=1.5)
            #ax[i/3,i%3].add_artist(circoff)    
            #ax[i/3,i%3].plot(cx,cy,'x',markersize=9,color='red')

        
        add_customcolorbar(f,cf,'mm',pos='right')
        figname=FIGPATH+'Fx_elliptic_mtn_'+oroh+'.pdf'
        if isavefig: f.savefig(figname,dpi=mydpi,bbox_inches='tight')
           
           
           
    iplt_elliptic_mtn_sel = 0
    if iplt_elliptic_mtn_sel:
        import matplotlib.patches as ptch
        BASEw = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_1km/postprocessing/composites/TOT_PREC/'
        evalrad = 30.
        oroh = 'h1000ax'
        exps = filter(lambda x: x.startswith(oroh),os.listdir(BASEw))
        
        nexps = len(exps)
        f,ax=plt.subplots(1,2,sharex=True,sharey=True,figsize=(8,4))
         
        c = 256
        
        
        expgeom = {'h1000ax10ay20':(10,20),
                 'h1000ax20ay10':(20,10)}
        
        levs=[1,5,10,15,20,25,30] #np.linspace(1.,100,8)
        
        tckw = 20
        tcks = np.arange(c-2*tckw,c+2*tckw+1,tckw)
        dr = 2.5*tckw
        for i in range(2):
            expn = expgeom.keys()[i]
            print expn
            fln = BASEw+expn+'/60_homo/ensmean.nc'
            nc = Dataset(fln)
            
            tpfld = np.sum(nc.variables['TOT_PREC'][:,nhalo:-nhalo,nhalo:-nhalo],axis=0)/5.
            cf=ax[i].contourf(tpfld,cmap=viridis_inv,levels=levs,extend='both')
            
            nc.close()
            h,w = expgeom[expn]
            ell = ptch.Ellipse((256,256),h*2.,w*2.,alpha=1,color='k',fill=False)
            ax[i].add_artist(ell)
            
            ax[i].grid(True)
            ax[i].set_xticks(tcks)
            ax[i].set_yticks(tcks)
            ax[i].set_xticklabels(tcks-256)
            ax[i].set_yticklabels(tcks-256)
            ax[i].set_xlim([c-dr,c+dr])
            ax[i].set_ylim([c-dr,c+dr])
            ax[i].set_xlabel('x in km')
            if i==0: ax[i].set_ylabel('y in km')
            strw = evalrad/3.
            centr=256
            tpstrip = tpfld[centr-strw:centr+strw,:]
            tpdist = np.mean(tpstrip,axis=0)
            cmov = np.argmax(tpdist)
            cx = X[centr,cmov]
            cy = Y[centr,cmov]-1
            
            
            R = np.sqrt((X-cx)**2+(Y-cy)**2)
            mask = R<evalrad
            tp = np.mean(tpfld[mask])
            tpstr = '%.1f'%tp
            if 1:
                #fln = BASE+'/ALLVAR_3D/'+expn+'/60_homo/interp_ensmean_day_3d_d0d5.nc'
                BASEtmp = '/net/o3/hymet/adeli/project_B2/512x512_7Kkmwind05_1km/postprocessing/composites/ALLVAR_3D/'
                fln = BASEtmp+expn+'/60_homo/interp_ensmean_day_3d_d0d5.nc'
                
                nc = Dataset(fln)
                print fln
                U=nc.variables['Uinterp'][20,-6,nhalo:-nhalo,nhalo:-nhalo]#-np.mean(nc.variables['Uinterp'][20,-6,nhalo:-nhalo,nhalo:-nhalo])
                V=nc.variables['Vinterp'][20,-6,nhalo:-nhalo,nhalo:-nhalo]
                
                #X=nc.variables['lon'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                #Y=nc.variables['lat'][nhalo:-nhalo,nhalo:-nhalo]*CONV
                print X.shape
                st = 6. #int(a/5)
                Q = ax[i].quiver(X[::st,::st],Y[::st,::st],U[::st,::st],V[::st,::st],
                                    units='inches',scale=10, color='grey',pivot='mid',angles='xy',
                                    headwidth=4,width=0.02)    
                if i==0: ax[0].quiverkey(Q,1.08,0.5,2,r'$1 \frac{m}{s}$',labelsep=0.05,
                     labelpos='N')
            
            #ax[i].set_title(tpstr+' mm/d',fontsize=12,fontweight='bold')            
            circoff = plt.Circle((cx,cy),evalrad,color='red',fill=False,linestyle='dashed',linewidth=1.5)
            ax[i].add_artist(circoff)    
            ax[i].plot(cx,cy,'x',markersize=9,color='red')
            ax[i].set_aspect('equal')
        
        add_customcolorbar(f,cf,'rain amount mm/d',pos='right',orient='vertical')
        ax[i].set_aspect('equal')
        figname=FIGPATH+'Fx_elliptic_mtn_'+oroh+'.pdf'
        f.savefig(figname,dpi=mydpi,bbox_inches='tight')

    
    iplt_rdepend_TP = False
    if iplt_rdepend_TP:
        rranges = np.arange(0,100,1)
        rvals = (rranges[1:]+rranges[0:-1])/2.
        tpvals = np.copy(rvals)
        npoints = np.copy(tpvals)
        
        scalevals = np.zeros([5,7])
        
        # separate by h and a
        fa,axa=plt.subplots(1,na,sharey=True)
        fh,axh=plt.subplots(1,nh,sharey=True)
        
        for i in range(nh):
            for j in range(na):
                expn = expnames[i,j]
                hcol,acol = colorsbyh[i,j],colorsbya[i,j]
                print expn
                basepath = BASETP+expn+'/60_homo/'

                #print fln
                
                #iterate over ens?
                
                fln = basepath+'ensmean.nc'
                ncfl = Dataset(fln)
                var = np.mean(ncfl.variables['TOT_PREC'][1:,:],axis=0)
                ncfl.close()
                for ri in range(len(rranges)-1):
                    rhi = rranges[ri+1]
                    rlo = rranges[ri]
                    mask = np.logical_and(R>=rlo,R<rhi)
                    npoints[ri] = len(var[mask])
                    tpvals[ri] = np.mean(var[mask])

                # print the mean
                tpvalth = 1.5
                
                scalevals[i,j] = np.min(np.where(tpvals<tpvalth))


                                 
                axh[i].plot(rvals,tpvals,color=hcol,linewidth=2,label=expn)
                axa[j].plot(rvals,tpvals,color=acol,linewidth=2,label=expn)
                #axh[i].legend()
                axh[i].grid()
                #axa[j].legend()
                axa[j].grid()
                
                # check the ensmembers
                ipltensmem = False
                if ipltensmem:
                    ensnames = filter(lambda x: x.endswith('daysum_TP.nc'), 
                                      os.listdir(basepath))
                    for ensmem in ensnames:
                        fln = basepath +ensmem
                        ncfl=Dataset(fln)
                        var = np.mean(ncfl.variables['TOT_PREC'][1:,:],axis=0)
                        ncfl.close()
                        for ri in range(len(rranges)-1):
                            rhi = rranges[ri+1]
                            rlo = rranges[ri]
                            mask = np.logical_and(R>=rlo,R<rhi)
                            npoints[ri] = len(var[mask])
                            tpvals[ri] = np.mean(var[mask])
                    
                        #plt.plot(rvals,tpvals)              
                        axh[i].plot(rvals,tpvals,color=hcol,
                                    linewidth=1,label=expn,alpha=0.5)
                        axa[j].plot(rvals,tpvals,color=acol,
                                    linewidth=1,label=expn,alpha=0.5)
                                    
        fa.savefig(FIGPATH+'testa.pdf')
        fh.savefig(FIGPATH+'testh.pdf')
                
                # integral over ensmembers                
                
                

    		
