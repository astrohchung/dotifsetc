
# coding: utf-8

# In[ ]:


#!jupyter nbconvert --to script dotifs_util.ipynb


# In[ ]:


import numpy as np
from scipy import constants
import scipy.interpolate
import math
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
#np.set_printoptions(threshold=np.inf)
np.set_printoptions(threshold=1000)


# In[ ]:


def mag2flux(mag, zero_pt=21.1, ABwave=None):
    if ABwave != None:
        return 10.**(-0.4*(mag+2.406+5*np.log10(ABwave)))
    return 10.**(-0.4*(mag+zero_pt))


# In[ ]:


def planck(wave, temp):
    w=wave/1e8
    c1 =  3.7417749e-5  #=2*!DPI*h*c*c   
    c2 = 1.4387687    # =h*c/k
    val=c2/w/temp
    bbflux=c1/(w**5 * (math.e**val-1))*1e-8
    return bbflux


# In[ ]:


def flux2bpmag(flam, wave, filtertrans_input, filterwave=None,
               errflag=None, flam_err=None, bpmag_err=None, itpkind='linear'):

    mintrans=0.00011
    filtertrans=np.copy(filtertrans_input)
    if filterwave != None:
        itpfunc=scipy.interpolate.interp1d(filterwave, filtertrans, kind=itpkind)
        filtertrans=itpfunc(wave)
        
    filtertrans=(filtertrans >= mintrans)*filtertrans
    
    nwave=len(wave)

    if len(filtertrans) != nwave:
        errflag=1
        return
    
    constc=constants.c


    flux=np.sum(flam*filtertrans)

    
    refflam=np.ones(nwave)*3631e-23/(wave**2.)*constc*1e10
    refflux=np.sum(refflam*filtertrans)

    if flam_err != None:
        err=(np.sum((flam_err*filtertrans)**2))**0.5
        bpmag_err=-2.5/np.log(10)/flux*err

    bpmag=-2.5*np.log10(flux/refflux)
    return bpmag

