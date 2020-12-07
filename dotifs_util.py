#!/usr/bin/env python
# coding: utf-8

# In[15]:


#!jupyter nbconvert --to python dotifs_util.ipynb


# In[2]:


import numpy as np
from scipy import constants
import scipy.interpolate
import math
#np.set_printoptions(threshold=np.inf)
np.set_printoptions(threshold=1000)


# In[11]:


def mag2flux(mag, zero_pt=21.1, ABwave=None):
    if np.any(ABwave != None):
         return 10.**(-0.4*(mag+2.406+5*np.log10(ABwave)))
    return 10.**(-0.4*(mag+zero_pt))


# In[ ]:


def flux2mag(flux, zero_pt=21.1, ABwave=None):
    if np.any(ABwave != None):
        return np.log10(flux)/(-0.4)-2.406-5*np.log10(ABwave)
    return np.log10(flux)/(-0.4)-zero_pt


# In[4]:


def planck(wave, temp):
    w=wave/1e8
    c1 =  3.7417749e-5  #=2*!DPI*h*c*c   
    c2 = 1.4387687    # =h*c/k
    val=c2/w/temp
    bbflux=c1/(w**5 * (math.e**val-1))*1e-8
    return bbflux


# In[5]:


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


# In[ ]:




