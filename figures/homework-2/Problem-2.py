#!/usr/bin/env python
# coding: utf-8

# # Import Packages

# - First, let's run the cell below to import packages and the last cell with [Helper Functions](#helper).
# - Back from Helper Function
# <a id='helper_back'></a>

# In[135]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi
from scipy.io import loadmat


# ## (a)
# Please check the budiansly function

# ## (b)

# In[136]:


K = 36.6
mu = 45
a = 0.05
ncracks_per_volume = np.linspace(0, 1000, num=1000)
Kw = 0.01
muw = 0.0

crack_length1 = 0.1
K_dry1, mu_dry1, porosity1 = budiansly(crack_length1/2, K, mu, a, ncracks_per_volume)
K2HSu_1, K2HSl_1, U2HSu_1, U2HSl_1 = hashin_shtrikman(porosity1, Kw, muw, K, mu)

crack_length2 = 0.2
K_dry2, mu_dry2, porosity2 = budiansly(crack_length2/2, K, mu, a, ncracks_per_volume)
K2HSu_2, K2HSl_2, U2HSu_2, U2HSl_2 = hashin_shtrikman(porosity2, Kw, muw, K, mu)

fig, axs = plt.subplots(1,2, figsize=(12, 6))
axs[0].plot(ncracks_per_volume, K_dry1, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=1mm')
axs[0].plot(ncracks_per_volume, K2HSu_1, linestyle='--', linewidth=2, alpha=1, label='HS Upper')
axs[0].plot(ncracks_per_volume, K2HSl_1, linestyle='--', linewidth=2, alpha=1, label='HS Lower')
axs[0].plot(ncracks_per_volume, K_dry2, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=2mm')
axs[0].legend(loc='upper right',fontsize=10, shadow=False) 
axs[0].set_xlabel('N/Vbulk') 
axs[0].set_ylabel('Bulk modulus (GPa)')
axs[0].set_xlim(0,1000)
axs[0].set_ylim(0, )

axs[1].plot(ncracks_per_volume, mu_dry1, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=1mm')
axs[1].plot(ncracks_per_volume, U2HSu_2, linestyle='--', linewidth=2, alpha=1, label='HS Upper')
axs[1].plot(ncracks_per_volume, U2HSl_2, linestyle='--', linewidth=2, alpha=1, label='HS Lower')
axs[1].plot(ncracks_per_volume, mu_dry2, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=2mm')
axs[1].legend(loc='upper right',fontsize=10, shadow=False) 
axs[1].set_xlabel('N/Vbulk') 
axs[1].set_ylabel('Shear modulus (GPa)')
axs[1].set_xlim(0,1000)
axs[1].set_ylim(0, )
fig.savefig('./p2-b.pdf', dpi=800, format='pdf')


# ## (c)

# In[137]:


K = 36.6
mu = 45
a = 0.05
ncracks_per_volume = np.linspace(0, 20000, num=1000)
Kw = 0.01
muw = 0.0

crack_length1 = 0.1
K_dry, mu_dry, porosity = budiansly(crack_length1/2, K, mu, a, ncracks_per_volume)
K2HSu, K2HSl, U2HSu, U2HSl = hashin_shtrikman(porosity, Kw, muw, K, mu)

fig, axs = plt.subplots(1,1, figsize=(12, 6))
axs.plot(porosity, K_dry, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=1mm')
axs.plot(porosity, K2HSu, linestyle='--', linewidth=2, alpha=1, label='HS Upper')
axs.plot(porosity, K2HSl, linestyle='--', linewidth=2, alpha=1, label='HS Lower')
axs.legend(loc='upper right',fontsize=10, shadow=False) 
axs.set_xlabel('Porosity') 
axs.set_ylabel('Bulk modulus (GPa)')
axs.set_xlim(0, 0.4)
axs.set_ylim(0, )
fig.savefig('./p2-c.pdf', dpi=800, format='pdf')


# ## (d)

# In[138]:


K = 36.6
mu = 45
a = 0.05
ncracks_per_volume = np.linspace(0, 1000, num=1000)
Kw = 0.01
muw = 0.0
rho = 2.65
crack_length2 = 0.2
K_dry, mu_dry, porosity = budiansly(crack_length2/2, K, mu, a, ncracks_per_volume)
K2HSu, K2HSl, U2HSu, U2HSl = hashin_shtrikman(porosity, Kw, muw, K, mu)

vp_budiansly = np.power((K_dry+(4*mu_dry/3)) / rho , 1/2)
vp_HSu = np.power((K2HSu+(4*U2HSu/3)) / rho , 1/2)
vp_HSl = np.power((K2HSl+(4*U2HSl/3)) / rho , 1/2)

fig, axs = plt.subplots(1,1, figsize=(12, 6))
axs.plot(porosity, vp_budiansly, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=2mm')
axs.plot(porosity, vp_HSu, linestyle='--', linewidth=2, alpha=1, label='HS Upper')
axs.plot(porosity, vp_HSl, linestyle='--', linewidth=2, alpha=1, label='HS Lower')
axs.legend(loc='upper right',fontsize=10, shadow=False) 
axs.set_xlabel('Porosity') 
axs.set_ylabel('P wave velocity (km/s))')
axs.set_xlim(0, )
axs.set_ylim(0, )
fig.savefig('./p2-d.pdf', dpi=800, format='pdf')


# ## (e)

# In[139]:


Mat = loadmat("Data.mat")
Data = Mat['Data']

# Select the useful data
Data30Vp=[]; Data30Porosity=[]; Data30Fraction=[]
Data40Vp=[]; Data40Porosity=[]; Data40Fraction=[]
for i in range(0,np.size(Mat['Data'],0)):
    if Data[i,0]==30:
        Data30Vp.append(Data[i,7])
        Data30Porosity.append(Data[i,3])
        Data30Fraction.append(Data[i,4])
    elif Data[i,0]==40:
        Data40Vp.append(Data[i,7])
        Data40Porosity.append(Data[i,3])
        Data40Fraction.append(Data[i,4])

# Divide data into 2 parts: clean sandstone and shaley sandstone
CleanVp=[]; CleanPorosity=[]; CleanFraction=[]
ShaleyVp=[]; ShaleyPorosity=[]; ShaleyFraction=[]
for i in range(0, np.size(Data30Vp)):
    if Data30Fraction[i] < 0.15:
        CleanVp.append(Data30Vp[i])
        CleanPorosity.append(Data30Porosity[i])
        CleanFraction.append(Data30Fraction[i])
    else:
        ShaleyVp.append(Data30Vp[i])
        ShaleyPorosity.append(Data30Porosity[i])
        ShaleyFraction.append(Data30Fraction[i])
for i in range(0, np.size(Data40Vp)):
    if Data40Fraction[i] < 0.15:
        CleanVp.append(Data40Vp[i])
        CleanPorosity.append(Data40Porosity[i])
        CleanFraction.append(Data40Fraction[i])
    else:
        ShaleyVp.append(Data40Vp[i])
        ShaleyPorosity.append(Data40Porosity[i])
        ShaleyFraction.append(Data40Fraction[i])    
CleanVp = np.array(CleanVp)
CleanPorosity = np.array(CleanPorosity)
CleanFraction = np.array(CleanFraction)
ShaleyVp = np.array(ShaleyVp)
ShaleyPorosity = np.array(ShaleyPorosity)
ShaleyFraction = np.array(ShaleyFraction) 

# Plot
fig, axs = plt.subplots(1,1, figsize=(12, 6))
axs.plot(porosity, vp_budiansly, linestyle='--', linewidth=2, alpha=1, label='Budiansky c=2mm')
axs.plot(porosity, vp_HSu, linestyle='--', linewidth=2, alpha=1, label='HS Upper')
axs.plot(porosity, vp_HSl, linestyle='--', linewidth=2, alpha=1, label='HS Lower')
axs.plot(CleanPorosity, CleanVp,'.',markersize=15, color='red', alpha=0.5, label='Data.mat')

axs.legend(loc='upper right',fontsize=10, shadow=False) 
axs.set_xlabel('Porosity') 
axs.set_ylabel('P wave velocity (km/s))')
axs.set_xlim(0, )
axs.set_ylim(0, )
fig.savefig('./p2-e.pdf', dpi=800, format='pdf')


# ## (f)

# In[140]:


a =  [0.05, 0.10, 0.15, 0.2, 0.25, 0.3] # np.linspace(0, 1, num=1000)
K = 36.6
mu = 45
rho = 2.65
ncracks_per_volume = np.linspace(0.01, 1000, num=1000)
crack_length2 = 0.2
fig, axs = plt.subplots(1,1, figsize=(12, 6))
for i in a:
    K_dry, mu_dry, porosity = budiansly(crack_length2/2, K, mu, i, ncracks_per_volume)
    vp_budiansly = np.power((K_dry+(4*mu_dry/3)) / rho , 1/2)
    axs.plot(porosity, vp_budiansly, linestyle='--', linewidth=2, alpha=1, label='Aspect ratio = '+str(i))
    
axs.plot(CleanPorosity, CleanVp,'.',markersize=5, color='black', alpha=0.5, label='Data.mat')
axs.legend(loc='upper right',fontsize=10, shadow=False) 
axs.set_xlabel('Porosity') 
axs.set_ylabel('P wave velocity (km/s)')
axs.set_xlim(0,1)
axs.set_ylim(0, )
fig.savefig('./p2-f.pdf', dpi=800, format='pdf')


# ## Helper Functions
# <a id='helper'></a>

# In[141]:


def budiansly(crack_radius, K, mu, a, ncracks_per_volume):
    """
    Calculate the Budiansky Self-Consistency model.
    
    Parameters
    ----------
        crack_radius : float
            Crack's radius in cm.
        K : float
            Bulk modulus.
        mu : float
            Shear modulus.
        a : float 
            Aspect ratio.
        ncracks_per_volume : float or array
            The number of cracks per unit volume
            
    Returns
    -------
        k_dry : float or array
            Dry bulk modulus.
        mu_dry : float or array
            Dry shear modulus.
        porosity: float or array
            Porosity.
        
    """
    e = ncracks_per_volume*crack_radius**3  # the crack density parameter, which is defined as the number of cracks per unit volume times the crack radius cubed.
    v = ((3*K)-(2*mu))/(2*(3*K+mu))
    v_sc = v*(1-((16/9)*e))
    Kdry = K*(1-(((16/9)*((1-v_sc**2)/(1-(2*v_sc))))*e))
    mudry= mu*(1-(((32/45)*(((1-v_sc)*(5-v_sc))/(2-v_sc)))*e))
    porosity = 4*pi*a*e/3

    return Kdry, mudry, porosity


def hashin_shtrikman(F1, k1, u1, k2, u2):
    """
    Hashin-Shtrikman bounds for a mixture of two constituents.
    
    Parameters
    ----------
        F1 : array
            Array of N volume fractions of first minerals, the value is between 0 and 1 .
        k1 : float
            Bulk modulus of first minerals.
        u1 : float
            Shear modulus of first minerals.
        k2 : float
            Bulk modulus of second minerals.
        u2 : float
            Shear modulus of second minerals.
    
    Returns
    -------
        Ku : array
            Hashin-Shtrikman Bulk modulus' upper bound.
        Kl : array
            Hashin-Shtrikman Bulk modulus' lower bound
        Uu : array
            Hashin-Shtrikman Shear modulus' upper bound
        Ul : array
            Hashin-Shtrikman Shear modulus' lower bound
        
    """
    c = 4/3;
    F2 = 1-F1
    kmx = max(k1,k2)
    kmn = min(k1,k2)
    umx = max(u1,u2)
    umn = min(u1,u2)
    
    Ku = []; Kl = [];
    Uu = []; Ul = [];
    for i in range(0,len(F1)):
        F = np.array([F1[i], F2[i]])
        K = np.array([k1, k2])
        U = np.array([u1, u2])
        Kupper = 1/np.sum(F/(K+c*umx))-c*umx     # HS K upper bound
        Klow = 1/np.sum(F/(K+c*umn))-c*umn       # HS K lower bound
        etamx = umx*(9*kmx+8*umx)/(kmx+2*umx)/6
        etamn = umn*(9*kmn+8*umn)/(kmn+2*umn)/6
        Uupper = 1/np.sum(F/(U+etamx))-etamx     # HS U upper bound
        if min(U) != 0:
            Ulow = 1/np.sum(F/(U+etamn))-etamn   # HS U lower bound        
        else:
            Ulow = 0
        if F1[i] == 0:
            Ulow = 1/np.sum(F2[i]/(u2+etamn)) - etamn            
        Ku.append(Kupper); Kl.append(Klow); 
        Uu.append(Uupper); Ul.append(Ulow); 
                      
    return np.array(Ku), np.array(Kl), np.array(Uu), np.array(Ul)




# -[Turn back!](#helper_back)
