
#%%
import numpy as np
from math import pi



#%%
def BulkSatFluidSub(Kdry,por,Kmin,Kfl2):
    """
%BulkSat1 Calculates the saturated bulk modulus after fluid substitution
%
% --Input--
% 1) Kdry: Dry rock bulk modulus (an array of n samples)
% 2) por:  Porosity of the rock (an array of n samples)
% 3) Kmin: Bulk modulus of the mineral grains (an array of n samples)
% 4) Kfl2: Bulk modulus of the substituting fluid (an array of n samples)
%
% --Output--
% 1) Ksat2: Saturated bulk modulus of after the fluid substitution (an
%           array of n samples)
% -------------------------------------------------------------------------
    """
    Ksat2 = Kdry + (1 - Kdry/Kmin)**2 / (por/Kfl2 + (1-por)/Kmin - Kdry/(Kmin**2))

    return Ksat2


def PatchySatWhiteDutta(Kdry,MUdry,Kmin,RHOmin,
    POR,PERM,SATfl1,SATfl2,RHOfl1,RHOfl2,Kfl1,Kfl2,VISCfl1,VISCfl2,r,fq):
    """
%PatchySatWhiteDutta 
%
% Input Arguments                                    Units
% ----------------------------------------------------------------
% 1.   KDry      = Dry bulk modulus                    (Pa)     
% 2.   MUDry     = Dry shear modulus                   (Pa)
% 3.   KMin      = Grain bulk modulus                  (Pa)
% 4.   RHOMin    = Grain Bulk density                  (kg/m^3)
% 5.   POR       = porosity                            (fraction)
% 6.   PERM      = absolute  permeability              (m^2)
% 7.   SATfl1    = Saturation of fluid 1               (fraction)
% 8.   SATfl2    = Saturation of fluid 2               (fraction)
% 9.   RHOfl1    = Density of fluid 1                  (kg/m^3)
% 10.  RHOfl2    = Density of fluid 2                  (kg/m^3)
% 11.  Kfl1      = Bulk moduli of fluid 1              (Pa)
% 12.  Kfl2      = Bulk moduli of fluid 2              (Pa)
% 13.  VISCfl1   = Viscosity of fluid 1                (Pa.s)
% 14.  VISCfl1   = Viscosity of fluid 2                (Pa.s)  
% 15.  r         = patch radius                        (m)
% 16.  fq        = frequencies to calculate properties (hz)
%
%
% Output Variables
% ------------------------------------------------------------------
% 1. Bulk      = bulk modulus of resulting composite
% 2. Vp        = p-wave velocity of resulting composite
% 3. Attn      = attenuation 
% 4. RHO       = density of resulting composite
%--------------------------------------------------------------------------

    """
    w    = 2*pi*fq # Angular Frequency
    b    = r/(SATfl1**(1/3)) # Outer Shell radius
    
#  Calculations for Complex Bulk
    K1   = BulkSatFluidSub(Kdry,POR,Kmin,Kfl1) # Bulk of full Sat. fluid 1.
    K2   = BulkSatFluidSub(Kdry,POR,Kmin,Kfl2) # Bulk of full Sat. fluid 2.
    R1   = (K1 - Kdry)*(3*K2 + 4*MUdry) \
        / ((1 - Kdry/Kmin)*(K2*(3*K1 + 4*MUdry) + 4*MUdry*(K1 - K2)*SATfl1))
    R2   = (K2 - Kdry)*(3*K1 + 4*MUdry) \
        / ((1 - Kdry/Kmin)*(K2*(3*K1 + 4*MUdry) + 4*MUdry*(K1 - K2)*SATfl1))
    KA1  = (POR/Kfl1 + (1-POR)/Kmin - Kdry/(Kmin**2))**(-1)
    KA2  = (POR/Kfl2 + (1-POR)/Kmin - Kdry/(Kmin**2))**(-1)
    KE1  = KA1*(1 - ((Kfl1*(1 - K1/Kmin)*(1 - Kdry/Kmin)) \
        / (POR*K1*(1 - Kfl1/Kmin))))
    KE2  = KA2*(1 - ((Kfl2*(1 - K2/Kmin)*(1 - Kdry/Kmin)) \
        / (POR*K2*(1 - Kfl2/Kmin))))
    a1   = np.sqrt( 1j*w*VISCfl1/(PERM*KE1) )
    a2   = np.sqrt(1j*w*VISCfl2/(PERM*KE2))
    Z1   = (r*VISCfl1*(1 - np.exp(-2*a1*r))) \
        / (PERM*((a1*r - 1) + (a1*r + 1)*np.exp(-2*a1*r)))
    Z2     = (a2*b+1) + (a2*b-1)*np.exp(2*(b-r)*a2)
    Z2     = -Z2*VISCfl2*r/PERM
    Z2     = Z2/(((a2*b+1)*(a2*r-1))-((a2*b-1)*(a2*r+1)*np.exp(2*a2*(b-r))))
    Q1   = (1 - Kdry/Kmin)*KA1/K1
    Q2   = (1 - Kdry/Kmin)*KA2/K2
    W    = (3*r*r*(R1 - R2)*(-Q1 + Q2))/((b**3)*1j*w*(Z1 + Z2))
    Kinf = (K2*(3*K1 + 4*MUdry) + 4*MUdry*(K1 - K2)*SATfl1) \
        / ((3*K1 + 4*MUdry) - 3*(K1 - K2)*SATfl1)
    Bulk = Kinf/(1 - Kinf*W) # Complex Bulk
    
    #  Calculations for P-wave Velocity and Attenuation
    RHO  = SATfl1*((1 - POR)*RHOmin + POR*RHOfl1) + \
           SATfl2*((1 -POR)*RHOmin + POR*RHOfl2)
    M    = Bulk + 4*MUdry/3
    thta = np.arctan(np.imag(M)/np.real(M))
    Vp   = np.sqrt(np.abs(M)/RHO)/np.cos(thta/2) # P-wave Velocity
    Attn = w*np.tan(thta/2)/Vp # Attenuation
    Bulk = np.abs(Bulk)

    return Bulk,Vp,Attn,RHO


# #%%
# Kdry=1
# MUdry=1
# Kmin=0.3
# RHOmin=1
# POR=0.2
# PERM=1
# SATfl1=1
# SATfl2=1
# RHOfl1=1
# RHOfl2=1
# Kfl1=1
# Kfl2=1
# VISCfl1=1
# VISCfl2=1
# r=1
# fq=1


# Bulk,Vp,Attn,RHO=PatchySatWhiteDutta(Kdry,MUdry,Kmin,RHOmin,
#     POR,PERM,SATfl1,SATfl2,RHOfl1,RHOfl2,Kfl1,Kfl2,VISCfl1,VISCfl2,r,fq)
# %%
