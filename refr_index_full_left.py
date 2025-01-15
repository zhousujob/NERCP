"""
waveproperties_mod.refr_index_full

**Description**:
_____________________________________________________________________________________________________________________

Routine to calculate the the refractive index and the wave numbers
_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________

**Inputs**:
_____________________________________________________________________________________________________________________

theta: wave normal angle in rad

S:Stix S parameter

P:Stix P parameter

R:Stix R parameter

L:Stix L parameter

D:Stix D parameter

w_wave_arg: wave frequency
______________________________________________________________________________________________________________________
_______________________________________________________________________________________________________________________

**Outputs**:
_____________________________________________________________________________________________________________________

R_tmp: refractive index
______________________________________________________________________
________________________________________________________________________________________________________________________

**Reference**:
_____________________________________________________________________________________________________________________

Kimura, I. (1966). Effects of ions on whistler-mode ray tracing. Radio Science 1, 269–283
______________________________________________________________________________________________________________________
________________________________________________________________________________________________________________________

"""

import numpy as np
from WPIT.Environment_mod import const


def refr_index_full_left(theta_arg,w_wave_arg,S_arg,D_arg,P_arg,R_arg,L_arg):
    aRef_left,aRef_right=np.nan,np.nan
    A=S_arg*np.sin(theta_arg)*np.sin(theta_arg)+P_arg*np.cos(theta_arg)*np.cos(theta_arg)
    B=R_arg*L_arg*np.sin(theta_arg)*np.sin(theta_arg)+P_arg*S_arg*(1+np.cos(theta_arg)*np.cos(theta_arg))
    C=P_arg*R_arg*L_arg
    F=np.sqrt(B*B-4*A*C)
    eta_sq_plus=(B+F)/(2*A)
    eta_sq_minus=(B-F)/(2*A)
    # if eta_sq_plus>0:
    #     ref_ind=np.sqrt(eta_sq_plus)
    # else:
    #     ref_ind=np.sqrt(eta_sq_minus)

    refsq=(B-np.sqrt(B*B-4*A*C))/(2*A)
    ref_ind=np.sqrt(refsq)

    kappa=ref_ind*w_wave_arg/const.c_light
    kappa_par=kappa*np.cos(theta_arg)
    kappa_per=kappa*np.sin(theta_arg)
    
    plr_plus=-(eta_sq_plus-S_arg)/D_arg
    plr_minus=-(eta_sq_minus-S_arg)/D_arg
    
    # if plr_plus>0 and plr_minus<0:
    #     aRef_left=np.sqrt(eta_sq_plus)
    #     aRef_right=np.sqrt(eta_sq_minus)

    # elif plr_plus<0 and plr_minus>0:
    #     aRef_left=np.sqrt(eta_sq_minus)
    #     aRef_right=np.sqrt(eta_sq_plus)
    
    
    if (plr_plus>0 and plr_minus<0) or (plr_plus>plr_minus>0) or (plr_minus<plr_plus<0):
        aRef_left=np.sqrt(eta_sq_plus)
        aRef_right=np.sqrt(eta_sq_minus)

    else:
        
        aRef_left=np.sqrt(eta_sq_minus)
        aRef_right=np.sqrt(eta_sq_plus)
        
        
        
        
        
        
        
    kappa=aRef_left*w_wave_arg/const.c_light
    kappa_par=kappa*np.cos(theta_arg)
    kappa_per=kappa*np.sin(theta_arg)
    
    
    
    
    
    
    
    return aRef_left,aRef_right,eta_sq_plus,eta_sq_minus,kappa,kappa_par,kappa_per