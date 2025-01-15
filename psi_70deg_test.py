# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 08:16:35 2025

@author: zhousu
"""

#import required packages
import numpy as np
import os
import sys
import matplotlib.pyplot as plt


#Define WPIT package location
current_dir =  os.path.abspath(os.path.dirname('__file__'))
fpath = os.path.abspath(current_dir + "/..")
#Add WPIT path to Python path
sys.path.append(fpath)

#import WPIT modules
import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
import WPIT.WPI_mod.EMIC_ion_mod as wpi

### Simulation parameters


L_shell=4 #L shell of simulation
lamdaeq=np.deg2rad(0) #equatorial magnetic latitude (lamda=0 deg)
Beq =env.Bmag_dipole(L_shell,lamdaeq) #equatorial magnetic field strength


wce_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.me) #equatorial electron cyclotron frequency
wcHe_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mHe) #equatorial helium cyclotron frequency
wcH_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mH) #equatorial hydrogen cyclotron frequency
wcO_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mO) #equatorial oxygen cyclotron frequency

wpe_eq=15*wce_eq #equatorial plasma frequency
w_wave=0.96*wcHe_eq #wave frequency (96% of the equatorial helium cyclotron frequency)

ne_eq=(env.const.me*env.const.epsilon0*wpe_eq*wpe_eq)/(env.const.qe*env.const.qe) #calculate equatorial electron density from equatorial plasma frequency
nH_eq=0.77*ne_eq #77% hydrogen
nHe_eq=0.2*ne_eq #20% helium
nO_eq=0.03*ne_eq #3% oxygen

Byw0_sim=3*10**(-9)  # y-component of the wave magnetic field (3nT)

aeq0_deg=50   #initial equatorial pitch angle
aeq0=np.deg2rad(aeq0_deg) #convert pitch angle to rad

Ekev0=50 #initial energy keV

# lamda0_deg=3.5 # starting electron latitude
# lamda0=np.deg2rad(lamda0_deg) #convert latitude to rad

# lamda0_deg=20.185
# lamda0=np.deg2rad(lamda0_deg) #convert latitude to rad

#calculate the latitude of magnetic mirror point
aeq0_rad=np.deg2rad(aeq0_deg)
coef=np.array([1,0,0,0,0,3*(np.sin(aeq0_rad))**4,-4*(np.sin(aeq0_rad))**4])
result=np.roots(coef)
for i in range(len(result)):
    if result[i].real>=0 and result[i].imag==0:
        kappa=result[i].real
        aMlat_rad=np.arccos(np.sqrt(kappa))
        aMlat_deg=np.rad2deg(aMlat_rad)
lamda0_deg=aMlat_deg
lamda0=np.deg2rad(lamda0_deg) #convert latitude to rad







theta0_deg=10**-5  # wave normal angle

theta0_deg=70

theta0=np.deg2rad(theta0_deg) #convert wave normal angle to rad

eta0_deg=np.linspace(-180,180,48) #initial phases of electrons, 48 test electrons are used
eta0_deg=np.linspace(0,360,48) #initial phases of electrons, 48 test electrons are used



eta0=np.deg2rad(eta0_deg) #convert initial phases to rad

m_res=1 #WPI resonance number (0=Landau resonance)

### Integration parameters ##################################################################
Tgyro=(2*np.pi)/wcH_eq #proton gyro period
t=15 #simulation duration (s) 
h=Tgyro/50 #simulation stepsize
Nsteps=int(t/h) #number of simulation steps


#using WPIT.Environment_mod.aeq2alpha routine
alpha0=env.aeq2alpha(L_shell,lamda0,aeq0)

print('\u03B1:',np.rad2deg(alpha0))


#using WPIT.Environment_mod.initial_velocity routine
upar0,uper0,ppar0,pper0,gamma0=env.initial_velocity(Ekev0,alpha0,env.const.mH)

print('upar0:',upar0,'m/s')
print('uper0:',uper0,'m/s')
print('ppar0:',ppar0,'Ns')
print('pper0:',pper0,'Ns')
print('gamma0:',gamma0)

#calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
B0 =env.Bmag_dipole(L_shell,lamda0)
#calculate electron gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
wce_0=env.omega_cyclotron(B0,env.const.qe,env.const.me)
#calculate helium gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
wcHe_0=env.omega_cyclotron(B0,env.const.qe,env.const.mHe)
#calculate hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
wcH_0=env.omega_cyclotron(B0,env.const.qe,env.const.mH)
#calculate oxygen gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
wcO_0=env.omega_cyclotron(B0,env.const.qe,env.const.mO)

#calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
S0,D0,P0,R0,L0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B0)
#calculate the refractive index and the wavenumber at ion's initial location 
# using WPIT.WaveProperties_mod.refr_index_full routine
eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0,w_wave,S0,P0,R0,L0)

#calculate resonant velocity and energy at ion's initial location using 
# WPIT.WaveProperties_mod.resonant_velocity routine
v_para_res0, v_per_res0, v_tot_res0, E_res0,gamma_res0=wave.resonant_velocity(m_res,w_wave,kappaz0,wce_0,alpha0,env.const.mH)

#calculate the gradient of hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.dwc_ds routine
dwcds0=env.dwc_ds(wcH_0,lamda0,L_shell)
#calculate the gradient of the magnetic field strength at ion's initial location using WPIT.Environment_mod.dB_ds routine
dBdz0=env.dB_ds(B0,lamda0,L_shell)
#calculate the wave component amplitudes at ion's initial location 
# using WPIT.WaveProperties_mod.wave_amplitudes_bell routine
Bxw0, Byw0, Bzw0, Exw0, Eyw0, Ezw0=wave.wave_amplitudes_bell(mu0,P0,D0,S0,Byw0_sim,theta0)
#calculate wpi parameters at electron's initial location using WPIT.WPI_mod.EMIC_ion_mod.wpi_params routine
beta0,BwR0,BwL0,EwR0,EwL0,pwR0,pwL0,wR0,wL0=wpi.wpi_params(pper0,kappax0,env.const.qi,env.const.mH,B0,Exw0,Eyw0,Bxw0,Byw0,gamma0)
beta0_sim=0

#calulcate electorn plasma frequency
wpe_0=env.omega_plasma(ne_eq,env.const.qe,env.const.me)
#calulcate helium plasma frequency
wpHe_0=env.omega_plasma(nHe_eq,env.const.qe,env.const.mHe)
#calulcate hydrogen plasma frequency
wpH_0=env.omega_plasma(nH_eq,env.const.qe,env.const.mH)
#calulcate oxygen plasma frequency
wpO_0=env.omega_plasma(nO_eq,env.const.qe,env.const.mO)
dwcds=env.dwc_ds(wcH_0,lamda0,L_shell)

#calculate initial parameters for the investigation of non-linear interactions
C0_0=wpi.nonlinear_C0(ppar0,kappaz0,m_res,gamma0,env.const.qi,env.const.mH,wcH_0,Ezw0)
C1p_0=wpi.nonlinear_C1p(pper0,ppar0,kappaz0,m_res,env.const.qi,env.const.mH,gamma0,wR0,EwR0,wcH_0)
C1m_0=wpi.nonlinear_C1m(pper0,ppar0,kappaz0,m_res,env.const.qi,env.const.mH,gamma0,wL0,EwL0,wcH_0)
thet_0,wtrsq_0=wpi.nonlinear_theta(C0_0,C1p_0,C1m_0,m_res,beta0)
dkpar_dt0=0
H_0=wpi.nonlinear_H(pper0,ppar0,kappaz0,gamma0,m_res,env.const.mH,wcH_0,dkpar_dt0,dwcds,0)
S_0=wpi.nonlinear_S(H_0,wtrsq_0)

deta_dt0=wpi.detadt(-ppar0,m_res,wcH_0,gamma0,kappaz0,env.const.mH,w_wave)


#define empty arrays to be filled during the simulation run
pperrk_su=np.zeros((len(eta0),Nsteps+1))
pparrk_su=np.zeros((len(eta0),Nsteps+1))
etark_su=np.zeros((len(eta0),Nsteps+1))
nurk_su=np.zeros((len(eta0),Nsteps+1))
lamdark_su=np.zeros((len(eta0),Nsteps+1))
timerk_su=np.zeros((len(eta0),Nsteps+1))
uperrk_su=np.zeros((len(eta0),Nsteps+1))
uparrk_su=np.zeros((len(eta0),Nsteps+1))
zetark_su=np.zeros((len(eta0),Nsteps+1))
alphark_su=np.zeros((len(eta0),Nsteps+1))
alpha2rk_su=np.zeros((len(eta0),Nsteps+1))
aeqrk_su=np.zeros((len(eta0),Nsteps+1))
aeq2rk_su=np.zeros((len(eta0),Nsteps+1))
aeq3rk_su=np.zeros((len(eta0),Nsteps+1))
Exw_out_su=np.zeros((len(eta0),Nsteps+1))
Eyw_out_su=np.zeros((len(eta0),Nsteps+1))
Ezw_out_su=np.zeros((len(eta0),Nsteps+1))
Bxw_out_su=np.zeros((len(eta0),Nsteps+1))
Byw_out_su=np.zeros((len(eta0),Nsteps+1))
Bzw_out_su=np.zeros((len(eta0),Nsteps+1))
Bw_out_su=np.zeros((len(eta0),Nsteps+1))
Ew_out_su=np.zeros((len(eta0),Nsteps+1))
vresz_out_su=np.zeros((len(eta0),Nsteps+1))
Eres_out_su=np.zeros((len(eta0),Nsteps+1))
gammares_out_su=np.zeros((len(eta0),Nsteps+1))
mu_adiabatic_out_su=np.zeros((len(eta0),Nsteps+1))
mu_out_su=np.zeros((len(eta0),Nsteps+1))
deta_dt_out_su=np.zeros((len(eta0),Nsteps+1))
B_earth_out_su=np.zeros((len(eta0),Nsteps+1))
S_stix_out_su=np.zeros((len(eta0),Nsteps+1))
D_stix_out_su=np.zeros((len(eta0),Nsteps+1))
P_stix_out_su=np.zeros((len(eta0),Nsteps+1))
R_stix_out_su=np.zeros((len(eta0),Nsteps+1))
L_stix_out_su=np.zeros((len(eta0),Nsteps+1))
kappa_out_su=np.zeros((len(eta0),Nsteps+1))
kx_out=np.zeros((len(eta0),Nsteps+1))
kz_out=np.zeros((len(eta0),Nsteps+1))
wh_out_su=np.zeros((len(eta0),Nsteps+1))
dwce_ds_out_su=np.zeros((len(eta0),Nsteps+1))
gamma_out_su=np.zeros((len(eta0),Nsteps+1))
gamma2_out_su=np.zeros((len(eta0),Nsteps+1))

C0_out=np.zeros((len(eta0),Nsteps+1))
C1p_out=np.zeros((len(eta0),Nsteps+1))
C1m_out=np.zeros((len(eta0),Nsteps+1))
thet_out=np.zeros((len(eta0),Nsteps+1))
wtrsq_out=np.zeros((len(eta0),Nsteps+1))
dkpar_dtout=np.zeros((len(eta0),Nsteps+1))
H_out=np.zeros((len(eta0),Nsteps+1))
S_out=np.zeros((len(eta0),Nsteps+1))
detadt_out=np.zeros((len(eta0),Nsteps+1))


Phi_out_su=np.zeros((len(eta0),Nsteps+1))
E_kin_su=np.zeros((len(eta0),Nsteps+1))
E_kin_out=np.zeros((len(eta0),Nsteps+1))
u_par_out_su=np.zeros((len(eta0),Nsteps+1))
u_per_out_su=np.zeros((len(eta0),Nsteps+1))

for k in range(0,len(eta0)):

    #initialize the parameters
    pperrk_su[k,0]=pper0
    pparrk_su[k,0]=-ppar0 #'-' beacause the ion move towards the south pole
    etark_su[k,0]=eta0[k]
    nurk_su[k,0]=deta_dt0
    detadt_out[k,0]=deta_dt0
    lamdark_su[k,0]=lamda0
    timerk_su[k,0]=0
    uperrk_su[k,0]=uper0
    uparrk_su[k,0]=-upar0 #'-' beacause the ion move towards the south pole
    zetark_su[k,0]=0
    alphark_su[k,0]=alpha0
    alpha2rk_su[k,0]=alpha0
    aeqrk_su[k,0]=aeq0
    aeq2rk_su[k,0]=aeq0
    aeq3rk_su[k,0]=aeq0
    Exw_out_su[k,0]=Exw0
    Eyw_out_su[k,0]=Eyw0
    Ezw_out_su[k,0]=Ezw0
    Bxw_out_su[k,0]=Bxw0
    Byw_out_su[k,0]=Byw0
    Bzw_out_su[k,0]=Bzw0
    vresz_out_su[k,0]=v_para_res0
    Eres_out_su[k,0]=E_res0
    gammares_out_su[k,0]=gamma_res0
    mu_out_su[k,0]=mu0
    S_stix_out_su[k,0]=S0
    D_stix_out_su[k,0]=D0
    P_stix_out_su[k,0]=P0
    R_stix_out_su[k,0]=R0
    L_stix_out_su[k,0]=L0
    kappa_out_su[k,0]=kappaz0
    kx_out[k,0]=kappax0
    kz_out[k,0]=kappaz0
    gamma_out_su[k,0]=gamma0
    gamma2_out_su[k,0]=gamma0
    E_kin_su[k,0]=Ekev0*1.602176487E-16
    E_kin_out[k,0]=Ekev0*1.602176487E-16
    u_par_out_su[k,0]=-upar0 #'-' beacause the ion move towards the south pole
    u_per_out_su[k,0]=uper0
    C0_out[k,0]=C0_0
    C1p_out[k,0]=C1p_0
    C1m_out[k,0]=C1m_0
    thet_out[k,0]=thet_0
    wtrsq_out[k,0]=wtrsq_0
    dkpar_dtout[k,0]=dkpar_dt0
    H_out[k,0]=H_0
    S_out[k,0]=S_0

    i=0
    
    while i<Nsteps:

#    ######################################################################################################
#    #First step of Runge Kutta
#    ######################################################################################################

        #calculate the magnetic field strength
        B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i])
        
        #calculate the Stix parameters
        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
        
        #calculate the refractive index and the wave numbers
        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)

        #calculate the Lorentz factor
        gamma_run=((E_kin_su[k,i])/(env.const.mH*env.const.c_light*env.const.c_light))+1
        
        #calculate the hydrogen cyclotron frequency
        wcH_run=env.omega_cyclotron(B_run,env.const.qe,env.const.mH)

        #calculate the gradient of the hydrogen cyclotron frequency
        dwce_ds_run=env.dwc_ds(wcH_run,lamdark_su[k,i],L_shell)
        
        #calculate the gradient of the magnetic field strength
        dBdz_run=env.dB_ds(B_run,lamdark_su[k,i],L_shell)
    
        #define the wave field
        #waves are present only to the northern hemisphere
        if (np.rad2deg(lamdark_su[k,i]))<0:
            Byw0_s=0
        else:
            Byw0_s=Byw0_sim   
    
        #calculate the amplitudes of the wave magnetic and electric field components
        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(mu_run,P_run,D_run,S_run,Byw0_s,theta0)
        #kx_run=0
        
        #calculate the wpi parameters
        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i],kx_run,env.const.qi,env.const.mH,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)

        #calculate time derivative of the fyrofrequency
        dwcdt_run=wpi.dwcdt(pparrk_su[k,i],env.const.mH,gamma_run,dwce_ds_run)
        
        #calculate time derivative of the parallel wave number
        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h
        
        #calculate the Runge-Kutta coeeficients for each differential equation
        k1=wpi.dzdt(pparrk_su[k,i],gamma_run,env.const.mH)
        l1=wpi.dlamdadt(pparrk_su[k,i],lamdark_su[k,i],gamma_run,env.const.mH,L_shell)
        m1=wpi.dppardt(pperrk_su[k,i],etark_su[k,i],gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        n1=wpi.dpperdt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        o1=wpi.detadt(pparrk_su[k,i],m_res,wcH_run,gamma_run,kz_run,env.const.mH,w_wave)
        p1=wpi.dgammadt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run)
        q1=wpi.daeqdt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],aeqrk_su[k,i],Ezw_run,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        r1=wpi.dEkdt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
        s1=wpi.dalphadt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],Ezw_run,m_res,env.const.qi,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        
#    ######################################################################################################
#    #Second step of Runge Kutta
#    ######################################################################################################

        #calculate the magnetic field strength
        B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i]+0.5*h*l1)
        
        #calculate the Stix parameters
        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
        
        #calculate the refractive index and the wave numbers
        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)

        #calculate the Lorentz factor
        gamma_run=((E_kin_su[k,i]+0.5*h*r1)/(env.const.mH*env.const.c_light*env.const.c_light))+1
        
        #calculate the hydrogen cyclotron frequency
        wcH_run=env.omega_cyclotron(B_run,env.const.qe,env.const.mH)
#         alphark_run=np.arctan(pperrk_su[k,i]+0.5*h*n1/pparrk_su[k,i]+0.5*h*m1)
    
        #calculate the gradient of the hydrogen cyclotron frequency
        dwce_ds_run=env.dwc_ds(wcH_run,lamdark_su[k,i]+0.5*h*l1,L_shell)
        
        #calculate the gradient of the magnetic field strength
        dBdz_run=env.dB_ds(B_run,lamdark_su[k,i]+0.5*h*l1,L_shell)
    
        #define the wave field
        #waves are present only to the northern hemisphere
        if (np.rad2deg(lamdark_su[k,i]+0.5*h*l1))<0:
            Byw0_s=0
        else:
            Byw0_s=Byw0_sim   
    
        #calculate the amplitudes of the wave magnetic and electric field components
        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(mu_run,P_run,D_run,S_run,Byw0_s,theta0)
        #kx_run=0
        
        #calculate the wpi parameters
        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i]+0.5*h*n1,kx_run,env.const.qi,env.const.mH,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)

        #calculate time derivative of the fyrofrequency
        dwcdt_run=wpi.dwcdt(pparrk_su[k,i]+0.5*h*m1,env.const.mH,gamma_run,dwce_ds_run)
        
        #calculate time derivative of the parallel wave number
        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h 
        
        #calculate the Runge-Kutta coeeficients for each differential equation
        k2=wpi.dzdt(pparrk_su[k,i]+0.5*h*m1,gamma_run,env.const.mH)
        l2=wpi.dlamdadt(pparrk_su[k,i]+0.5*h*m1,lamdark_su[k,i]+0.5*h*l1,gamma_run,env.const.mH,L_shell)
        m2=wpi.dppardt(pperrk_su[k,i]+0.5*h*n1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        n2=wpi.dpperdt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        o2=wpi.detadt(pparrk_su[k,i]+0.5*h*m1,m_res,wcH_run,gamma_run,kz_run,env.const.mH,w_wave)
        p2=wpi.dgammadt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run)
        q2=wpi.daeqdt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,aeqrk_su[k,i]+0.5*h*q1,Ezw_run,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        r2=wpi.dEkdt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
        s2=wpi.dalphadt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,Ezw_run,m_res,env.const.qi,pwR_run,pwL_run,beta_run,wR_run,wL_run)

#    ######################################################################################################
#    #Third step of Runge Kutta
#    ######################################################################################################


        #calculate the magnetic field strength
        B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i]+0.5*h*l2)
        
        #calculate the Stix parameters
        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
        
        #calculate the refractive index and the wave numbers
        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)

        #calculate the Lorentz factor
        gamma_run=((E_kin_su[k,i]+0.5*h*r2)/(env.const.mH*env.const.c_light*env.const.c_light))+1
        
        #calculate the hydrogen cyclotron frequency
        wcH_run=env.omega_cyclotron(B_run,env.const.qe,env.const.mH)
#         alphark_run=np.arctan(pperrk_su[k,i]+0.5*h*n2/pparrk_su[k,i]+0.5*h*m2)
    
        #calculate the gradient of the hydrogen cyclotron frequency
        dwce_ds_run=env.dwc_ds(wcH_run,lamdark_su[k,i]+0.5*h*l2,L_shell)
        
        #calculate the gradient of the magnetic field strength
        dBdz_run=env.dB_ds(B_run,lamdark_su[k,i]+0.5*h*l2,L_shell)

        #define the wave field
        #waves are present only to the northern hemisphere
        if (np.rad2deg(lamdark_su[k,i]+0.5*h*l2))<0:
            Byw0_s=0
        else:
            Byw0_s=Byw0_sim   
    
        #calculate the amplitudes of the wave magnetic and electric field components
        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(mu_run,P_run,D_run,S_run,Byw0_s,theta0)
        #kx_run=0
        
        #calculate the wpi parameters
        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i]+0.5*h*n2,kx_run,env.const.qi,env.const.mH,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)

        #calculate time derivative of the fyrofrequency
        dwcdt_run=wpi.dwcdt(pparrk_su[k,i]+0.5*h*m2,env.const.mH,gamma_run,dwce_ds_run)
        
        #calculate time derivative of the parallel wave number
        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h
        
        #calculate the Runge-Kutta coeeficients for each differential equation
        k3=wpi.dzdt(pparrk_su[k,i]+0.5*h*m2,gamma_run,env.const.mH)
        l3=wpi.dlamdadt(pparrk_su[k,i]+0.5*h*m2,lamdark_su[k,i]+0.5*h*l2,gamma_run,env.const.mH,L_shell)
        m3=wpi.dppardt(pperrk_su[k,i]+0.5*h*n2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        n3=wpi.dpperdt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        o3=wpi.detadt(pparrk_su[k,i]+0.5*h*m2,m_res,wcH_run,gamma_run,kz_run,env.const.mH,w_wave)
        p3=wpi.dgammadt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run)
        q3=wpi.daeqdt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,aeqrk_su[k,i]+0.5*h*q2,Ezw_run,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        r3=wpi.dEkdt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
        s3=wpi.dalphadt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,Ezw_run,m_res,env.const.qi,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        
#    ######################################################################################################
#    #Fourth step of Runge Kutta
#    ######################################################################################################


        #calculate the magnetic field strength
        B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i]+h*l3)
        
        #calculate the Stix parameters
        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
        
        #calculate the refractive index and the wave numbers
        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)

        #calculate the Lorentz factor
        gamma_run=((E_kin_su[k,i]+h*r3)/(env.const.mH*env.const.c_light*env.const.c_light))+1
        
        #calculate the hydrogen cyclotron frequency
        wcH_run=env.omega_cyclotron(B_run,env.const.qe,env.const.mH)
#         alphark_run=np.arctan(pperrk_su[k,i]+0.5*h*n3/pparrk_su[k,i]+0.5*h*m3)
    
        #calculate the gradient of the hydrogen cyclotron frequency
        dwce_ds_run=env.dwc_ds(wcH_run,lamdark_su[k,i]+h*l3,L_shell)
        
        #calculate the gradient of the magnetic field strength
        dBdz_run=env.dB_ds(B_run,lamdark_su[k,i]+h*l3,L_shell)

        #define the wave field
        #waves are present only to the northern hemisphere
        if (np.rad2deg(lamdark_su[k,i]+h*l3))<0:
            Byw0_s=0
        else:
            Byw0_s=Byw0_sim   
    
        #calculate the amplitudes of the wave magnetic and electric field components
        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(mu_run,P_run,D_run,S_run,Byw0_s,theta0)
        #kx_run=0
        
        #calculate the wpi parameters
        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i]+h*n3,kx_run,env.const.qi,env.const.mH,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)

        #calculate time derivative of the fyrofrequency
        dwcdt_run=wpi.dwcdt(pparrk_su[k,i]+h*m3,env.const.mH,gamma_run,dwce_ds_run)
        
        #calculate time derivative of the parallel wave number
        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h
        
        #calculate the Runge-Kutta coeeficients for each differential equation
        k4=wpi.dzdt(pparrk_su[k,i]+h*m3,gamma_run,env.const.mH)
        l4=wpi.dlamdadt(pparrk_su[k,i]+h*m3,lamdark_su[k,i]+h*l3,gamma_run,env.const.mH,L_shell)
        m4=wpi.dppardt(pperrk_su[k,i]+h*n3,etark_su[k,i]+h*o3,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        n4=wpi.dpperdt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
        o4=wpi.detadt(pparrk_su[k,i]+h*m3,m_res,wcH_run,gamma_run,kz_run,env.const.mH,w_wave)
        p4=wpi.dgammadt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run)
        q4=wpi.daeqdt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,aeqrk_su[k,i]+h*q3,Ezw_run,gamma_run,m_res,env.const.qi,env.const.mH,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        r4=wpi.dEkdt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,gamma_run,m_res,env.const.qi,env.const.mH,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
        s4=wpi.dalphadt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,Ezw_run,m_res,env.const.qi,pwR_run,pwL_run,beta_run,wR_run,wL_run)
        
#    ######################################################################################################
#    Calculate parameters at final step
#    ######################################################################################################
        
        zetark_su[k,i+1]=zetark_su[k,i]+(h/6)*(k1+2*k2+2*k3+k4)
        pparrk_su[k,i+1]=pparrk_su[k,i]+(h/6)*(m1+2*m2+2*m3+m4)
        pperrk_su[k,i+1]=pperrk_su[k,i]+(h/6)*(n1+2*n2+2*n3+n4)
        etark_su[k,i+1]=(etark_su[k,i]+(h/6)*(o1+2*o2+2*o3+o4))
        lamdark_su[k,i+1]=lamdark_su[k,i]+(h/6)*(l1+2*l2+2*l3+l4)
        alphark_su[k,i+1]=alphark_su[k,i]+(h/6)*(s1+2*s2+2*s3+s4)
        aeqrk_su[k,i+1]=aeqrk_su[k,i]+(h/6)*(q1+2*q2+2*q3+q4)
        E_kin_su[k,i+1]=E_kin_su[k,i]+(h/6)*(r1+2*r2+2*r3+r4)
        gamma_out_su[k,i+1]=gamma_out_su[k,i]+(h/6)*(p1+2*p2+2*p3+p4)
        detadt_out[k,i+1]=(1/6)*(o1+2*o2+2*o3+o4)
        u_par_out_su[k,i+1]=pparrk_su[k,i+1]/(gamma_run*env.const.mH)
        u_per_out_su[k,i+1]=pperrk_su[k,i+1]/(gamma_run*env.const.mH)  
        
#    ######################################################################################################
#    Non-linear investigation
#    ######################################################################################################
        
        C0_run=wpi.nonlinear_C0(pparrk_su[k,i+1],kz_run,m_res,gamma_run,env.const.qe,env.const.mH,wcH_run,Ezw_run)
        C1p_run=wpi.nonlinear_C1p(pperrk_su[k,i+1],pparrk_su[k,i+1],kz_run,m_res,env.const.qe,env.const.mH,gamma_run,wR_run,EwR_run,wcH_run)
        C1m_run=wpi.nonlinear_C1m(pperrk_su[k,i+1],pparrk_su[k,i+1],kz_run,m_res,env.const.qe,env.const.mH,gamma_run,wL_run,EwL_run,wcH_run)
        thet_run,wtrsq_run= wpi.nonlinear_theta(C0_run,C1p_run,C1m_run,m_res,beta_run)
        H_run=wpi.nonlinear_H(pperrk_su[k,i+1],pparrk_su[k,i+1],kz_run,gamma_run,m_res,env.const.mH,wcH_run,dkpardt_run,dwce_ds_run,0)
        S_run=wpi.nonlinear_S(H_run,wtrsq_run)   
        
        C0_out[k,i+1]=C0_run
        C1p_out[k,i+1]=C1p_run
        C1m_out[k,i+1]=C1m_run
        thet_out[k,i+1]=thet_run
        wtrsq_out[k,i+1]=wtrsq_run
        H_out[k,i+1]=H_run
        S_out[k,i+1]=S_run
        
        #output wave numbers
        kx_out[k,i+1]=kx_run
        kz_out[k,i+1]=kz_run
        
        #outpout wave fields
        Exw_out_su[k,i+1]=Exw_run
        Eyw_out_su[k,i+1]=Eyw_run
        Ezw_out_su[k,i+1]=Ezw_run
        Bxw_out_su[k,i+1]=Bxw_run
        Byw_out_su[k,i+1]=Byw_run
        Bzw_out_su[k,i+1]=Bzw_run

        #output first adiabatic inveriant
        mu_out_su[k,i+1]=mu_run
         
        #calculate simulation time
        i=i+1
        timerk_su[k,i]=timerk_su[k,i-1]+h
        kappa_out_su[k,i]=kz_run
        print(timerk_su[k,i]) 
        
        
fonts=15
eta0_sim=eta0
from matplotlib import cm

norm = cm.colors.Normalize(vmin=eta0_sim.min(), vmax=eta0_sim.max())
cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap.set_array([])
cmap2 = cm.ScalarMappable(norm=norm, cmap=cm.jet)
cmap2.set_array([])

fig, ax = plt.subplots(figsize=(8,6))
s=5


for r in range(0,len(eta0_sim)):
    ax.plot(np.rad2deg(lamdark_su[r,:-1]),np.rad2deg(aeqrk_su[r,:-1]),c=cmap.to_rgba(eta0_sim[r]))
    
ax.grid(alpha=.3)
ax.set_xlim(0,lamda0_deg)
ax.set_xlabel('Latitude (deg) \n Equatorial pitch angle as a function of magnetic latitude. \n Reproduction of Figure 3f of Su et al. (2014)',fontsize=12)
ax.set_title('Equatorially-mapped pitch-angle',fontsize=fonts)
ax.set_ylabel(r'$\alpha_{eq}$ (deg)',fontsize=fonts)

ax.set_ylim(49.7,50.3)
cbar=fig.colorbar(cmap, ticks=[-np.pi,-np.pi/2,0,np.pi/2,np.pi])
cbar.ax.set_yticklabels([r'$-\pi$',r'$-\pi/2$', '0',r'$\pi/2$', r'$-\pi$'])  # horizontal colorbar
cbar.set_label(r'$\eta_0$ (rad)', rotation=270,labelpad=15,fontsize=fonts)

#uncomment for saving the figure
# plt.tight_layout()
# plt.savefig('Figure_3f_of_Su_et_al_2014')

plt.show()

fig, ax = plt.subplots(figsize=(8,6))
s=5


for r in range(0,len(eta0_sim)):
    ax.plot(np.rad2deg(lamdark_su[r,:-1]),E_kin_su[r,:-1]/1.602176487E-16,c=cmap.to_rgba(eta0_sim[r]))

ax.grid(alpha=.3)
ax.set_xlim(0,lamda0_deg)
ax.set_xlabel('Latitude (deg) \n Total energy as a function of magnetic latitude \n Reproduction of Figure 3g of Su et al. (2014)',fontsize=12)
ax.set_ylabel('Ek (keV)',fontsize=fonts)
ax.set_title('Energy',fontsize=fonts)
ax.set_ylim(49.9,50.1)
cbar=fig.colorbar(cmap, ticks=[-np.pi,-np.pi/2,0,np.pi/2,np.pi])
cbar.ax.set_yticklabels([r'$-\pi$',r'$-\pi/2$', '0',r'$\pi/2$', r'$-\pi$'])  # horizontal colorbar
cbar.set_label(r'$\eta_0$ (rad)', rotation=270,labelpad=15,fontsize=fonts)

#uncomment for saving the figure
# plt.tight_layout()
# plt.savefig('Figure_3g_of_Su_et_al_2014')

plt.show()




dalpha=[]
dEkin=[]

from matplotlib import cm

aeq = np.ma.array(aeqrk_su, mask=np.isnan(aeqrk_su))
Ekin = np.ma.array(E_kin_su, mask=np.isnan(E_kin_su))
for r in range(0,len(eta0_sim)):
    ls=np.max(np.nonzero(aeq[r,:]))
    els=np.max(np.nonzero(Ekin[r,:]))
    dalphaf=np.rad2deg(aeq[r,ls])-aeq0_deg
    dEkinf=(Ekin[r,ls]/1.602176487E-16)-Ekev0
    dalpha.append(dalphaf)
    dEkin.append(dEkinf)








fig, ax = plt.subplots(figsize=(8,6))
s=5

colors=eta0_sim[:]

ax.grid(alpha=.3)
ax.plot(np.rad2deg(eta0_sim),dalpha)
ax.scatter(np.rad2deg(eta0_sim),dalpha,marker='o',facecolors='none', edgecolors=cmap.to_rgba(eta0_sim[:]),s=100)
ax.set_xlim(0,360)
ax.set_xlabel('$\eta_0$ (deg) \n Total pitch-angle change S \n Reproduction of Figure 4c of Su et al. (2014)',fontsize=12)
ax.set_ylabel(r'$\Delta \alpha_{eq}$ (deg)',fontsize=fonts)
ax.set_title('Total pitch-angle change',fontsize=fonts)
ax.set_ylim(-0.3,0.3)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)


#uncomment for saving the figure
plt.tight_layout()
plt.savefig('Figure_4c_of_Su_et_al_2014')



plt.show()

fig, ax = plt.subplots(figsize=(8,6))
s=5

colors=eta0_sim[:]
# for r in range(0,len(eta0)):
#     ax.plot(np.rad2deg(lamda[r,:-1]),np.rad2deg(aeq[r,:-1]),c=cmap.to_rgba(eta0[r]))
ax.grid(alpha=.3)
ax.plot(np.rad2deg(eta0_sim),dEkin)
ax.scatter(np.rad2deg(eta0_sim),dEkin,marker='o',facecolors='none', edgecolors=cmap.to_rgba(eta0_sim[:]),s=100)
ax.set_xlim(0,360)
ax.set_xlabel('$\eta_0$ (deg)\n Total pitch-angle change S \n Reproduction of Figure 4d of Su et al. (2014)',fontsize=12)
ax.set_ylabel(r'$\Delta E_{kin}}$ (keV)',fontsize=fonts)
ax.set_title('Total energy change',fontsize=fonts)
ax.set_ylim(-0.1,0.1)
ticks=np.arange(0,2*np.pi,1)
ticks=np.arange(0,2*np.pi,1)
cbar=fig.colorbar(cmap, ticks=ticks)
cbar.set_label('gyro-phase (rad)', rotation=270,labelpad=15)

#uncomment for saving the figure
plt.tight_layout()
plt.savefig('Figure_4d_of_Su_et_al_2014')



plt.show()





