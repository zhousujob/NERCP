# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:48:16 2023

@author: zhousu
"""


# Latitudinal dependence of nonlinear interaction between electromagnetic ion cyclotron wave and terrestrial ring current ions
# Su, Zhu, Xiao, Zheng, Zhang, Liu, Shen, Wang, Wang (2014)
# Publication link: https://aip.scitation.org/doi/10.1063/1.4880036
# code source: https://github.com/stourgai/WPIT/blob/main/WPIT_tests/Latitudinal%20dependence%20of%20nonlinear%20interaction%20between%20electromagnetic%20ion%20cyclotron%20wave%20and%20terrestrial%20ring%20current%20ions.ipynb

# This script is modified from E:\mywork\python_work\examples\wave_particle_interaction\example_replot_SuZP_2014POP.py




#import required packages
import numpy as np
import os,numlib,refr_index_full_left
import sys
import matplotlib.pyplot as plt
import matplotlib,space_science
import scipy.io as sio

#import WPIT modules
import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
import WPIT.WPI_mod.EMIC_ion_mod as wpi


#Ne=210   #number density
emic_intensity=[0.1*1e-9,1*1e-9,3*1e-9,3*1e-9,3*1e-9,3*1e-9,3*1e-9,3*1e-9,3*1e-9]
freq_norm=[0.96,0.96,0.96,0.7,0.8,0.9,0.96,0.96,0.96]
wna=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,10,40,70]
par_number=[6,7,8]
Energy=np.array([50])  #initial energy keV
pitch_angle=np.arange(1,89.7,0.5)
#pitch_angle=[38]
#Nsteps=1000000




density_ratio=15
L_shell=4 #L shell of simulation
#aeq0_deg=53
lat_cutoff=8
simu_time=200  #seconds
eta0_deg=np.linspace(0,360,49) #initial phases of electrons, 48 test electrons are used
eta0=np.deg2rad(eta0_deg) #convert initial phases to rad
m_res=1 #WPI resonance number (0=Landau resonance)
k_p_par=-1
particle_mass=env.const.mH
particle_charge=env.const.qi
ion_com=[77,20,3]


#sys.exit()
save_dir=r'E:\mywork\my_drafts\nonlinear_affected_by_EMIC_parameters\data/'


LineWidth=0.8;LabelSize=12;TickWidth=1.6;TickLen=6
alpha_grid=0.36
#figure, axis lim
# xlim_s_lat=[0,30]
# ylim_s_lat=[-4,0]

xlim_eta_lat=[0,15]
xticks_eta_lat=range(0,15,2)

ylim_s_lat=[-4,0]


xlim_aeq_lat=[-20,20]
ylim_aeq_lat=[50,80]


ylim_daeq_eta=[-10,20]

ylim_Ek_lat=[48,54]
ylim_dEk_eta=[-1,3]


#panel_label=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)']
panel_label=['(a)','(f)','(k)',\
             '(b)','(g)','(l)',\
                 '(c)','(h)','(m)',\
                     '(d)','(i)','(n)',\
                         '(e)','(j)','(o)']
plt.close('all')
num_row=5
num_column=3
myfig=plt.figure(figsize=(17,15))
hf=numlib.axes_create(row=num_row, column=num_column,top=0.04,bottom=0.05,left=0.05,right=0.06,gap_row=0.04,gap_column=0.06)


lamdaeq=np.deg2rad(0) #equatorial magnetic latitude (lamda=0 deg)
Beq =env.Bmag_dipole(L_shell,lamdaeq) #equatorial magnetic field strength


wce_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.me) #equatorial electron cyclotron frequency
wcHe_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mHe) #equatorial helium cyclotron frequency
wcH_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mH) #equatorial hydrogen cyclotron frequency
wcO_eq=env.omega_cyclotron(Beq,env.const.qe,env.const.mO) #equatorial oxygen cyclotron frequency

#wpe_eq=np.sqrt(aNe*env.const.me/(Beq**2*env.const.epsilon0))*wce_eq
#w_wave=0.96*wcH_eq #wave frequency (96% of the equatorial helium cyclotron frequency)


wpe_eq=density_ratio*wce_eq #equatorial plasma frequency
ne_eq=(env.const.me*env.const.epsilon0*wpe_eq*wpe_eq)/(env.const.qe*env.const.qe) #calculate equatorial electron density from equatorial plasma frequency

nH_eq=(ion_com[0]/100)*ne_eq #77% hydrogen
nHe_eq=(ion_com[1]/100)*ne_eq #20% helium
nO_eq=(ion_com[2]/100)*ne_eq #3% oxygen


#sys.exit()
#in the following, find the cut-off frequency for left-hand EMIC waves
lat_dip=np.arange(0,89.9,0.01)


aeq_start=np.zeros((len(wna),len(Energy),len(pitch_angle),len(eta0)))*np.nan
aeq_end=aeq_start.copy()
ek_start=aeq_start.copy()
ek_end=aeq_start.copy()

transmit_time=aeq_start.copy()



#sys.exit()

for iPar in par_number:
    
    theta0_nh_deg=wna[iPar]  # wave normal angle
    theta0_nh_rad=np.deg2rad(theta0_nh_deg) #convert wave normal angle to rad

    theta0_sh_deg=180-theta0_nh_deg  # wave normal angle
    theta0_sh_rad=np.deg2rad(theta0_sh_deg) #convert wave normal angle to rad
    
    freq_ratio=freq_norm[iPar]
    wave_freq=freq_ratio*wcHe_eq

    w_wave=wave_freq
    wave_frequency=wave_freq
    
    Byw0_sim=emic_intensity[iPar]  # y-component of the wave magnetic field (3nT)
    
    
    
    for iEk in range(len(Energy)):
    
        Ekev0=Energy[iEk]
        for iAeq in range(len(pitch_angle)):
            
            aeq0_deg=pitch_angle[iAeq]
    
            #calculate the latitude of magnetic mirror point
            aeq0_rad=np.deg2rad(aeq0_deg)
            coef=np.array([1,0,0,0,0,3*(np.sin(aeq0_rad))**4,-4*(np.sin(aeq0_rad))**4])
            result=np.roots(coef)
            for iRes in range(len(result)):
                if result[iRes].real>=0 and result[iRes].imag==0:
                    kappa=result[iRes].real
                    aMlat_rad=np.arccos(np.sqrt(kappa))
                    aMlat_deg=np.rad2deg(aMlat_rad)
            lamda0_deg=aMlat_deg
            #lamda0_deg=20
            
            aRef_plus=np.ones(len(lat_dip))*np.nan  #refractive index
            aRef_minus=aRef_plus.copy()
            aRef_left=aRef_plus.copy()
            aRef_right=aRef_plus.copy()
        
            aLat=aRef_plus.copy() 
            plr_plus=aRef_plus.copy()
            plr_minus=aRef_plus.copy()
            aD_stix=aRef_plus.copy()
            #in the following, calculate the latitude where cutoff frequency occurs 
            w_wave=wave_freq

            for iLat in range(len(lat_dip)):        
                aLat_deg=lat_dip[iLat]
                aLat_rad=np.deg2rad(aLat_deg)
                aBfield=env.Bmag_dipole(L_shell,aLat_rad)
                S0,D0,P0,R0,L0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, aBfield)
                eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0_nh_rad,w_wave,S0,P0,R0,L0)
                plr_plus[iLat]=-(eta_sq_plus0-S0)/D0
                plr_minus[iLat]=-(eta_sq_minus0-S0)/D0
                ref_left0,ref_right0,eta_sq_plus0,eta_sq_minus0,kappa0,kappaz0,kappax0=refr_index_full_left.refr_index_full_left(theta0_nh_rad,w_wave,S0,D0,P0,R0,L0)
        
                aRef_left[iLat]=ref_left0
                aRef_right[iLat]=ref_right0
                aD_stix[iLat]=D0
            fp_index=np.where(np.isnan(aRef_left))
            lat_cutoff_deg=lat_dip[fp_index[0][0]-1]  #lat_cutoff_deg is the last point where the waves can propagate
            bD_stix=np.array([aD_stix[i]*aD_stix[i+1] for i in range(len(aD_stix)-1)])
            fp_D=np.where(bD_stix<0)
            lat_crossover_deg=lat_dip[fp_D[0][0]-1]  #this is the latitude of crossover frequency
            #sys.exit()
            
            
                    
            # 1. Define simulation parameters
            # Here we define all the initial parameters of the simulation in respect with the particle and the wave
                
            #aeq0=np.deg2rad(aeq0_deg) #convert pitch angle to rad
            lamda0_mirr_rad=np.deg2rad(lamda0_deg) #convert latitude to rad 
            
            ### Integration parameters ##################################################################
            Tgyro=(2*np.pi)/wcH_eq #proton gyro period
            t=simu_time #simulation duration (s) 
            h=Tgyro/100 #simulation stepsize
            Nsteps=int(t/h) #number of simulation steps
            #sys.exit()
            
            #2. Find initial electron's local pitch angle
            
            #using WPIT.Environment_mod.aeq2alpha routine
            alpha0=env.aeq2alpha(L_shell,lamda0_mirr_rad,aeq0_rad)
            
            print('\u03B1:',np.rad2deg(alpha0))
            #3. Find initial momentum, velocity and lorentz factor
            
            #using WPIT.Environment_mod.initial_velocity routine
            upar0,uper0,ppar0,pper0,gamma0=env.initial_velocity(Ekev0,alpha0,particle_mass)
            
            print('upar0:',upar0,'m/s')
            print('uper0:',uper0,'m/s')
            print('ppar0:',ppar0,'Ns')
            print('pper0:',pper0,'Ns')
            print('gamma0:',gamma0)
                
            #calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
            B0 =env.Bmag_dipole(L_shell,lamda0_mirr_rad)
            
            #define empty arrays to be filled during the simulation run
            pperrk_su=np.zeros((len(eta0),Nsteps+1))*np.nan
            pparrk_su=pperrk_su.copy()
            #sys.exit()
            etark_su=pperrk_su.copy()
            nurk_su=pperrk_su.copy()
            lamdark_su=pperrk_su.copy()
            timerk_su=pperrk_su.copy()
            uperrk_su=pperrk_su.copy()
            uparrk_su=pperrk_su.copy()
            #zetark_su=pperrk_su.copy()
            alphark_su=pperrk_su.copy()
            #alpha2rk_su=np.zeros((len(eta0),Nsteps+1))
            aeqrk_su=pperrk_su.copy()
            #sys.exit()
            #aeq2rk_su=np.zeros((len(eta0),Nsteps+1))
            #aeq3rk_su=np.zeros((len(eta0),Nsteps+1))
            # Exw_out_su=pperrk_su.copy()
            # Eyw_out_su=pperrk_su.copy()
            # Ezw_out_su=pperrk_su.copy()
            # Bxw_out_su=pperrk_su.copy()
            # Byw_out_su=pperrk_su.copy()
            # Bzw_out_su=pperrk_su.copy()
            # Bw_out_su=pperrk_su.copy()
            # Ew_out_su=pperrk_su.copy()
            #vresz_out_su=pperrk_su.copy()
            #Eres_out_su=pperrk_su.copy()
            #gammares_out_su=pperrk_su.copy()
            #mu_adiabatic_out_su=pperrk_su.copy()
            mu_out_su=pperrk_su.copy()
            #deta_dt_out_su=pperrk_su.copy()
            # B_earth_out_su=pperrk_su.copy()
            # S_stix_out_su=pperrk_su.copy()
            # D_stix_out_su=pperrk_su.copy()
            # P_stix_out_su=pperrk_su.copy()
            # R_stix_out_su=pperrk_su.copy()
            # L_stix_out_su=pperrk_su.copy()
            kappa_out_su=pperrk_su.copy()
            #kx_out=pperrk_su.copy()
            #kz_out=pperrk_su.copy()
            #wh_out_su=pperrk_su.copy()
            #dwce_ds_out_su=pperrk_su.copy()
            #gamma_out_su=pperrk_su.copy()
            #gamma2_out_su=pperrk_su.copy()
            
            # C0_out=pperrk_su.copy()
            # C1p_out=pperrk_su.copy()
            # C1m_out=pperrk_su.copy()
            # thet_out=pperrk_su.copy()
            # wtrsq_out=pperrk_su.copy()
            dkpar_dtout=pperrk_su.copy()
            #H_out=pperrk_su.copy()
            #S_out=pperrk_su.copy()
            detadt_out=pperrk_su.copy()
            
            
            #Phi_out_su=pperrk_su.copy()
            E_kin_su=pperrk_su.copy()
            E_kin_out=pperrk_su.copy()
            #u_par_out_su=pperrk_su.copy()
            #u_per_out_su=pperrk_su.copy()
            B_out=pperrk_su.copy()
            #sys.exit()
            #in the following, find the latitude of linear resonance
            lamda_arr=np.arange(lamda0_deg,0,-0.01)
            lamda_arr=np.arange(0,lamda0_deg,0.01)
        
            
            aResCon=np.zeros(len(lamda_arr))*np.nan
            bResCon=aResCon.copy()
            #sys.exit()
            for iLamda in range(len(lamda_arr)):
                aLamda_rad=np.deg2rad(lamda_arr[iLamda])
                aAlpha=env.aeq2alpha(L_shell,aLamda_rad,aeq0_rad)
                aBfield=env.Bmag_dipole(L_shell,aLamda_rad)
                aUpar,aUper,aPpar,aPper,aGamma=env.initial_velocity(Ekev0,aAlpha,particle_mass)
                aWcyc=env.omega_cyclotron(aBfield,env.const.qe,particle_mass)
                #calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
                aS,aD,aP,aR,aL=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, aBfield)
                #calculate the refractive index and the wavenumber at ion's initial location 
                # using WPIT.WaveProperties_mod.refr_index_full routine
                aEta_sq_plus,aEta_sq_minus,aMu,aKappa,aKappaz,aKappax=wave.refr_index_full(theta0_nh_rad,w_wave,aS,aP,aR,aL)
                aRef_left,aRef_right,aEta_sq_plus,aEta_sq_minus,aKappa,aKappaz,aKappax=\
                    refr_index_full_left.refr_index_full_left(theta0_nh_rad,w_wave,aS,aD,aP,aR,aL)
        
                aResCon[iLamda]=m_res*aWcyc/aGamma+(aKappaz*k_p_par*aPpar)/(aGamma*particle_mass)-w_wave
                      
            for k in range(0,len(eta0)):
                lamdark_su[k,0]=lamda0_mirr_rad
                #print('k is:'+str(k))
                #if k==2:
                #    sys.exit()
        
                #if np.rad2deg(lamdark_su[k,0])<lat_crossover_deg:
                i=0
                theta0=theta0_nh_rad
                #print(np.rad2deg(lamdark_su[k,i])) 
                print('iPar='+str(iPar))                   
                print('Ek='+str(Ekev0)+' keV,  '+'aeq0='+str(aeq0_deg)+' degree'+', eta0='+str(round(np.rad2deg(eta0[k]),2))+' degree')
                #sys.exit()
                #initialize parameters            
                #using WPIT.Environment_mod.aeq2alpha routine
                alpha0=env.aeq2alpha(L_shell,lamda0_mirr_rad,aeq0_rad)
                
                #3. Find initial momentum, velocity and lorentz factor
                #using WPIT.Environment_mod.initial_velocity routine
                upar0,uper0,ppar0,pper0,gamma0=env.initial_velocity(Ekev0,alpha0,particle_mass)
                
                #4. Calculate all the initial parameters
                
                #calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
                B0 =env.Bmag_dipole(L_shell,lamda0_mirr_rad)
                #calculate electron gyrofrequency at ion's initial location using WPIT.Environment_mod.omega_cyclotron routine
                wcyc_0=env.omega_cyclotron(B0,particle_charge,particle_mass)
        
                
                #calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
                S0,D0,P0,R0,L0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B0)
                #calculate the refractive index and the wavenumber at ion's initial location 
                # using WPIT.WaveProperties_mod.refr_index_full routine
                
                eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0,w_wave,S0,P0,R0,L0)
                #use the revised function "refr_index_full_left"
                ref_left0,ref_right0,eta_sq_plus0,eta_sq_minus0,kappa0,kappaz0,kappax0=refr_index_full_left.refr_index_full_left(theta0,w_wave,S0,D0,P0,R0,L0)
                #calculate resonant velocity and energy at ion's initial location using 
                # WPIT.WaveProperties_mod.resonant_velocity routine
                v_para_res0, v_per_res0, v_tot_res0, E_res0,gamma_res0=wave.resonant_velocity(m_res,w_wave,kappaz0,wcyc_0,alpha0,particle_mass)
                
                #calculate the gradient of hydrogen gyrofrequency at ion's initial location using WPIT.Environment_mod.dwc_ds routine
                dwcds0=env.dwc_ds(wcyc_0,lamda0_mirr_rad,L_shell)
                #calculate the gradient of the magnetic field strength at ion's initial location using WPIT.Environment_mod.dB_ds routine
                dBdz0=env.dB_ds(B0,lamda0_mirr_rad,L_shell)
                #calculate the wave component amplitudes at ion's initial location 
                # using WPIT.WaveProperties_mod.wave_amplitudes_bell routine
                Bxw0, Byw0, Bzw0, Exw0, Eyw0, Ezw0=wave.wave_amplitudes_bell(ref_left0,P0,D0,S0,Byw0_sim,theta0)
                #calculate wpi parameters at electron's initial location using WPIT.WPI_mod.EMIC_ion_mod.wpi_params routine
                beta0,BwR0,BwL0,EwR0,EwL0,pwR0,pwL0,wR0,wL0=wpi.wpi_params(pper0,kappax0,particle_charge,particle_mass,B0,Exw0,Eyw0,Bxw0,Byw0,gamma0)
                beta0_sim=0
                
                dwcds=env.dwc_ds(wcyc_0,lamda0_mirr_rad,L_shell)
                
                #calculate initial parameters for the investigation of non-linear interactions
                C0_0=wpi.nonlinear_C0(ppar0,kappaz0,m_res,gamma0,particle_charge,particle_mass,wcyc_0,Ezw0)
                C1p_0=wpi.nonlinear_C1p(pper0,ppar0,kappaz0,m_res,particle_charge,particle_mass,gamma0,wR0,EwR0,wcyc_0)
                C1m_0=wpi.nonlinear_C1m(pper0,ppar0,kappaz0,m_res,particle_charge,particle_mass,gamma0,wL0,EwL0,wcyc_0)
                thet_0,wtrsq_0=wpi.nonlinear_theta(C0_0,C1p_0,C1m_0,m_res,beta0)
                dkpar_dt0=0
                
                #revised by Su Zhou
                lamda0_rad=lamda0_mirr_rad
                lamda0_rad_2=lamda0_rad+0.00001
                #2. Find initial electron's local pitch angle
                
                #using WPIT.Environment_mod.aeq2alpha routine
                #aeq0_rad=aeq0
                alpha0_2=env.aeq2alpha(L_shell,lamda0_rad_2,aeq0_rad)
                
                #3. Find initial momentum, velocity and lorentz factor
                #using WPIT.Environment_mod.initial_velocity routine
                upar0_2,uper0_2,ppar0_2,pper0_2,gamma0_2=env.initial_velocity(Ekev0,alpha0_2,particle_mass)
                
                #4. Calculate all the initial parameters
                #calculate magnetic field strength at ion's initial location using WPIT.Environment_mod.Bmag_dipole routine
                B0_2 =env.Bmag_dipole(L_shell,lamda0_rad_2)
                
                #calculate the Stix parameters at ion's initial location using WPIT.WaveProperties_mod.stix_parameters routine
                S0_2,D0_2,P0_2,R0_2,L0_2=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B0_2)
                #calculate the refractive index and the wavenumber at ion's initial location 
                # using WPIT.WaveProperties_mod.refr_index_full routine
                
                ref_left0_2,ref_right0_2,eta_sq_plus0_2,eta_sq_minus0_2,kappa0_2,kappaz0_2,kappax0_2=\
                    refr_index_full_left.refr_index_full_left(theta0,w_wave,S0_2,D0_2,P0_2,R0_2,L0_2)
        
                
                dkdl=(kappaz0_2-kappaz0)/(lamda0_rad_2-lamda0_rad)
                dkdz=dkdl*1/(env.const.Re*L_shell*np.cos(lamda0_rad)*np.sqrt(1+3*np.sin(lamda0_rad)**2))
                dkpardt0=(ppar0/(gamma0*particle_mass))*np.cos(theta0)*dkdz
                dkpar_dt0=dkpardt0
                
                H_0=wpi.nonlinear_H(pper0,ppar0,kappaz0,gamma0,m_res,particle_mass,wcyc_0,dkpar_dt0,dwcds,0)
                S_0=wpi.nonlinear_S(H_0,wtrsq_0)
                
                
                deta_dt0=wpi.detadt(-ppar0,m_res,wcyc_0,gamma0,kappaz0,particle_mass,w_wave)
                
                if np.rad2deg(lamda0_rad)>lat_crossover_deg:
                    kappaz0,kappax0=np.nan,np.nan
                    w_wave=np.nan
                    Exw0,Eyw0,Ezw0=np.nan,np.nan,np.nan
                    Bxw0,Byw0,Bzw0=np.nan,np.nan,np.nan
                    v_para_res0,E_res0,gamma_res0=np.nan,np.nan,np.nan
                    mu0=np.nan
                    S0,D0,P0,R0,L0=np.nan,np.nan,np.nan,np.nan,np.nan
                    C0_0,C1p_0,C1m_0=np.nan,np.nan,np.nan
                    thet_0,wtrsq_0=np.nan,np.nan
                    dkpar_dt0=np.nan
                    H_0,S_0=np.nan,np.nan
                    deta_dt0=wpi.detadt(-ppar0,m_res,wcyc_0,gamma0,0,particle_mass,0)
        
                #sys.exit()
                #initialize the parameters
                pperrk_su[k,0]=pper0
                pparrk_su[k,0]=-ppar0 #'-' beacause the ion move towards the south pole
                etark_su[k,0]=eta0[k]
                nurk_su[k,0]=deta_dt0
                detadt_out[k,0]=deta_dt0
                lamdark_su[k,0]=lamda0_mirr_rad
                timerk_su[k,0]=0
                uperrk_su[k,0]=uper0
                uparrk_su[k,0]=-upar0 #'-' beacause the ion move towards the south pole
                #zetark_su[k,0]=0
                alphark_su[k,0]=alpha0
                #alpha2rk_su[k,0]=alpha0
                aeqrk_su[k,0]=aeq0_rad
                #aeq2rk_su[k,0]=aeq0
                #aeq3rk_su[k,0]=aeq0
                # Exw_out_su[k,0]=Exw0
                # Eyw_out_su[k,0]=Eyw0
                # Ezw_out_su[k,0]=Ezw0
                # Bxw_out_su[k,0]=Bxw0
                # Byw_out_su[k,0]=Byw0
                # Bzw_out_su[k,0]=Bzw0
                #vresz_out_su[k,0]=v_para_res0
                #Eres_out_su[k,0]=E_res0
                #gammares_out_su[k,0]=gamma_res0
                mu_out_su[k,0]=mu0
                # S_stix_out_su[k,0]=S0
                # D_stix_out_su[k,0]=D0
                # P_stix_out_su[k,0]=P0
                # R_stix_out_su[k,0]=R0
                # L_stix_out_su[k,0]=L0
                kappa_out_su[k,0]=kappaz0
                # kx_out[k,0]=kappax0
                # kz_out[k,0]=kappaz0
                #gamma_out_su[k,0]=gamma0
                #gamma2_out_su[k,0]=gamma0
                E_kin_su[k,0]=Ekev0*1.602176487E-16
                E_kin_out[k,0]=Ekev0*1.602176487E-16
                # u_par_out_su[k,0]=-upar0 #'-' beacause the ion move towards the south pole
                # u_per_out_su[k,0]=uper0
                # C0_out[k,0]=C0_0
                # C1p_out[k,0]=C1p_0
                # C1m_out[k,0]=C1m_0
                # thet_out[k,0]=thet_0
                # wtrsq_out[k,0]=wtrsq_0
                dkpar_dtout[k,0]=dkpar_dt0
                #H_out[k,0]=H_0
                #S_out[k,0]=S_0
                B_out[k,0]=B0
                
                #start simulation below the latitude of crossover frequency
                #while i<Nsteps:
                for i in range(Nsteps):
                    #print(i)
                    #print(np.rad2deg(lamdark_su[k,i]))
                    #B_mirr =env.Bmag_dipole(L_shell,lamda0_mirr_rad)
                    
                    
                    
                    # wc_local=particle_charge*B_out[k,i]/particle_mass
                    # Tgyro_local=(2*np.pi)/wc_local
                    # h=Tgyro_local/50
            
            #    ######################################################################################################
            #    #First step of Runge Kutta
            #    ######################################################################################################
                    #waves are present only to the northern hemisphere
                    if lamdark_su[k,i]>=0:
                        theta0=theta0_nh_rad
                    else:
                        theta0=theta0_sh_rad
        
                        
                    Byw0_s=Byw0_sim
                    if i>(Nsteps*0.9):
                        raise ValueError('the proton cannot reache the south magnetic mirror point')
                    
                    if lamdark_su[k,i]<0 and pparrk_su[k,i]>0:  #run only for the northern hemisphere
                        #print('i='+str(i))
                        break
                    
                
                 
        
                    #calculate the magnetic field strength
                    B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i])
                    gamma_run=((E_kin_su[k,i])/(particle_mass*env.const.c_light*env.const.c_light))+1
                    #calculate the hydrogen cyclotron frequency
                    wcyc_run=env.omega_cyclotron(B_run,particle_charge,particle_mass)
            
                    #calculate the gradient of the hydrogen cyclotron frequency
                    dwce_ds_run=env.dwc_ds(wcyc_run,lamdark_su[k,i],L_shell)
                    
                    #calculate the gradient of the magnetic field strength
                    dBdz_run=env.dB_ds(B_run,lamdark_su[k,i],L_shell)
        
                    #sys.exit()
                    if abs(np.rad2deg(lamdark_su[k,i]))<=lat_crossover_deg:
                        w_wave=wave_frequency
        
                        #calculate the Stix parameters
                        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
                        
                        #calculate the refractive index and the wave numbers
                        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)
                        ref_left_run,ref_right_run,eta_sq_plus,eta_sq_minus,kapp_run,kz_run,kx_run=refr_index_full_left.refr_index_full_left(theta0,w_wave,S_run,D_run,P_run,R_run,L_run)
                        
                        #calculate the amplitudes of the wave magnetic and electric field components
                        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(ref_left_run,P_run,D_run,S_run,Byw0_s,theta0)
                        #kx_run=0
                        
                        #calculate the wpi parameters
                        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i],kx_run,particle_charge,particle_mass,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)
                        
                        #calculate time derivative of the fyrofrequency
                        dwcdt_run=wpi.dwcdt(pparrk_su[k,i],particle_mass,gamma_run,dwce_ds_run)
                        
                        #calculate time derivative of the parallel wave number
                        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h
                    else:
                        Ezw_run,beta_run,wR_run,wL_run=0,0,0,0
                        pwR_run,pwL_run=0,0
                        kz_run,w_wave=0,0
                        EwR_run,EwL_run=0,0
                        kapp_run=0
                        
                    #calculate the Runge-Kutta coeeficients for each differential equation
                    k1=wpi.dzdt(pparrk_su[k,i],gamma_run,particle_mass)
                    l1=wpi.dlamdadt(pparrk_su[k,i],lamdark_su[k,i],gamma_run,particle_mass,L_shell)
                    m1=wpi.dppardt(pperrk_su[k,i],etark_su[k,i],gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    n1=wpi.dpperdt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    o1=wpi.detadt(pparrk_su[k,i],m_res,wcyc_run,gamma_run,kz_run,particle_mass,w_wave)
                    p1=wpi.dgammadt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run)
                    q1=wpi.daeqdt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],aeqrk_su[k,i],Ezw_run,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run)
                    r1=wpi.dEkdt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
                    s1=wpi.dalphadt(pperrk_su[k,i],pparrk_su[k,i],etark_su[k,i],Ezw_run,m_res,particle_charge,pwR_run,pwL_run,beta_run,wR_run,wL_run)
            
            #    ######################################################################################################
            #    #Second step of Runge Kutta
            #    ######################################################################################################
                    
                    #calculate the magnetic field strength
                    B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i]+0.5*h*l1)
                    gamma_run=((E_kin_su[k,i]+0.5*h*r1)/(particle_mass*env.const.c_light*env.const.c_light))+1
                    #calculate the hydrogen cyclotron frequency
                    wcyc_run=env.omega_cyclotron(B_run,particle_charge,particle_mass)
            #         alphark_run=np.arctan(pperrk_su[k,i]+0.5*h*n1/pparrk_su[k,i]+0.5*h*m1)
                
                    #calculate the gradient of the hydrogen cyclotron frequency
                    dwce_ds_run=env.dwc_ds(wcyc_run,lamdark_su[k,i]+0.5*h*l1,L_shell)
                    
                    #calculate the gradient of the magnetic field strength
                    dBdz_run=env.dB_ds(B_run,lamdark_su[k,i]+0.5*h*l1,L_shell)
        
                    if abs(np.rad2deg(lamdark_su[k,i]))<=lat_crossover_deg:
                        w_wave=wave_frequency
        
                        #calculate the Stix parameters
                        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
                        
                        #calculate the refractive index and the wave numbers
                        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)
                        ref_left_run,ref_right_run,eta_sq_plus,eta_sq_minus,kapp_run,kz_run,kx_run=refr_index_full_left.refr_index_full_left(theta0,w_wave,S_run,D_run,P_run,R_run,L_run)
                
                        #calculate the amplitudes of the wave magnetic and electric field components
                        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(ref_left_run,P_run,D_run,S_run,Byw0_s,theta0)
                        #kx_run=0
                        
                        #calculate the wpi parameters
                        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i]+0.5*h*n1,kx_run,particle_charge,particle_mass,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)
                
                        #calculate time derivative of the fyrofrequency
                        dwcdt_run=wpi.dwcdt(pparrk_su[k,i]+0.5*h*m1,particle_mass,gamma_run,dwce_ds_run)
                        
                        #calculate time derivative of the parallel wave number
                        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h 
                    else:
                        Ezw_run,beta_run,wR_run,wL_run=0,0,0,0
                        pwR_run,pwL_run=0,0
                        kz_run,w_wave=0,0
                        EwR_run,EwL_run=0,0
                        kapp_run=0
                        
                    #calculate the Runge-Kutta coeeficients for each differential equation
                    k2=wpi.dzdt(pparrk_su[k,i]+0.5*h*m1,gamma_run,particle_mass)
                    l2=wpi.dlamdadt(pparrk_su[k,i]+0.5*h*m1,lamdark_su[k,i]+0.5*h*l1,gamma_run,particle_mass,L_shell)
                    m2=wpi.dppardt(pperrk_su[k,i]+0.5*h*n1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    n2=wpi.dpperdt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    o2=wpi.detadt(pparrk_su[k,i]+0.5*h*m1,m_res,wcyc_run,gamma_run,kz_run,particle_mass,w_wave)
                    p2=wpi.dgammadt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run)
                    q2=wpi.daeqdt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,aeqrk_su[k,i]+0.5*h*q1,Ezw_run,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run)
                    r2=wpi.dEkdt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
                    s2=wpi.dalphadt(pperrk_su[k,i]+0.5*h*n1,pparrk_su[k,i]+0.5*h*m1,etark_su[k,i]+0.5*h*o1,Ezw_run,m_res,particle_charge,pwR_run,pwL_run,beta_run,wR_run,wL_run)
            
            #    ######################################################################################################
            #    #Third step of Runge Kutta
            #    ######################################################################################################
            
                    #calculate the magnetic field strength
                    B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i]+0.5*h*l2)
                    #calculate the Lorentz factor
                    gamma_run=((E_kin_su[k,i]+0.5*h*r2)/(particle_mass*env.const.c_light*env.const.c_light))+1
                    
                    #calculate the hydrogen cyclotron frequency
                    wcyc_run=env.omega_cyclotron(B_run,particle_charge,particle_mass)
            #         alphark_run=np.arctan(pperrk_su[k,i]+0.5*h*n2/pparrk_su[k,i]+0.5*h*m2)
                
                    #calculate the gradient of the hydrogen cyclotron frequency
                    dwce_ds_run=env.dwc_ds(wcyc_run,lamdark_su[k,i]+0.5*h*l2,L_shell)
                    
                    #calculate the gradient of the magnetic field strength
                    dBdz_run=env.dB_ds(B_run,lamdark_su[k,i]+0.5*h*l2,L_shell)
                    
                    if abs(np.rad2deg(lamdark_su[k,i]))<=lat_crossover_deg:
                        w_wave=wave_frequency
                    
                        #calculate the Stix parameters
                        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
                        
                        #calculate the refractive index and the wave numbers
                        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)
                        ref_left_run,ref_right_run,eta_sq_plus,eta_sq_minus,kapp_run,kz_run,kx_run=refr_index_full_left.refr_index_full_left(theta0,w_wave,S_run,D_run,P_run,R_run,L_run)
            
                        #calculate the amplitudes of the wave magnetic and electric field components
                        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(ref_left_run,P_run,D_run,S_run,Byw0_s,theta0)
                        #kx_run=0
                        
                        #calculate the wpi parameters
                        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i]+0.5*h*n2,kx_run,particle_charge,particle_mass,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)
                       
                        #calculate time derivative of the fyrofrequency
                        dwcdt_run=wpi.dwcdt(pparrk_su[k,i]+0.5*h*m2,particle_mass,gamma_run,dwce_ds_run)
                        
                        #calculate time derivative of the parallel wave number
                        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h
                    else:
                        Ezw_run,beta_run,wR_run,wL_run=0,0,0,0
                        pwR_run,pwL_run=0,0
                        kz_run,w_wave=0,0
                        EwR_run,EwL_run=0,0
                        kapp_run=0
                        
                    
                    #calculate the Runge-Kutta coeeficients for each differential equation
                    k3=wpi.dzdt(pparrk_su[k,i]+0.5*h*m2,gamma_run,particle_mass)
                    l3=wpi.dlamdadt(pparrk_su[k,i]+0.5*h*m2,lamdark_su[k,i]+0.5*h*l2,gamma_run,particle_mass,L_shell)
                    m3=wpi.dppardt(pperrk_su[k,i]+0.5*h*n2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    n3=wpi.dpperdt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    o3=wpi.detadt(pparrk_su[k,i]+0.5*h*m2,m_res,wcyc_run,gamma_run,kz_run,particle_mass,w_wave)
                    p3=wpi.dgammadt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run)
                    q3=wpi.daeqdt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,aeqrk_su[k,i]+0.5*h*q2,Ezw_run,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run)
                    r3=wpi.dEkdt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
                    s3=wpi.dalphadt(pperrk_su[k,i]+0.5*h*n2,pparrk_su[k,i]+0.5*h*m2,etark_su[k,i]+0.5*h*o2,Ezw_run,m_res,particle_charge,pwR_run,pwL_run,beta_run,wR_run,wL_run)
            
            #    ######################################################################################################
            #    #Fourth step of Runge Kutta
            #    ######################################################################################################
                     
            
                    #calculate the magnetic field strength
                    B_run=env.Bmag_dipole(L_shell,lamdark_su[k,i]+h*l3)
                    #calculate the Lorentz factor
                    gamma_run=((E_kin_su[k,i]+h*r3)/(particle_mass*env.const.c_light*env.const.c_light))+1
                    
                    #calculate the hydrogen cyclotron frequency
                    wcyc_run=env.omega_cyclotron(B_run,particle_charge,particle_mass)
            #         alphark_run=np.arctan(pperrk_su[k,i]+0.5*h*n3/pparrk_su[k,i]+0.5*h*m3)
                
                    #calculate the gradient of the hydrogen cyclotron frequency
                    dwce_ds_run=env.dwc_ds(wcyc_run,lamdark_su[k,i]+h*l3,L_shell)
                    
                    #calculate the gradient of the magnetic field strength
                    dBdz_run=env.dB_ds(B_run,lamdark_su[k,i]+h*l3,L_shell)
                    if abs(np.rad2deg(lamdark_su[k,i]))<=lat_crossover_deg:
                        w_wave=wave_frequency
                    
                        #calculate the Stix parameters
                        S_run,D_run,P_run,R_run,L_run=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, B_run)
                        
                        #calculate the refractive index and the wave numbers
                        eta_sq_plus,eta_sq_minus,mu_run,kapp_run,kz_run,kx_run=wave.refr_index_full(theta0,w_wave,S_run,P_run,R_run,L_run)
                        ref_left_run,ref_right_run,eta_sq_plus,eta_sq_minus,kapp_run,kz_run,kx_run=refr_index_full_left.refr_index_full_left(theta0,w_wave,S_run,D_run,P_run,R_run,L_run)
            
                        #calculate the amplitudes of the wave magnetic and electric field components
                        Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(ref_left_run,P_run,D_run,S_run,Byw0_s,theta0)
                        #kx_run=0
                        
                        #calculate the wpi parameters
                        beta_run,BwR_run,BwL_run,EwR_run,EwL_run,pwR_run,pwL_run,wR_run,wL_run=wpi.wpi_params(pperrk_su[k,i]+h*n3,kx_run,particle_charge,particle_mass,B_run,Exw_run,Eyw_run,Bxw_run,Byw_run,gamma_run)
                
                        #calculate time derivative of the fyrofrequency
                        dwcdt_run=wpi.dwcdt(pparrk_su[k,i]+h*m3,particle_mass,gamma_run,dwce_ds_run)
                        
                        #calculate time derivative of the parallel wave number
                        dkpardt_run=(kappa_out_su[k,i]-kappa_out_su[k,i-1])/h
                    else:
                        Ezw_run,beta_run,wR_run,wL_run=0,0,0,0
                        pwR_run,pwL_run=0,0
                        kz_run,w_wave=0,0
                        EwR_run,EwL_run=0,0
                        kapp_run=0
                        
                    
                    #calculate the Runge-Kutta coeeficients for each differential equation
                    k4=wpi.dzdt(pparrk_su[k,i]+h*m3,gamma_run,particle_mass)
                    l4=wpi.dlamdadt(pparrk_su[k,i]+h*m3,lamdark_su[k,i]+h*l3,gamma_run,particle_mass,L_shell)
                    m4=wpi.dppardt(pperrk_su[k,i]+h*n3,etark_su[k,i]+h*o3,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    n4=wpi.dpperdt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run,B_run,dBdz_run)
                    o4=wpi.detadt(pparrk_su[k,i]+h*m3,m_res,wcyc_run,gamma_run,kz_run,particle_mass,w_wave)
                    p4=wpi.dgammadt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run)
                    q4=wpi.daeqdt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,aeqrk_su[k,i]+h*q3,Ezw_run,gamma_run,m_res,particle_charge,particle_mass,pwR_run,pwL_run,beta_run,wR_run,wL_run)
                    r4=wpi.dEkdt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,gamma_run,m_res,particle_charge,particle_mass,Ezw_run,beta_run,EwR_run,EwL_run,w_wave,kapp_run)
                    s4=wpi.dalphadt(pperrk_su[k,i]+h*n3,pparrk_su[k,i]+h*m3,etark_su[k,i]+h*o3,Ezw_run,m_res,particle_charge,pwR_run,pwL_run,beta_run,wR_run,wL_run)
                    
                  
            #    ######################################################################################################
            #    Calculate parameters at final step
            #    ######################################################################################################
                    
                    #zetark_su[k,i+1]=zetark_su[k,i]+(h/6)*(k1+2*k2+2*k3+k4)
                    pparrk_su[k,i+1]=pparrk_su[k,i]+(h/6)*(m1+2*m2+2*m3+m4)
                    pperrk_su[k,i+1]=pperrk_su[k,i]+(h/6)*(n1+2*n2+2*n3+n4)
                    etark_su[k,i+1]=(etark_su[k,i]+(h/6)*(o1+2*o2+2*o3+o4))
                    lamdark_su[k,i+1]=lamdark_su[k,i]+(h/6)*(l1+2*l2+2*l3+l4)
                    alphark_su[k,i+1]=alphark_su[k,i]+(h/6)*(s1+2*s2+2*s3+s4)
                    aeqrk_su[k,i+1]=aeqrk_su[k,i]+(h/6)*(q1+2*q2+2*q3+q4)
                    E_kin_su[k,i+1]=E_kin_su[k,i]+(h/6)*(r1+2*r2+2*r3+r4)
                    #gamma_out_su[k,i+1]=gamma_out_su[k,i]+(h/6)*(p1+2*p2+2*p3+p4)
                    detadt_out[k,i+1]=(1/6)*(o1+2*o2+2*o3+o4)
                    #u_par_out_su[k,i+1]=pparrk_su[k,i+1]/(gamma_run*particle_mass)
                    #u_per_out_su[k,i+1]=pperrk_su[k,i+1]/(gamma_run*particle_mass)  
                    
            #    ######################################################################################################
            #    Non-linear investigation
            #    ######################################################################################################
                    if abs(np.rad2deg(lamdark_su[k,i]))<=lat_crossover_deg:
                        pass
                    else:
                        dkpardt_run=np.nan
                        
                    C0_run=wpi.nonlinear_C0(pparrk_su[k,i+1],kz_run,m_res,gamma_run,particle_charge,particle_mass,wcyc_run,Ezw_run)
                    C1p_run=wpi.nonlinear_C1p(pperrk_su[k,i+1],pparrk_su[k,i+1],kz_run,m_res,particle_charge,particle_mass,gamma_run,wR_run,EwR_run,wcyc_run)
                    C1m_run=wpi.nonlinear_C1m(pperrk_su[k,i+1],pparrk_su[k,i+1],kz_run,m_res,particle_charge,particle_mass,gamma_run,wL_run,EwL_run,wcyc_run)
                    thet_run,wtrsq_run= wpi.nonlinear_theta(C0_run,C1p_run,C1m_run,m_res,beta_run)
                    H_nonlinear_run=wpi.nonlinear_H(pperrk_su[k,i+1],pparrk_su[k,i+1],kz_run,gamma_run,m_res,particle_mass,wcyc_run,dkpardt_run,dwce_ds_run,0)
                    S_nonlinear_run=wpi.nonlinear_S(H_nonlinear_run,wtrsq_run) 
                    
                    
                    if abs(np.rad2deg(lamdark_su[k,i]))<=lat_crossover_deg:
                        pass
                    else:
                        C0_run,C1p_run,C1m_run=np.nan,np.nan,np.nan
                        thet_run,wtrsq_run=np.nan,np.nan
                        H_nonlinear_run,S_nonlinear_run=np.nan,np.nan
                        kx_run,kz_run=np.nan,np.nan
                        Exw_run,Eyw_run,Bxw_run=np.nan,np.nan,np.nan
                        Bxw_run,Byw_run,Bzw_run=np.nan,np.nan,np.nan
                        mu_run=np.nan
                    
                    # C0_out[k,i+1]=C0_run
                    # C1p_out[k,i+1]=C1p_run
                    # C1m_out[k,i+1]=C1m_run
                    # thet_out[k,i+1]=thet_run
                    # wtrsq_out[k,i+1]=wtrsq_run
                    # H_out[k,i+1]=H_nonlinear_run
                    # S_out[k,i+1]=S_nonlinear_run
                    
                    #output wave numbers
                    #kx_out[k,i+1]=kx_run
                    #kz_out[k,i+1]=kz_run
                    
                    # #outpout wave fields
                    # Exw_out_su[k,i+1]=Exw_run
                    # Eyw_out_su[k,i+1]=Eyw_run
                    # Ezw_out_su[k,i+1]=Ezw_run
                    # Bxw_out_su[k,i+1]=Bxw_run
                    # Byw_out_su[k,i+1]=Byw_run
                    # Bzw_out_su[k,i+1]=Bzw_run
            
                    #output first adiabatic inveriant
                    mu_out_su[k,i+1]=mu_run
                     
                    #calculate simulation time
                    #i=i+1
                    timerk_su[k,i+1]=timerk_su[k,i]+h
                    kappa_out_su[k,i+1]=kapp_run
                    B_out[k,i+1]=B_run
                    #if i==10801 and int(np.rad2deg(eta0[k]))==127:
                    #    sys.exit()
                #calculate delta alpha-eq and delta Ek
                fp=(~np.isnan(lamdark_su[k,:]))&(pparrk_su[k,:]<0)
                aLatitude_deg=np.rad2deg(lamdark_su[k,fp])

                aTime=timerk_su[k,fp]
                aAlpha_eq=aeqrk_su[k,fp]
                aEk=E_kin_su[k,fp]
                aPpar_mome=pparrk_su[k,fp]
                
                transmit_time[iPar,iEk,iAeq,k]=aTime[-1]-aTime[0]
                aeq_start[iPar,iEk,iAeq,k]=aAlpha_eq[0]
                aeq_end[iPar,iEk,iAeq,k]=aAlpha_eq[-1]
                ek_start[iPar,iEk,iAeq,k]=aEk[0]
                ek_end[iPar,iEk,iAeq,k]=aEk[-1]


                
                #delta_aeq[iPar,iEk,iAeq,k]=aAlpha_eq[-1]-aAlpha_eq[0]
                #delta_ek[iPar,iEk,iAeq,k]=aEk[-1]-aEk[0]
                if aLatitude_deg[-1]>0:
                    raise ValueError('the proton does not complete a bounce')
                if transmit_time[iPar,iEk,iAeq,k]<0.1:
                    raise ValueError('transmit time is too short')

            del(pparrk_su,pperrk_su,etark_su,lamdark_su,alphark_su,\
                aeqrk_su,E_kin_su,detadt_out)
            del(mu_out_su,timerk_su,kappa_out_su,B_out)
            
            
            
            



data={'transmit_time':transmit_time,'aeq_start':aeq_start,'aeq_end':aeq_end,\
      'ek_start':ek_start,'ek_end':ek_end,\
      'par_number':par_number,'emic_intensity':emic_intensity,'freq_norm':freq_norm,'wna':wna,\
          'Energy':Energy,'pitch_angle':pitch_angle,'eta0_deg':eta0_deg}
save_name='test_particle_diffusion_coef_678.mat'
sio.savemat(save_dir+save_name,data)

    
    


    
   
    
    
 
    
    

































