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

#import WPIT modules
import WPIT.Environment_mod as env
import WPIT.WaveProperties_mod as wave
import WPIT.WPI_mod.EMIC_ion_mod as wpi


#Ne=210   #number density
wna=[1e-3,1e-3,1e-3]
freq_norm=[0.7,0.8,0.9]
emic_intensity=[3*10**(-9),3*10**(-9),3*10**(-9)]

Energy=[10,10,10]  #initial energy keV
pitch_angle=[29,29,29]

density_ratio=15
L_shell=4 #L shell of simulation
#aeq0_deg=53
lat_cutoff_0=0.01
lat_cutoff_1=4
lat_cutoff_2=8




simu_time=30  #seconds
eta0_deg=np.linspace(0,360,49) #initial phases of electrons, 48 test electrons are used
eta0=np.deg2rad(eta0_deg) #convert initial phases to rad
m_res=1 #WPI resonance number (0=Landau resonance)
k_p_par=-1
particle_mass=env.const.mH
particle_charge=env.const.qi
ion_com=[77,20,3]


#sys.exit()
save_dir=r'E:\mywork\my_drafts\nonlinear_affected_by_EMIC_parameters\figs/'
save=True


LineWidth=0.8;LabelSize=15;TickWidth=1.6;TickLen=6
alpha_grid=0.4
#figure, axis lim
# xlim_s_lat=[0,30]
# ylim_s_lat=[-4,0]

xlim_eta_lat=[0,15]

xticks_eta_lat=range(0,16,2)

ylim_s_lat=[-4,0]


xlim_aeq_lat=[0,15]

ylim_aeq_lat_0=[0,50]
ylim_aeq_lat_1=[10,60]
ylim_aeq_lat_2=[10,80]



ylim_daeq_eta_0=[-20,10]
ylim_daeq_eta_1=[-20,30]
ylim_daeq_eta_2=[-20,40]


ylim_Ek_lat_0=[9,11]
ylim_Ek_lat_1=[9,12]
ylim_Ek_lat_2=[9,12]

ylim_dEk_eta_0=[-0.5,0.5]
ylim_dEk_eta_1=[-1,2]
ylim_dEk_eta_2=[-1,2]



#panel_label=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)']
panel_label=['(a)','(f)','(k)',\
             '(b)','(g)','(l)',\
                 '(c)','(h)','(m)',\
                     '(d)','(i)','(n)',\
                         '(e)','(j)','(o)']
plt.close('all')
num_row=5
num_column=3
myfig=plt.figure(figsize=(18,15))
hf=numlib.axes_create(row=num_row, column=num_column,top=0.04,bottom=0.04,left=0.05,right=0.06,gap_row=0.037,gap_column=0.06)


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




#sys.exit()
for iPar in range(len(wna)):
    
    theta0_nh_deg=wna[iPar]  # wave normal angle
    theta0_nh_rad=np.deg2rad(theta0_nh_deg) #convert wave normal angle to rad

    theta0_sh_deg=180-theta0_nh_deg  # wave normal angle
    theta0_sh_rad=np.deg2rad(theta0_sh_deg) #convert wave normal angle to rad
    
    Ekev0=Energy[iPar]
    aeq0_deg=pitch_angle[iPar]
    
    freq_ratio=freq_norm[iPar]
    wave_freq=freq_ratio*wcHe_eq

    w_wave=wave_freq
    wave_frequency=wave_freq
    
    Byw0_sim=emic_intensity[iPar]  # y-component of the wave magnetic field (3nT)

    

    
    
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
    for i in range(len(lat_dip)):        
        aLat_deg=lat_dip[i]
        aLat_rad=np.deg2rad(aLat_deg)
        aBfield=env.Bmag_dipole(L_shell,aLat_rad)
        S0,D0,P0,R0,L0=wave.stix_parameters(w_wave, ne_eq, nH_eq, nHe_eq, nO_eq, aBfield)
        eta_sq_plus0,eta_sq_minus0,mu0,kappa0,kappaz0,kappax0=wave.refr_index_full(theta0_nh_rad,w_wave,S0,P0,R0,L0)
        plr_plus[i]=-(eta_sq_plus0-S0)/D0
        plr_minus[i]=-(eta_sq_minus0-S0)/D0
        ref_left0,ref_right0,eta_sq_plus0,eta_sq_minus0,kappa0,kappaz0,kappax0=refr_index_full_left.refr_index_full_left(theta0_nh_rad,w_wave,S0,D0,P0,R0,L0)

        aRef_left[i]=ref_left0
        aRef_right[i]=ref_right0
        aD_stix[i]=D0
    fp_index=np.where(np.isnan(aRef_left))
    lat_cutoff_deg=lat_dip[fp_index[0][0]-1]  #lat_cutoff_deg is the last point where the waves can propagate
    bD_stix=np.array([aD_stix[i]*aD_stix[i+1] for i in range(len(aD_stix)-1)])
    fp_D=np.where(bD_stix<0)
    lat_crossover_deg=lat_dip[fp_D[0][0]]  #this is the latitude of crossover frequency
    #sys.exit()
    
    
            
    # 1. Define simulation parameters
    # Here we define all the initial parameters of the simulation in respect with the particle and the wave
        
    aeq0=np.deg2rad(aeq0_deg) #convert pitch angle to rad
    lamda0_mirr_rad=np.deg2rad(lamda0_deg) #convert latitude to rad    
    
    ### Integration parameters ##################################################################
    Tgyro=(2*np.pi)/wcH_eq #proton gyro period
    t=simu_time #simulation duration (s) 
    h=Tgyro/50 #simulation stepsize
    Nsteps=int(t/h) #number of simulation steps
    
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
    zetark_su=pperrk_su.copy()
    alphark_su=pperrk_su.copy()
    #alpha2rk_su=np.zeros((len(eta0),Nsteps+1))
    aeqrk_su=pperrk_su.copy()
    #aeq2rk_su=np.zeros((len(eta0),Nsteps+1))
    #aeq3rk_su=np.zeros((len(eta0),Nsteps+1))
    Exw_out_su=pperrk_su.copy()
    Eyw_out_su=pperrk_su.copy()
    Ezw_out_su=pperrk_su.copy()
    Bxw_out_su=pperrk_su.copy()
    Byw_out_su=pperrk_su.copy()
    Bzw_out_su=pperrk_su.copy()
    Bw_out_su=pperrk_su.copy()
    Ew_out_su=pperrk_su.copy()
    vresz_out_su=pperrk_su.copy()
    Eres_out_su=pperrk_su.copy()
    gammares_out_su=pperrk_su.copy()
    mu_adiabatic_out_su=pperrk_su.copy()
    mu_out_su=pperrk_su.copy()
    deta_dt_out_su=pperrk_su.copy()
    B_earth_out_su=pperrk_su.copy()
    S_stix_out_su=pperrk_su.copy()
    D_stix_out_su=pperrk_su.copy()
    P_stix_out_su=pperrk_su.copy()
    R_stix_out_su=pperrk_su.copy()
    L_stix_out_su=pperrk_su.copy()
    kappa_out_su=pperrk_su.copy()
    kx_out=pperrk_su.copy()
    kz_out=pperrk_su.copy()
    wh_out_su=pperrk_su.copy()
    dwce_ds_out_su=pperrk_su.copy()
    gamma_out_su=pperrk_su.copy()
    gamma2_out_su=pperrk_su.copy()
    
    C0_out=pperrk_su.copy()
    C1p_out=pperrk_su.copy()
    C1m_out=pperrk_su.copy()
    thet_out=pperrk_su.copy()
    wtrsq_out=pperrk_su.copy()
    dkpar_dtout=pperrk_su.copy()
    H_out=pperrk_su.copy()
    S_out=pperrk_su.copy()
    detadt_out=pperrk_su.copy()
    
    
    Phi_out_su=pperrk_su.copy()
    E_kin_su=pperrk_su.copy()
    E_kin_out=pperrk_su.copy()
    u_par_out_su=pperrk_su.copy()
    u_per_out_su=pperrk_su.copy()
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
        aeq0_rad=aeq0
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
        zetark_su[k,0]=0
        alphark_su[k,0]=alpha0
        #alpha2rk_su[k,0]=alpha0
        aeqrk_su[k,0]=aeq0_rad
        #aeq2rk_su[k,0]=aeq0
        #aeq3rk_su[k,0]=aeq0
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
        
        #start simulation below the latitude of crossover frequency
        while i<Nsteps:
    
    #    ######################################################################################################
    #    #First step of Runge Kutta
    #    ######################################################################################################
            #waves are present only to the northern hemisphere
            if lamdark_su[k,i]>=0:
                theta0=theta0_nh_rad
            else:
                theta0=theta0_sh_rad

                
            Byw0_s=Byw0_sim
            
            if (np.rad2deg(lamdark_su[k,i]))<0.01:  #run only for the northern hemisphere
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
            u_par_out_su[k,i+1]=pparrk_su[k,i+1]/(gamma_run*particle_mass)
            u_per_out_su[k,i+1]=pperrk_su[k,i+1]/(gamma_run*particle_mass)  
            
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
            
            C0_out[k,i+1]=C0_run
            C1p_out[k,i+1]=C1p_run
            C1m_out[k,i+1]=C1m_run
            thet_out[k,i+1]=thet_run
            wtrsq_out[k,i+1]=wtrsq_run
            H_out[k,i+1]=H_nonlinear_run
            S_out[k,i+1]=S_nonlinear_run
            
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
                
        
        
            
    # Plots of the results
    if iPar==0:
        lat_cutoff=lat_cutoff_0
    elif iPar==1:
        lat_cutoff=lat_cutoff_1
    else:
        lat_cutoff=lat_cutoff_2

        
    # Define colormaps for ploting
    
    fonts=15
    eta0_sim=eta0
    from matplotlib import cm
    
    norm = cm.colors.Normalize(vmin=eta0_sim.min(), vmax=eta0_sim.max())
    #cmap = cm.ScalarMappable(norm=norm, cmap=cm.nipy_spectral)
    cmap = cm.ScalarMappable(norm=norm, cmap=cm.jet)

    #cmap = cm.ScalarMappable(norm=norm, cmap=cm.hsv)

    cmap.set_array([])
    cmap2 = cm.ScalarMappable(norm=norm, cmap=cm.viridis)
    cmap2.set_array([])
    
    

    # figure 0, Wave-ion phase vs magnetic latitude
    axes_num=int(iPar%num_column)

    from matplotlib.axes._axes import _log as matplotlib_axes_logger
    matplotlib_axes_logger.setLevel('ERROR')
    import warnings
    warnings.simplefilter("ignore")
    
    last=1200
    #detadt_out
    for r in range(0,len(eta0_sim)-1):
        aDetadt=detadt_out[r,:]
        ResNum=0
        for iDe in range(len(aDetadt)-1):
            bDetadt=aDetadt[iDe]*aDetadt[iDe+1]
            if bDetadt<0: ResNum=ResNum+1
        if ResNum>=2: # the occurrence of phase trapping
            print('Energy: '+str(Ekev0)+' keV')
    
            print('eta0: '+str(np.rad2deg(eta0_sim[r])))
            print('Resonance times: '+str(ResNum))
            print('----------------')
            #for r in range(0,len(eta0_sim)):
            aLamda=np.rad2deg(lamdark_su[r,:-1])
            aEta=(etark_su[r,:-1]-np.pi)%(2*np.pi)
            eta_index=[]
            for i in range(len(aEta)-1):
                deta_eta=aEta[i]-aEta[i+1]
                if abs(deta_eta)>1:
                    eta_index.append(i)
            if len(eta_index)==0:
                cLamda=aLamda.copy()
                cEta=aEta/np.pi
                hf[axes_num].plot(cLamda,cEta,lw=LineWidth,c=cmap.to_rgba(eta0_sim[r]))
            else:
                if eta_index[0]>1: eta_index=[-1]+eta_index
                if eta_index[-1]<len(aEta): eta_index.append(len(aEta))
                for j in range(len(eta_index)-1):
                    cLamda=aLamda[eta_index[j]+1:eta_index[j+1]+1]
                    cEta=aEta[eta_index[j]+1:eta_index[j+1]+1]/np.pi
                    hf[axes_num].plot(cLamda,cEta,lw=LineWidth,c=cmap.to_rgba(eta0_sim[r]))
        else: #without phase trapping
            aLamda=np.rad2deg(lamdark_su[r,:-1])
            aEta=(etark_su[r,:-1]-np.pi)%(2*np.pi)
            eta_index=[]
            fp=aLamda>=lat_cutoff
            bLamda=aLamda[fp]
            bEta=aEta[fp]
            for i in range(len(bEta)-1):
                deta_eta=bEta[i]-bEta[i+1]
                if abs(deta_eta)>1:
                    eta_index.append(i)
            if len(eta_index)==0:
                cLamda=bLamda.copy()
                cEta=bEta/np.pi
                hf[axes_num].plot(cLamda,cEta,lw=LineWidth,c=cmap.to_rgba(eta0_sim[r]))
            else:
                if eta_index[0]>1: eta_index=[-1]+eta_index
                if eta_index[-1]<len(bEta): eta_index.append(len(bEta))
                for j in range(len(eta_index)-1):
                    cLamda=bLamda[eta_index[j]+1:eta_index[j+1]+1]
                    cEta=bEta[eta_index[j]+1:eta_index[j+1]+1]/np.pi
                    hf[axes_num].plot(cLamda,cEta,lw=LineWidth,c=cmap.to_rgba(eta0_sim[r]))
        
        
    
    hf[axes_num].grid(alpha=alpha_grid)
    hf[axes_num].set_xlim(xlim_eta_lat)
    hf[axes_num].set_xticks(xticks_eta_lat)

    hf[axes_num].set_ylim(0,2)
    hf[axes_num].set_yticks([0,0.5,1,1.5,2])
    hf[axes_num].set_yticklabels(['0',r'$\rm\pi/2$', r'$\rm\pi$',r'$\rm3\pi/2$',r'$\rm2\pi$'])
    hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
    hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)
    
    hf[axes_num].text(x=0.5,y=-0.19,s=r'Latitude $\rm\lambda\ (deg)$',fontsize=LabelSize,ha='center',transform=hf[axes_num].transAxes)
    hf[axes_num].set_ylabel(r'$\eta$'+' (rad)',fontsize=LabelSize)
    #hf[axes_num].text(x=0.5,y=1.05,s=r'$\rm\omega=$'+str(freq_ratio[iPar])+r'$\rm\Omega_{H^+}$'+r', $\rmE_k$='+str(Ekev0)+r' keV, $\rm\alpha_{eq}$='+str(aeq0_deg)+r'$\rm^\degree$',fontsize=LabelSize+1,ha='center',transform=hf[axes_num].transAxes)
    #hf[axes_num].text(x=0.5,y=1.04,s=r'$\rm\omega=$'+str(freq_ratio[iFreq]),fontsize=LabelSize+1,ha='center',transform=hf[axes_num].transAxes)
    hf[axes_num].text(x=0.5,y=1.05,s=r'$\rm\omega$='+str(freq_ratio)+r'$\rm\Omega_{He\_eq}$',fontsize=LabelSize+2,ha='center',transform=hf[axes_num].transAxes)


    if axes_num==0:
        
        bar_axes=myfig.add_axes([0.12,0.1,0.1,0.1])
        cbar=myfig.colorbar(cmap, ticks=[0,np.pi/2,np.pi,3/2*np.pi,2*np.pi],cax=bar_axes)
        cbar.ax.set_yticklabels(['0',r'$\rm\pi/2$', r'$\rm\pi$',r'$\rm3\pi/2$',r'$\rm2\pi$'],fontsize=LabelSize)  # horizontal colorbar
        
        
        cbar.set_label(r'$\eta_0$ (rad)', rotation=90,labelpad=-5,fontsize=LabelSize)
        pos=hf[axes_num].get_position()
        
        bar_position=[0.947,pos.y0,0.01,(pos.y1-pos.y0)]
        cbar.ax.set_position(bar_position)
    #plot linear resonance latitude
    #plot linear resonance latitude
    bResCon=[aResCon[i]*aResCon[i+1] for i in range(len(aResCon)-1)]
    bResCon.append(aResCon[-1]**2)
    fp=[i for i in range(len(bResCon)) if bResCon[i]<0]
    if len(fp)>0:
        bResLat=lamda_arr[fp[0]]
        aYlim=hf[axes_num].get_ylim()
        y_LinRes=np.linspace(aYlim[0],aYlim[1],100)
        x_LinRes=np.ones(len(y_LinRes))*bResLat
        hf[axes_num].plot(x_LinRes,y_LinRes,color='k',lw=2,ls='--')
        #hf1[axes_num].set_ylim(aYlim)
    
    
    

    
    
    
    
    
    
    
    
    
    
    #Equatorial pitch angle vs magnetic latitude
    axes_num=int(iPar%num_column)+num_column
    
    for r in range(0,len(eta0_sim)-1):
        fp=lamdark_su[r,:]>0
        hf[axes_num].plot(np.rad2deg(lamdark_su[r,fp]),np.rad2deg(aeqrk_su[r,fp]),lw=LineWidth,c=cmap.to_rgba(eta0_sim[r]))
        
    hf[axes_num].grid(alpha=alpha_grid)
    hf[axes_num].set_xlim(xlim_aeq_lat)
    if iPar==0:
        ylim_aeq_lat=ylim_aeq_lat_0
    elif iPar==1:
        ylim_aeq_lat=ylim_aeq_lat_1
    else:
        ylim_aeq_lat=ylim_aeq_lat_2
    
    hf[axes_num].set_ylim(ylim_aeq_lat)

    
    hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
    hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)

    hf[axes_num].text(x=0.5,y=-0.19,s=r'Latitude $\rm\lambda\ (deg)$',fontsize=LabelSize,ha='center',transform=hf[axes_num].transAxes)
    
    hf[axes_num].text(x=-0.15,y=0.5,s=r'$\rm\alpha_{eq}$ (deg)',fontsize=LabelSize,ha='center',va='center',rotation='90',transform=hf[axes_num].transAxes)

    #hf2[axes_num].text(x=0.5,y=1.04,s=r'$\rmN_e=$'+str(Ne)+r'$\rm\ cm^{-3}$'+', Energy='+str(Ekev0)+' keV',fontsize=LabelSize+1,ha='center',transform=hf1[axes_num].transAxes)
    if axes_num==3:
        
        bar_axes=myfig.add_axes([0.11,0.1,0.1,0.1])
        cbar=myfig.colorbar(cmap, ticks=[0,np.pi/2,np.pi,3/2*np.pi,2*np.pi],cax=bar_axes)
        cbar.ax.set_yticklabels(['0',r'$\rm\pi/2$', r'$\rm\pi$',r'$\rm3\pi/2$',r'$\rm2\pi$'],fontsize=LabelSize)  # horizontal colorbar
        
        cbar.set_label(r'$\eta_0$ (rad)', rotation=90,labelpad=-5,fontsize=LabelSize)
        pos=hf[axes_num].get_position()
        
        bar_position=[0.947,pos.y0,0.01,(pos.y1-pos.y0)]
        cbar.ax.set_position(bar_position)
    
    #plot linear resonance latitude
    if 'bResLat' in locals():
        aYlim=hf[axes_num].get_ylim()
        y_LinRes=np.linspace(aYlim[0],aYlim[1],100)
        x_LinRes=np.ones(len(y_LinRes))*bResLat
        hf[axes_num].plot(x_LinRes,y_LinRes,color='k',lw=2,ls='--')
        hf[axes_num].set_ylim(aYlim)
    
    
    #Calcualte total energy and pitch angle change
    
    dalpha=[]
    dEkin=[]
    
    from matplotlib import cm
    
    aeq = np.ma.array(aeqrk_su, mask=np.isnan(aeqrk_su))
    Ekin = np.ma.array(E_kin_su, mask=np.isnan(E_kin_su))
    for r in range(0,len(eta0_sim)-1):
        ls=np.max(np.nonzero(aeq[r,:]))
        els=np.max(np.nonzero(Ekin[r,:]))
        dalphaf=np.rad2deg(aeq[r,ls])-aeq0_deg
        dEkinf=(Ekin[r,ls]/1.602176487E-16)-Ekev0
        dalpha.append(dalphaf)
        dEkin.append(dEkinf)
        
    #plot Delta alpha_eq vs initial wave-electorn phase
    axes_num=int(iPar%num_column)+num_column*3
    
    colors=eta0_sim[:]
    
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),dalpha,color='k',lw=LineWidth)
    hf[axes_num].scatter(np.rad2deg(eta0_sim[:-1]),dalpha,marker='o',lw=2,facecolors='none', edgecolors=cmap.to_rgba(eta0_sim[:-1]),s=35)
    #plot average value
    dalpha_avg=np.mean(dalpha)
    avg_value=np.ones(len(dalpha))*dalpha_avg
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),avg_value,color='k',ls=':',lw=LineWidth+1)
    #plot standard deviation
    dalpha_std=np.std(dalpha)
    print('--------------')
    print('Energy: '+str(Ekev0)+' keV')
    print('Standard deviation of dalpha: '+str(dalpha_std))
    print('--------------')
    
    avg_plus=np.ones(len(dalpha))*(dalpha_avg+dalpha_std)
    avg_minus=np.ones(len(dalpha))*(dalpha_avg-dalpha_std)

    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),avg_plus,color='k',ls='--',lw=LineWidth+1)
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),avg_minus,color='k',ls='--',lw=LineWidth+1)

    hf[axes_num].grid(alpha=alpha_grid)

    hf[axes_num].set_xlim(0,360)
    if iPar==0:
        ylim_daeq_eta=ylim_daeq_eta_0
    elif iPar==1:
        ylim_daeq_eta=ylim_daeq_eta_1
    else:
        ylim_daeq_eta=ylim_daeq_eta_2
        
    hf[axes_num].set_ylim(ylim_daeq_eta)
    #hf[axes_num].set_ylim(ylim_daeq_eta)
    
    hf[axes_num].set_xticks(np.array([0,np.pi/2,np.pi,3*np.pi/2,2*np.pi])*180/np.pi)
    hf[axes_num].set_xticklabels(['0',r'$\rm\pi/2$', r'$\rm\pi$',r'$\rm3\pi/2$',r'$\rm2\pi$'])
    hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
    hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)

    hf[axes_num].text(x=0.5,y=-0.19,s=r'$\eta_0$ (deg)',fontsize=LabelSize,ha='center',transform=hf[axes_num].transAxes)

    #hf[axes_num].set_xlabel(r'$\eta_0$ (deg)',fontsize=LabelSize)
    hf[axes_num].text(x=-0.15,y=0.5,s=r'$\Delta \alpha_{\rmeq}$ (deg)',fontsize=LabelSize,ha='center',va='center',rotation='90',transform=hf[axes_num].transAxes)

    
    #plot Ek vs magnetic latitude
    axes_num=int(iPar%num_column)+num_column*2
    
    for r in range(0,len(eta0_sim)-1):
        fp=lamdark_su[r,:]>0
        hf[axes_num].plot(np.rad2deg(lamdark_su[r,fp]),E_kin_su[r,fp]/1.602176487E-16,c=cmap.to_rgba(eta0_sim[r]),lw=LineWidth)
        
    hf[axes_num].grid(alpha=alpha_grid)
    hf[axes_num].set_xlim(xlim_aeq_lat)
    if iPar==0:
        ylim_Ek_lat=ylim_Ek_lat_0
    elif iPar==1:
        ylim_Ek_lat=ylim_Ek_lat_1
    else:
        ylim_Ek_lat=ylim_Ek_lat_2
        
    hf[axes_num].set_ylim(ylim_Ek_lat)
    #hf[axes_num].set_ylim(ylim_Ek_lat)
    hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
    hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)
    
    hf[axes_num].text(x=0.5,y=-0.19,s=r'Latitude $\rm\lambda\ (deg)$',fontsize=LabelSize,ha='center',transform=hf[axes_num].transAxes)
    
    hf[axes_num].text(x=-0.15,y=0.5,s=r'$\rmE_k$ (keV)',fontsize=LabelSize,ha='center',va='center',rotation='90',transform=hf[axes_num].transAxes)
    
    #hf2[axes_num].text(x=0.5,y=1.04,s=r'$\rmN_e=$'+str(Ne)+r'$\rm\ cm^{-3}$'+', Energy='+str(Ekev0)+' keV',fontsize=LabelSize+1,ha='center',transform=hf1[axes_num].transAxes)
    #hf[axes_num].text(x=0.5,y=1.04,s=r'$\rm\omega_{pe}=$'+str(density_ratio)+r'$\rm\Omega_e$'+', Energy='+str(Ekev0)+' keV',fontsize=LabelSize+1,ha='center',transform=hf[axes_num].transAxes)
    if axes_num==6:
        
        bar_axes=myfig.add_axes([0.11,0.1,0.1,0.1])
        cbar=myfig.colorbar(cmap, ticks=[0,np.pi/2,np.pi,3/2*np.pi,2*np.pi],cax=bar_axes)
        cbar.ax.set_yticklabels(['0',r'$\rm\pi/2$', r'$\rm\pi$',r'$\rm3\pi/2$',r'$\rm2\pi$'],fontsize=LabelSize)  # horizontal colorbar
        
        cbar.set_label(r'$\eta_0$ (rad)', rotation=90,labelpad=-5,fontsize=LabelSize)
        pos=hf[axes_num].get_position()
        
        bar_position=[0.947,pos.y0,0.01,(pos.y1-pos.y0)]
        cbar.ax.set_position(bar_position)
    
    #plot linear resonance latitude
    if 'bResLat' in locals():
        aYlim=hf[axes_num].get_ylim()
        y_LinRes=np.linspace(aYlim[0],aYlim[1],100)
        x_LinRes=np.ones(len(y_LinRes))*bResLat
        hf[axes_num].plot(x_LinRes,y_LinRes,color='k',lw=2,ls='--')
        hf[axes_num].set_ylim(aYlim)
        
    #plot delta-Ek  versus latitude
    axes_num=int(iPar%num_column)+num_column*4
    
    colors=eta0_sim[:]
    
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),dEkin,color='k',lw=LineWidth)
    hf[axes_num].scatter(np.rad2deg(eta0_sim[:-1]),dEkin,marker='o',lw=2,facecolors='none', edgecolors=cmap.to_rgba(eta0_sim[:-1]),s=35)
    #plot average value
    dEk_avg=np.mean(dEkin)
    avg_value=np.ones(len(dEkin))*dEk_avg
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),avg_value,color='k',ls=':',lw=LineWidth+1)
    #plot standard deviation
    dEk_std=np.std(dEkin)
    print('--------------')
    print('Energy: '+str(Ekev0)+' keV')
    print('Standard deviation of dEkin: '+str(dEk_std))
    print('--------------')
    
    avg_plus=np.ones(len(dEkin))*(dEk_avg+dEk_std)
    avg_minus=np.ones(len(dEkin))*(dEk_avg-dEk_std)
    
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),avg_plus,color='k',ls='--',lw=LineWidth+1)
    hf[axes_num].plot(np.rad2deg(eta0_sim[:-1]),avg_minus,color='k',ls='--',lw=LineWidth+1)
    
    hf[axes_num].grid(alpha=alpha_grid)
    
    hf[axes_num].set_xlim(0,360)
    if iPar==0:
        ylim_dEk_eta=ylim_dEk_eta_0
    elif iPar==1:
        ylim_dEk_eta=ylim_dEk_eta_1
    else:
        ylim_dEk_eta=ylim_dEk_eta_2

    hf[axes_num].set_ylim(ylim_dEk_eta)
    
    hf[axes_num].set_xticks(np.array([0,np.pi/2,np.pi,3*np.pi/2,2*np.pi])*180/np.pi)
    hf[axes_num].set_xticklabels(['0',r'$\rm\pi/2$', r'$\rm\pi$',r'$\rm3\pi/2$',r'$\rm2\pi$'])
    hf[axes_num].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    hf[axes_num].tick_params(which='major',labelsize=LabelSize-1,width=TickWidth,length=TickLen,tickdir='in',right=True,top=True)
    hf[axes_num].tick_params(which='minor',width=TickWidth*0.6,length=TickLen*0.5,tickdir='in',right=True,top=True)
    
    hf[axes_num].text(x=0.5,y=-0.19,s=r'$\eta_0$ (deg)',fontsize=LabelSize,ha='center',transform=hf[axes_num].transAxes)
    
    #hf[axes_num].set_xlabel(r'$\eta_0$ (deg)',fontsize=LabelSize)
    hf[axes_num].text(x=-0.15,y=0.5,s=r'$\rm\Delta E_{\rmk}$ (keV)',fontsize=LabelSize,ha='center',va='center',rotation='90',transform=hf[axes_num].transAxes)
    


    #sys.exit()

  
axes_num=0
hf[axes_num].annotate(text='phase bunching',color='k',fontsize=LabelSize+1,xy=(0.22,0.76),\
    xycoords='axes fraction',xytext=(0.28,0.9),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))
  
    
axes_num=1
hf[axes_num].annotate(text='phase bunching',color='k',fontsize=LabelSize+1,xy=(0.57,0.78),\
    xycoords='axes fraction',xytext=(0.15,0.9),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))

axes_num=1
hf[axes_num].annotate(text='phase trapping',color='k',fontsize=LabelSize+1,xy=(0.15,0.48),\
    xycoords='axes fraction',xytext=(0.06,0.13),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))
   
axes_num=2
hf[axes_num].annotate(text='phase bunching',color='k',fontsize=LabelSize+1,xy=(0.77,0.95),\
    xycoords='axes fraction',xytext=(0.17,0.8),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))

axes_num=2
hf[axes_num].annotate(text='phase trapping',color='k',fontsize=LabelSize+1,xy=(0.27,0.36),\
    xycoords='axes fraction',xytext=(0.23,0.12),va='center',textcoords='axes fraction',arrowprops=dict(color='k',headwidth=9,width=1))
    

#set panel label for figure 1
#panel_num=0
for iP in range(len(hf)):
    hf[iP].text(x=0.01,y=0.89,s=panel_label[iP],color='k',weight='bold',bbox=dict(facecolor='w',edgecolor='None',alpha=0.35),fontsize=LabelSize+3,transform=hf[iP].transAxes)


  

fig_name='trajectory_EMIC_freq'
if save is True:
    myfig.savefig(save_dir+fig_name+'.png')
    myfig.savefig(save_dir+fig_name+'.tiff')


    
    


    
   
    
    
 
    
    

































